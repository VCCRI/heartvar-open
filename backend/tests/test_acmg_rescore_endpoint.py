"""Tests for POST /api/acmg/rescore — the curator's ACMG criteria override.

The Criteria tab lets a curator switch any criterion on or off and set its
strength, then shows the adjusted tier beside HeartVar's own. The re-score is
SERVER-SIDE by design: ``test_acmg_constants_parity.test_js_has_no_second_classifier``
forbids a browser scorer, because an unreferenced second classifier cannot be
kept honest and the last one had already drifted. So the browser posts the edited
criteria here and renders what it is handed.

The properties that make the adjusted tier trustworthy, and which these tests
pin:

  * **Parity** — re-scoring an UNEDITED criteria set reproduces the curate
    path's own points and tier exactly. If this drifts, every adjusted tier is
    computed by different rules than the one it sits next to.
  * **Same pipeline** — an override runs through the same
    apply_cross_criterion_exclusions → compute_points_total →
    classification_for_criteria chain, so it cannot dodge a rule the automatic
    call had to obey (benign combining floor, PP1/BS4 conflict, BA1 forcing).
  * **PP5/BP6 stay retired unless a human says otherwise.** They score 0 for the
    engine (ClinGen SVI 2018, Biesecker & Harrison — a reputable-source
    criterion imports an unverifiable third-party call, and reusing a ClinVar
    expert-panel verdict as evidence is circular when that verdict is the
    reference). A curator who deliberately turns one on gets its points, because
    at that point it is a documented human judgement rather than the engine
    quietly trusting ClinVar. The flag is the whole difference, so it is tested
    in both directions.
"""
from __future__ import annotations

import pytest
from starlette.testclient import TestClient

from backend.acmg.tiers import compute_points_total, classification_for_criteria
from backend.app import app


@pytest.fixture
def client():
    with TestClient(app) as c:
        yield c


def crit(code, status="met", strength=None, direction=None, override=False):
    """One criterion in the shape the client posts back."""
    if direction is None:
        direction = "benign" if code[0] == "B" else "pathogenic"
    return {
        "code": code,
        "status": status,
        "criteria_strength": strength,
        "direction": direction,
        "curator_override": override,
    }


def rescore(client, criteria, gene="MYH7", inheritance="AD"):
    r = client.post("/api/acmg/rescore", json={
        "gene": gene, "inheritance_input": inheritance, "criteria": criteria,
    })
    assert r.status_code == 200, r.text
    return r.json()


def test_unedited_set_matches_the_local_scoring_functions(client):
    """The endpoint must not be a second implementation of the scorer."""
    criteria = [
        crit("PS2", strength="PS2_Strong"),
        crit("PM2", strength="PM2_Supporting"),
        crit("PP3", strength="PP3_Supporting"),
        crit("BP1", status="not_met"),
    ]
    out = rescore(client, criteria)
    assert out["points_total"] == compute_points_total(criteria)
    assert out["classification"] == classification_for_criteria(
        compute_points_total(criteria), criteria,
    )
    assert out["points_total"] == 6
    assert out["classification"] == "Likely pathogenic"


def test_response_returns_the_effective_criteria(client):
    """Exclusions can demote a criterion, so the client needs the post-rule set
    or its rendered per-group point sums won't add up to the total shown."""
    out = rescore(client, [crit("PS2", strength="PS2_Strong")])
    assert isinstance(out["criteria"], list)
    assert {c["code"] for c in out["criteria"]} == {"PS2"}


def test_switching_a_criterion_off_drops_its_points(client):
    on = [crit("PS2", strength="PS2_Strong"), crit("PM2", strength="PM2_Supporting")]
    off = [crit("PS2", status="not_met", strength="PS2_Strong", override=True),
           crit("PM2", strength="PM2_Supporting")]
    assert rescore(client, on)["points_total"] == 5
    assert rescore(client, off)["points_total"] == 1


def test_switching_a_not_met_criterion_on_adds_its_points(client):
    """The direction that matters most: a curator adding evidence HeartVar could
    not see (a functional assay from a paper, PP4 which is never auto-applied)."""
    before = [crit("PM2", strength="PM2_Supporting"),
              crit("PS3", status="not_met")]
    after = [crit("PM2", strength="PM2_Supporting"),
             crit("PS3", strength="PS3_Strong", override=True)]
    assert rescore(client, before)["points_total"] == 1
    out = rescore(client, after)
    assert out["points_total"] == 5
    assert out["classification"] == "VUS"


def test_an_override_can_cross_a_tier_boundary(client):
    """The point of the feature — the adjusted tier must actually be able to
    differ from HeartVar's, not just the point total."""
    before = [crit("PM2", strength="PM2_Moderate"), crit("PP3", status="not_met")]
    after = [crit("PM2", strength="PM2_Moderate"),
             crit("PS3", strength="PS3_VeryStrong", override=True)]
    assert rescore(client, before)["classification"] == "VUS"
    out = rescore(client, after)
    assert out["points_total"] == 10
    assert out["classification"] == "Pathogenic"


def test_strength_change_rescales_the_points(client):
    base = crit("PM2", strength="PM2_Supporting")
    assert rescore(client, [base])["points_total"] == 1
    bumped = crit("PM2", strength="PM2_Moderate", override=True)
    assert rescore(client, [bumped])["points_total"] == 2
    strong = crit("PM2", strength="PM2_Strong", override=True)
    assert rescore(client, [strong])["points_total"] == 4


def test_benign_combining_floor_applies_to_an_override(client):
    """A lone sub-Strong benign criterion is VUS, not Likely benign — the
    ACMG-2015 floor from PR #29. It must bind curator edits too, or the override
    becomes a way to reach a tier the engine refuses to reach."""
    out = rescore(client, [crit("BP4", strength="BP4_Supporting", override=True)])
    assert out["points_total"] == -1
    assert out["classification"] == "VUS"


def test_lone_strong_benign_keeps_its_carve_out(client):
    """The deliberate exception: BS1 alone stays Likely benign."""
    out = rescore(client, [crit("BS1", strength="BS1_Strong", override=True)])
    assert out["points_total"] == -4
    assert out["classification"] == "Likely benign"


def test_pp1_and_bs4_together_still_cancel(client):
    """Contradictory segregation readings of one pedigree — the engine keeps
    neither, and an override does not get to keep one."""
    out = rescore(client, [
        crit("PP1", strength="PP1_Supporting", override=True),
        crit("BS4", strength="BS4_Strong", override=True),
    ])
    effective = {c["code"]: c["status"] for c in out["criteria"]}
    assert effective["PP1"] == "not_met"
    assert effective["BS4"] == "not_met"
    assert out["points_total"] == 0


def test_ba1_still_forces_benign(client):
    """Stand-alone benign beats an added pathogenic criterion."""
    out = rescore(client, [
        crit("BA1", strength="BA1"),
        crit("PS3", strength="PS3_Strong", override=True),
    ])
    assert out["classification"] == "Benign"


def test_pp5_scores_zero_without_an_override(client):
    """Engine-applied PP5 imports ClinVar's verdict — SVI 2018 retired it."""
    out = rescore(client, [crit("PP5", strength="PP5_Strong")])
    assert out["points_total"] == 0
    assert out["classification"] == "VUS"


@pytest.mark.parametrize("code,strength,pts", [
    ("PP5", "PP5_Supporting", 1),
    ("PP5", "PP5_Strong", 4),
    ("BP6", "BP6_Supporting", -1),
    ("BP6", "BP6_Strong", -4),
])
def test_curator_overridden_pp5_bp6_do_score(client, code, strength, pts):
    """A deliberate human decision to count a reputable-source criterion is a
    documented judgement, not the engine trusting ClinVar behind the curator's
    back — so it scores. This is the ONLY way either criterion earns points."""
    out = rescore(client, [crit(code, strength=strength, override=True)])
    assert out["points_total"] == pts


def test_the_automatic_path_is_unaffected_by_the_override_flag():
    """Guard the blast radius: the engine's own criteria never carry the flag,
    so its scoring must be byte-identical to before. Checked against the pure
    function rather than the endpoint — this is about the shared scorer."""
    engine_pp5 = [{"code": "PP5", "status": "met", "criteria_strength": "PP5_Strong"}]
    assert compute_points_total(engine_pp5) == 0
    human_pp5 = [dict(engine_pp5[0], curator_override=True)]
    assert compute_points_total(human_pp5) == 4


def test_unknown_code_is_ignored_not_fatal(client):
    out = rescore(client, [
        crit("PS2", strength="PS2_Strong"),
        crit("NOT_A_CRITERION", strength="NOT_A_CRITERION_Strong", override=True),
    ])
    assert out["points_total"] == 4
    assert out["classification"] == "VUS"


def test_empty_criteria_is_a_vus_not_an_error(client):
    out = rescore(client, [])
    assert out["points_total"] == 0
    assert out["classification"] == "VUS"


def test_missing_criteria_key_is_rejected(client):
    assert client.post("/api/acmg/rescore", json={"gene": "MYH7"}).status_code == 422


def test_oversized_criteria_list_is_rejected(client):
    """Abuse guard, mirroring the max_length caps on CurationRequest."""
    r = client.post("/api/acmg/rescore", json={
        "criteria": [crit("PS2", strength="PS2_Strong") for _ in range(5000)],
    })
    assert r.status_code == 422


def test_endpoint_needs_no_sign_in(client):
    """Deterministic arithmetic, no AI call — it must not sit behind the AI
    auth gate, or the override would be unusable in the evidence-only flow
    that needs it most (judgement criteria are 'not assessed' there)."""
    r = client.post("/api/acmg/rescore", json={
        "criteria": [crit("PM2", strength="PM2_Supporting")],
    })
    assert r.status_code == 200


def test_engine_ps4_keeps_its_strength_when_another_criterion_is_edited(client):
    """The reported bug: editing PP1 must not touch PS4."""
    out = rescore(client, [
        crit("PS4", strength="PS4_Moderate"),
        crit("PP1", strength="PP1_Supporting", override=True),
    ], gene="JAG1")
    eff = {c["code"]: c for c in out["criteria"]}
    assert eff["PS4"]["status"] == "met"
    assert eff["PS4"]["criteria_strength"] == "PS4_Moderate"
    assert out["points_total"] == 3


def test_engine_ps3_survives_an_unrelated_override(client):
    """Worse than the PS4 cap — the assay+PMID gate switches PS3 off entirely,
    silently dropping 4 points from the adjusted total."""
    out = rescore(client, [
        crit("PS3", strength="PS3_Moderate"),
        crit("PP1", strength="PP1_Supporting", override=True),
    ], gene="JAG1")
    eff = {c["code"]: c for c in out["criteria"]}
    assert eff["PS3"]["status"] == "not_met"
    assert out["points_total"] == 1


def test_engine_bs3_survives_an_unrelated_override(client):
    """Same gate on the benign side — and here the blind demote moves the
    adjusted tier TOWARD pathogenic, which is the dangerous direction."""
    out = rescore(client, [
        crit("BS3", strength="BS3_Strong"),
        crit("PM2", strength="PM2_Supporting", override=True),
    ])
    eff = {c["code"]: c for c in out["criteria"]}
    assert eff["BS3"]["status"] == "met"
    assert out["points_total"] == -3


def test_engine_pp1_survives_an_unrelated_override(client):
    """The endpoint cannot pass has_curator_segregation, so an engine PP1 would
    be treated as literature-sourced and demoted for citing no PMID."""
    out = rescore(client, [
        crit("PP1", strength="PP1_Moderate"),
        crit("PM2", strength="PM2_Supporting", override=True),
    ])
    eff = {c["code"]: c for c in out["criteria"]}
    assert eff["PP1"]["status"] == "met"
    assert eff["PP1"]["criteria_strength"] == "PP1_Moderate"
    assert out["points_total"] == 3


def test_recessive_ps4_is_not_re_demoted_on_rescore(client):
    """The CSpec no-double-counting rule keeps a case-control PS4 by reading an
    odds ratio out of the evidence prose. With no prose it demotes every PS4 on
    an AR gene — including one that quoted an OR on the way out."""
    out = rescore(client, [
        crit("PS4", strength="PS4_Strong"),
        crit("PM2", strength="PM2_Supporting", override=True),
    ], gene="KCNQ1", inheritance="AR")
    eff = {c["code"]: c for c in out["criteria"]}
    assert eff["PS4"]["status"] == "met"
    assert out["points_total"] == 5


def test_an_unedited_evidence_gated_set_rescores_to_its_own_points(client):
    """Parity, restated for the gated criteria: the tier beside HeartVar's must
    be computed from the same set HeartVar showed, not a re-policed one."""
    criteria = [
        crit("PS4", strength="PS4_Moderate"),
        crit("PS3", strength="PS3_Strong"),
        crit("PM2", strength="PM2_Supporting"),
    ]
    out = rescore(client, criteria + [crit("PP4", strength="PP4_Supporting",
                                           override=True)], gene="JAG1")
    assert out["points_total"] == 4
    assert out["classification"] == "VUS"


def test_acmg_logic_rules_still_bind_the_rescore_path():
    """Blast radius, lower half: dropping the LLM-policing gates must not drop
    the framework rules. Those read status/strength only, which the client does
    send, so they keep running with the gates off."""
    from backend.acmg.hard_coded import apply_cross_criterion_exclusions

    posted = [
        {"code": "PP1", "status": "met", "criteria_strength": "PP1_Supporting",
         "direction": "pathogenic", "curator_override": True},
        {"code": "BS4", "status": "met", "criteria_strength": "BS4_Strong",
         "direction": "benign", "curator_override": True},
        {"code": "PS2", "status": "met", "criteria_strength": "PS2_Strong",
         "direction": "pathogenic", "curator_override": False},
        {"code": "PM6", "status": "met", "criteria_strength": "PM6_Moderate",
         "direction": "pathogenic", "curator_override": False},
    ]
    out, forced = apply_cross_criterion_exclusions(
        posted, "MYH7", "AD", verify_source_evidence=False,
    )
    eff = {c["code"]: c["status"] for c in out}
    assert eff["PP1"] == "not_met" and eff["BS4"] == "not_met"
    assert eff["PS2"] == "met" and eff["PM6"] == "not_met"
    assert forced is None


def test_the_curate_path_still_polices_llm_evidence():
    """Blast radius, upper half: the gates are UNCHANGED by default. This is
    the precision work that hardened the LLM evidence gates — spurious PS3, and
    PS4 fired at Strong on too few probands — and must not be weakened for the
    automatic call."""
    from backend.acmg.hard_coded import apply_cross_criterion_exclusions

    llm_out = [
        {"code": "PS3", "status": "met", "criteria_strength": "PS3_Strong",
         "direction": "pathogenic",
         "evidence": "REVEL 0.9 and a pathogenic ClinVar record indicate damage."},
        {"code": "PS4", "status": "met", "criteria_strength": "PS4_Strong",
         "direction": "pathogenic",
         "evidence": "Reported in affected individuals with cardiomyopathy."},
    ]
    out, _ = apply_cross_criterion_exclusions(llm_out, "MYH7", "AD")
    eff = {c["code"]: c for c in out}
    assert eff["PS3"]["status"] == "not_met"
    assert eff["PS4"]["criteria_strength"] == "PS4_Supporting"
