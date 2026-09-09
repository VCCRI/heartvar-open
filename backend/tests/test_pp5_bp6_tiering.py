"""Unit tests for PP5/BP6 evaluation.

PP5/BP6 are RETIRED from scoring per ClinGen SVI (Biesecker & Harrison, Genet
Med 2018): HeartVar still derives the tiered PP5/BP6 *entry* from the proband's
exact-variant ClinVar record — status, review-star strength LABEL
(3-4★→Strong, 2★→Moderate, 1★→Supporting, 0★/conflicting→not applied), the
PP5-vs-BA1/BS1 conflict guard and the post-merge BP6-vs-strong-pathogenic
guard are all preserved so the curator sees the ClinVar assertion as context —
but the entry contributes ZERO points and never moves the classification.
These tests pin every branch of that logic AND the zero-point invariant.

No pytest dependency — runnable with pytest or directly
(``python -m backend.tests.test_pp5_bp6_tiering``).
"""
from __future__ import annotations

from backend.app import (
    _clinvar_pp5_bp6_criteria,
    _select_clinvar_assertion_record,
    _CLINVAR_STAR_STRENGTH,
    apply_cross_criterion_exclusions,
    compute_points_total,
    classification_for,
    tier_without_clinvar_assertion,
    HARD_CODED_CRITERIA_CODES,
    _CRITERION_NAMES,
)
from backend.clients.clinvar import _stars_for, _classify_tier, _hgvs_match


RS_4STAR = "practice guideline"
RS_3STAR = "reviewed by expert panel"
RS_2STAR = "criteria provided, multiple submitters, no conflicts"
RS_1STAR = "criteria provided, single submitter"
RS_CONFLICT = "criteria provided, conflicting classifications"
RS_0STAR = "no assertion criteria provided"


def _rec(sig, review_status, n_sub=1, accession="VCV000012345"):
    """A ClinVar record dict shaped like fetch_clinvar's records[] entries.
    `stars` is computed exactly as the client does (via _stars_for)."""
    return {
        "accession": accession,
        "clinical_significance": sig,
        "review_status": review_status,
        "stars": _stars_for(review_status),
        "number_submitters": n_sub,
    }


def _cv(*records, ok=True):
    return {"ok": ok, "found": bool(records), "records": list(records)}


def _pp5_bp6(clinvar_ev, ba1_met=False, bs1_met=False):
    entries = _clinvar_pp5_bp6_criteria(
        clinvar_ev, ba1_met=ba1_met, bs1_met=bs1_met
    )
    by_code = {e["code"]: e for e in entries}
    return by_code["PP5"], by_code["BP6"]


def _crit(code, status="met", strength=None, evidence=""):
    if strength is None:
        strength = code if status == "met" else None
    return {
        "code": code,
        "name": _CRITERION_NAMES.get(code, code),
        "status": status,
        "direction": "benign" if code.startswith("B") else "pathogenic",
        "criteria_strength": strength,
        "evidence": evidence,
        "source": "test",
    }


def test_pp5_star_to_strength():
    for rs, expect in [
        (RS_4STAR, "PP5_Strong"),
        (RS_3STAR, "PP5_Strong"),
        (RS_2STAR, "PP5_Moderate"),
        (RS_1STAR, "PP5_Supporting"),
    ]:
        pp5, _ = _pp5_bp6(_cv(_rec("Pathogenic", rs)))
        assert pp5["status"] == "met", rs
        assert pp5["criteria_strength"] == expect, rs


def test_pp5_zero_star_not_applied():
    pp5, _ = _pp5_bp6(_cv(_rec("Pathogenic", RS_0STAR)))
    assert pp5["status"] == "not_met"
    assert pp5["criteria_strength"] is None


def test_bp6_star_to_strength():
    for rs, expect in [
        (RS_3STAR, "BP6_Strong"),
        (RS_2STAR, "BP6_Moderate"),
        (RS_1STAR, "BP6_Supporting"),
    ]:
        _, bp6 = _pp5_bp6(_cv(_rec("Benign", rs)))
        assert bp6["status"] == "met", rs
        assert bp6["criteria_strength"] == expect, rs
        assert compute_points_total([bp6]) == 0, rs


def test_bp6_zero_star_not_applied():
    _, bp6 = _pp5_bp6(_cv(_rec("Benign", RS_0STAR)))
    assert bp6["status"] == "not_met"


def test_direction_routing_pathogenic_side():
    for sig in ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"):
        pp5, bp6 = _pp5_bp6(_cv(_rec(sig, RS_3STAR)))
        assert pp5["status"] == "met", sig
        assert bp6["status"] == "not_met", sig


def test_direction_routing_benign_side():
    for sig in ("Benign", "Likely benign", "Benign/Likely benign"):
        pp5, bp6 = _pp5_bp6(_cv(_rec(sig, RS_3STAR)))
        assert bp6["status"] == "met", sig
        assert pp5["status"] == "not_met", sig


def test_vus_record_applies_neither():
    pp5, bp6 = _pp5_bp6(_cv(_rec("Uncertain significance", RS_3STAR)))
    assert pp5["status"] == "not_met"
    assert bp6["status"] == "not_met"


def test_conflicting_excluded_by_tier_not_by_stars():
    rec = _rec("Conflicting classifications of pathogenicity", RS_CONFLICT)
    assert rec["stars"] == 1
    assert _classify_tier(rec["clinical_significance"]) == "VUS"
    pp5, bp6 = _pp5_bp6(_cv(rec))
    assert pp5["status"] == "not_met"
    assert bp6["status"] == "not_met"


def test_no_record_not_applied():
    pp5, bp6 = _pp5_bp6(_cv())
    assert pp5["status"] == "not_met"
    assert bp6["status"] == "not_met"
    assert "No exact-variant ClinVar assertion" in pp5["evidence"]


def test_lookup_unavailable_not_applied():
    pp5, bp6 = _pp5_bp6(_cv(ok=False))
    assert pp5["status"] == "not_met"
    assert "unavailable" in pp5["evidence"].lower()


def test_pp5_withheld_under_ba1():
    pp5, _ = _pp5_bp6(_cv(_rec("Pathogenic", RS_3STAR)), ba1_met=True)
    assert pp5["status"] == "not_met"
    assert pp5["criteria_strength"] is None
    assert "withheld" in pp5["evidence"].lower()
    assert "BA1" in pp5["evidence"]
    assert compute_points_total([pp5]) == 0


def test_pp5_withheld_under_bs1_only():
    pp5, _ = _pp5_bp6(_cv(_rec("Pathogenic", RS_3STAR)), bs1_met=True)
    assert pp5["status"] == "not_met"
    assert "BS1" in pp5["evidence"]


def test_bp6_not_suppressed_when_ba1_agrees():
    _, bp6 = _pp5_bp6(_cv(_rec("Benign", RS_3STAR)), ba1_met=True)
    assert bp6["status"] == "met"
    assert bp6["criteria_strength"] == "BP6_Strong"


def test_max_star_record_selected_over_higher_submitter():
    high_sub_1star = _rec("Pathogenic", RS_1STAR, n_sub=50, accession="VCV1")
    low_sub_3star = _rec("Pathogenic", RS_3STAR, n_sub=2, accession="VCV2")
    pp5, _ = _pp5_bp6(_cv(high_sub_1star, low_sub_3star))
    assert pp5["criteria_strength"] == "PP5_Strong"
    assert "VCV2" in pp5["evidence"]


def test_select_assertion_record_filters_direction():
    cv = _cv(_rec("Benign", RS_3STAR), _rec("Pathogenic", RS_2STAR))
    rec, stars, tier = _select_clinvar_assertion_record(cv, {"P", "LP"})
    assert tier == "P" and stars == 2
    rec, stars, tier = _select_clinvar_assertion_record(cv, {"B", "LB"})
    assert tier == "B" and stars == 3


def _run_exclusions(*entries, gene=None):
    """`gene` matters for PS3 since 2026-09-01: the PS3 strength ceiling comes
    from the gene's VCEP spec, and a gene with no spec (or None) is held at
    PS3_Supporting per Brnich 2019. A test that needs a Strong PS3 must name a
    gene whose spec publishes a Strong row."""
    criteria, _forced = apply_cross_criterion_exclusions(list(entries), gene=gene)
    return {c["code"]: c for c in criteria}


def test_bp6_demoted_by_pvs1():
    out = _run_exclusions(
        _crit("BP6", strength="BP6_Moderate"),
        _crit("PVS1", strength="PVS1"),
    )
    assert out["BP6"]["status"] == "not_met"
    assert "PVS1" in out["BP6"]["evidence"]


def test_bp6_demoted_by_strong_pathogenic_evidence():
    """The guard fires on any of PVS1 / PS1 / PS3 / PS4 at full Strong.

    Uses PS1 rather than PS3. Since 2026-09-01 PS3 CANNOT reach Strong: the
    Brnich 2019 ladder tops out at Moderate (11 classified variant controls is
    its MODERATE minimum, and the Strong rung lives in supplementary OddsPath
    tables we do not implement). So a PS3-based version of this test would assert
    an unreachable state — see
    test_ps3_can_no_longer_reach_strong_so_it_cannot_demote_bp6."""
    out = _run_exclusions(
        _crit("BP6", strength="BP6_Strong"),
        _crit("PS1", strength="PS1"),
        gene="MYH7",
    )
    assert out["BP6"]["status"] == "not_met"


def test_bp6_demoted_by_strong_ps4():
    ps4 = _crit("PS4", strength="PS4_Strong",
                evidence="case-control odds ratio 9.4 (95% CI 4.0-22)")
    ps4["facts"] = {"or_lower_ci": 4.0}
    out = _run_exclusions(_crit("BP6", strength="BP6_Moderate"), ps4)
    assert out["BP6"]["status"] == "not_met"


def test_bp6_kept_when_ps4_or_is_prose_only():
    """The downstream half of the PS4 case-control fix. An odds ratio asserted
    only in prose no longer holds PS4 at Strong, so it no longer knocks out a
    contradicting BP6 either. Benign-direction, and the right answer: an
    unverifiable OR should not silently delete reputable-source benign
    evidence."""
    out = _run_exclusions(
        _crit("BP6", strength="BP6_Moderate"),
        _crit("PS4", strength="PS4_Strong",
              evidence="case-control odds ratio 9.4 (95% CI 4.0-22)"),
    )
    assert out["BP6"]["status"] == "met"


def test_bp6_kept_when_pathogenic_only_moderate():
    out = _run_exclusions(
        _crit("BP6", strength="BP6_Moderate"),
        _crit("PS3", strength="PS3_Moderate"),
        _crit("PVS1", strength="PVS1_Moderate"),
    )
    assert out["BP6"]["status"] == "met"


def test_bp6_kept_when_only_pm_level_pathogenic():
    out = _run_exclusions(
        _crit("BP6", strength="BP6_Strong"),
        _crit("PM2", strength="PM2_Supporting"),
        _crit("PP3", strength="PP3"),
    )
    assert out["BP6"]["status"] == "met"


def test_pp5_strong_does_not_lift_vus():
    pp5, _ = _pp5_bp6(_cv(_rec("Pathogenic", RS_3STAR)))
    assert pp5["status"] == "met"
    criteria = [_crit("PM1", strength="PM1_Moderate"), pp5]
    total = compute_points_total(criteria)
    assert total == 2
    assert classification_for(total) == "VUS"


def test_bp6_strong_does_not_push_lb_to_benign():
    _, bp6 = _pp5_bp6(_cv(_rec("Benign", RS_3STAR)))
    assert bp6["status"] == "met"
    criteria = [_crit("BS1", strength="BS1_Strong"), bp6]
    total = compute_points_total(criteria)
    assert total == -4
    assert classification_for(total) == "Likely benign"


def test_pp5_moderate_points_value_is_zero():
    pp5, _ = _pp5_bp6(_cv(_rec("Likely pathogenic", RS_2STAR)))
    assert compute_points_total([pp5]) == 0


def test_transparency_string_names_source_and_stars():
    pp5, _ = _pp5_bp6(_cv(_rec("Pathogenic", RS_3STAR, accession="VCV000099999")))
    ev = pp5["evidence"]
    assert "ClinVar" in ev
    assert "VCV000099999" in ev
    assert "★" in ev
    assert "Strong" in ev


def test_pp5_bp6_remain_server_owned():
    assert "PP5" in HARD_CODED_CRITERIA_CODES
    assert "BP6" in HARD_CODED_CRITERIA_CODES


def test_one_star_supporting_is_a_single_switch():
    assert _CLINVAR_STAR_STRENGTH[1] == "Supporting"


def test_pp5_demoted_by_bs3_strong():
    out = _run_exclusions(
        _crit("PP5", strength="PP5_Strong"),
        _crit("BS3", strength="BS3_Strong",
              evidence="Functional patch-clamp assay shows benign; PMID:12345678"),
    )
    assert out["PP5"]["status"] == "not_met"
    assert "BS3" in out["PP5"]["evidence"]


def test_pp5_kept_when_bs3_lacks_pmid():
    out = _run_exclusions(
        _crit("PP5", strength="PP5_Strong"),
        _crit("BS3", strength="BS3_Strong", evidence="benign per ClinVar"),
    )
    assert out["PP5"]["status"] == "met"


def test_pp5_demoted_by_bare_bs4_default_strong():
    out = _run_exclusions(
        _crit("PP5", strength="PP5_Strong"),
        _crit("BS4", strength="BS4"),
    )
    assert out["PP5"]["status"] == "not_met"


def test_pp5_kept_when_benign_only_moderate():
    out = _run_exclusions(
        _crit("PP5", strength="PP5_Strong"),
        _crit("BS3", strength="BS3_Moderate"),
    )
    assert out["PP5"]["status"] == "met"


def test_pp5_kept_when_only_supporting_benign():
    out = _run_exclusions(
        _crit("PP5", strength="PP5_Moderate"),
        _crit("BP4", strength="BP4_Supporting"),
        _crit("BP7", strength="BP7"),
    )
    assert out["PP5"]["status"] == "met"


def test_cross_record_pathogenic_and_benign_suppresses_both():
    cv = _cv(
        _rec("Pathogenic", RS_3STAR, accession="VCV_P"),
        _rec("Benign", RS_3STAR, accession="VCV_B"),
    )
    pp5, bp6 = _pp5_bp6(cv)
    assert pp5["status"] == "not_met"
    assert bp6["status"] == "not_met"
    assert "Conflicting" in pp5["evidence"]
    assert "Conflicting" in bp6["evidence"]


def test_tier_without_clinvar_assertion_strips_pp5():
    pp5, _ = _pp5_bp6(_cv(_rec("Pathogenic", RS_3STAR)))
    criteria = [_crit("PM1", strength="PM1_Moderate"), pp5]
    assert compute_points_total(criteria) == 2
    assert classification_for(compute_points_total(criteria)) == "VUS"
    pts, tier = tier_without_clinvar_assertion(criteria)
    assert pts == 2
    assert tier == "VUS"


def test_tier_without_clinvar_assertion_respects_forced():
    _, bp6 = _pp5_bp6(_cv(_rec("Benign", RS_3STAR)))
    criteria = [_crit("BA1", strength="BA1"), bp6]
    _pts, tier = tier_without_clinvar_assertion(criteria, forced_classification="Benign")
    assert tier == "Benign"


def test_hgvs_match_is_case_insensitive():
    name = "NM_000257.4(MYH7):c.886G>A (p.Gly296Ser)"
    assert _hgvs_match(name, "c.886G>A") is True
    assert _hgvs_match(name, "c.886g>a") is True
    assert _hgvs_match(name, "c.886G>T") is False
    assert _hgvs_match(name, "c.88G>A") is False


if __name__ == "__main__":
    import sys
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as exc:
                failures += 1
                print(f"FAIL {name}: {exc}")
    sys.exit(1 if failures else 0)


def test_ps3_can_no_longer_reach_strong_so_it_cannot_demote_bp6():
    """A consequence of the Brnich 2019 ladder worth pinning rather than
    discovering. The BP6 guard only fires at full Strong/VeryStrong, and PS3's
    maximum earned strength is now Moderate, so PS3 can never trigger it.

    KNOWN UNDER-CALL, deliberate: 11 of the 26 spec'd genes DO publish a PS3
    Strong row, and we cannot reach it because the control count for Strong is in
    Brnich's supplementary OddsPath tables. Conservative direction. Reaching it
    would mean mapping each gene's own PS3 strength rows to assay types."""
    ps3 = _crit("PS3", strength="PS3_Strong",
                evidence="in-vitro functional assay shows loss of function (PMID 31234567)")
    ps3["facts"] = {"assay_variant_controls": 40}
    out = _run_exclusions(_crit("BP6", strength="BP6_Strong"), ps3, gene="MYH7")
    assert out["PS3"]["criteria_strength"] == "PS3_Moderate", (
        "40 controls is far past the 11-control Moderate bar, but Strong is not "
        "reachable without the OddsPath tables"
    )
    assert out["BP6"]["status"] == "met", (
        "a Moderate PS3 must not demote BP6 — only Strong/VeryStrong does"
    )


def test_pm1_collisions_are_resolved_after_the_applicability_gate():
    """Regression for a 2-point inconsistency BETWEEN PATHS.

    `_gate_criteria_applicability` does not only suppress — it ASSERTS PM1 to
    `met` for a variant inside a gene's published hotspot range or exon. That
    gate runs AFTER apply_cross_criterion_exclusions at every call site, so a
    PM1 the gate asserted never had the PM1/PM5 collision resolved. The
    deterministic path scored PM1 + PM5 together (+4) while the AI path, whose
    PM1 arrives before the exclusions, scored PM5 alone (+2).

    Observed on TNNI3 NM_000363.5:c.575G>A: ai_mode=none gave 6 points (Likely
    pathogenic), ai_mode=server gave 5 (VUS), same variant and phenotype. The
    no-AI answer violated the Cardiomyopathy CSpec's "PM5 should not be combined
    with PM1", so it was an OVER-call, and it caused 3 of the 7 variants where
    adding AI made the tier worse.

    The rules are now one function called at both positions."""
    from backend.acmg.hard_coded import apply_pm1_collision_rules

    def _pm(strength_pm1="PM1_Moderate"):
        return [
            {"code": "PM1", "status": "met", "criteria_strength": strength_pm1,
             "direction": "pathogenic", "evidence": "hotspot", "source": "unassessed"},
            {"code": "PM5", "status": "met", "criteria_strength": "PM5_Moderate",
             "direction": "pathogenic", "evidence": "same residue", "source": "hard_coded"},
        ]

    out = apply_pm1_collision_rules(_pm(), "TNNI3")
    got = {c["code"]: c["status"] for c in out}
    assert got == {"PM1": "not_met", "PM5": "met"}, got

    again = apply_pm1_collision_rules(out, "TNNI3")
    assert {c["code"]: c["status"] for c in again} == {"PM1": "not_met", "PM5": "met"}

    pm1 = next(c for c in out if c["code"] == "PM1")
    assert "tier unchanged" not in pm1["evidence"]
    assert "tier CAN fall" in pm1["evidence"]
