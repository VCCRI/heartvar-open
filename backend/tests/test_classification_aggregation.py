"""Unit tests for the tier-deciding aggregation functions in backend.app.

These pure functions turn the combined criteria set into the final
classification, so they are the most accuracy-critical code in the pipeline
and were previously untested:

  - ``_points_for``                  — strength string → signed Tavtigian points
  - ``compute_points_total``         — sum over met criteria
  - ``classification_for``           — points → ACMG tier (boundary values)
  - ``merge_hard_coded_and_ai``      — hard-coded wins, fill missing, canonical order
  - ``apply_cross_criterion_exclusions`` — the mutual-exclusion / guard rules

No pytest dependency — runnable with pytest or directly
(``python -m backend.tests.test_classification_aggregation``).
"""
from __future__ import annotations

from backend.app import (
    _points_for,
    compute_points_total,
    classification_for,
    merge_hard_coded_and_ai,
    build_no_ai_criteria,
    infer_supplementary_criteria,
    apply_cross_criterion_exclusions,
    _CRITERION_NAMES,
    AI_EVALUATED_CRITERIA_CODES,
    CANONICAL_CRITERIA_ORDER,
)


def _crit(code, status="met", strength=None, evidence="", direction=None):
    """Build a criterion dict in the shape the aggregation functions expect.
    A met criterion defaults its strength to the bare code; not_met → None."""
    if strength is None:
        strength = code if status == "met" else None
    return {
        "code": code,
        "name": _CRITERION_NAMES.get(code, code),
        "status": status,
        "direction": direction or ("benign" if code.startswith("B") else "pathogenic"),
        "criteria_strength": strength,
        "evidence": evidence,
        "source": "test",
    }


def test_points_for_tier_suffixes():
    assert _points_for("PVS1_Strong", "PVS1") == 4
    assert _points_for("PVS1_Moderate", "PVS1") == 2
    assert _points_for("PM2_Supporting", "PM2") == 1
    assert _points_for("PM3_Moderate", "PM3") == 2


def test_points_for_benign_tier_suffixes_are_negative():
    assert _points_for("BS1_Strong", "BS1") == -4
    assert _points_for("BP4_Supporting", "BP4") == -1
    assert _points_for("BS2_Strong", "BS2") == -4


def test_points_for_bare_codes():
    assert _points_for("PVS1", "PVS1") == 8
    assert _points_for("PS4", "PS4") == 4
    assert _points_for("PM1", "PM1") == 2
    assert _points_for("PP3", "PP3") == 1
    assert _points_for("BA1", "BA1") == -8
    assert _points_for("BS3", "BS3") == -4
    assert _points_for("BP7", "BP7") == -1


def test_points_for_empty_or_none_is_zero():
    assert _points_for(None, "PS1") == 0
    assert _points_for("", "PS1") == 0


def test_points_total_sums_only_met():
    criteria = [
        _crit("PVS1", "met", "PVS1"),
        _crit("PM2", "met", "PM2_Supporting"),
        _crit("BP4", "met", "BP4_Supporting"),
        _crit("PS1", "not_met"),
    ]
    assert compute_points_total(criteria) == 8


def test_points_total_empty_is_zero():
    assert compute_points_total([]) == 0
    assert compute_points_total(None) == 0


def test_points_total_net_benign():
    criteria = [
        _crit("BA1", "met", "BA1"),
        _crit("BP7", "met", "BP7_Supporting"),
    ]
    assert compute_points_total(criteria) == -9


def test_classification_boundaries():
    assert classification_for(10) == "Pathogenic"
    assert classification_for(11) == "Pathogenic"
    assert classification_for(9) == "Likely pathogenic"
    assert classification_for(6) == "Likely pathogenic"
    assert classification_for(5) == "VUS"
    assert classification_for(0) == "VUS"
    assert classification_for(-1) == "Likely benign"
    assert classification_for(-6) == "Likely benign"
    assert classification_for(-7) == "Benign"


_CANONICAL = (
    "PVS1", "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
    "BA1", "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
)


def test_merge_returns_all_28_in_canonical_order():
    merged = merge_hard_coded_and_ai([], [])
    assert [c["code"] for c in merged] == list(_CANONICAL)
    assert len(merged) == 28


def test_merge_fills_missing_ai_with_not_met_placeholder():
    merged = merge_hard_coded_and_ai([], [])
    ps3 = next(c for c in merged if c["code"] == "PS3")
    assert ps3["status"] == "not_met"
    assert ps3["source"] == "ai"
    assert ps3["criteria_strength"] is None
    assert "AI did not return this criterion" in ps3["evidence"]


def test_merge_does_not_blame_the_ai_for_a_silent_python_rule():
    """The five codes that stopped being sent to the model on 2026-09-08.

    Before the placeholder learned the difference, a PM1/PP1/BS4 the Python
    derivation could not produce came back stamped `source: "ai"` with "AI did
    not return this criterion". Both halves are false, because the model was
    never asked, and the string lands in the report the curator reads.
    """
    merged = {c["code"]: c for c in merge_hard_coded_and_ai([], [])}
    for code in ("PS1", "PM5", "PM1", "PP1", "BS4"):
        entry = merged[code]
        assert entry["status"] == "not_met"
        assert entry["source"] == "hard_coded", (
            f"{code} placeholder is still attributed to the AI"
        )
        assert "AI did not return" not in entry["evidence"], (
            f"{code} placeholder still blames the model for a Python rule"
        )


def test_merge_hard_coded_wins_over_overlapping_ai_code():
    hard_coded = [_crit("PM2", "met", "PM2_Supporting", "absent from gnomAD")]
    ai = [_crit("PM2", "not_met", None, "ai disagrees")]
    merged = merge_hard_coded_and_ai(hard_coded, ai)
    pm2 = next(c for c in merged if c["code"] == "PM2")
    assert pm2["status"] == "met"
    assert pm2["evidence"] == "absent from gnomAD"


def test_merge_keeps_ai_criterion_when_no_hard_coded_collision():
    ai = [_crit("PS1", "met", "PS1", "same-AA pathogenic in ClinVar")]
    merged = merge_hard_coded_and_ai([], ai)
    ps1 = next(c for c in merged if c["code"] == "PS1")
    assert ps1["status"] == "met"
    assert ps1["evidence"] == "same-AA pathogenic in ClinVar"


def test_ba1_forces_benign_and_zeros_pathogenic():
    criteria = [
        _crit("BA1", "met", "BA1"),
        _crit("PVS1", "met", "PVS1"),
        _crit("PM2", "met", "PM2_Supporting"),
    ]
    cleaned, forced = apply_cross_criterion_exclusions(criteria)
    assert forced == "Benign"
    by = {c["code"]: c for c in cleaned}
    assert by["PVS1"]["status"] == "not_met"
    assert by["PM2"]["status"] == "not_met"
    assert by["BA1"]["status"] == "met"


def test_ps2_precludes_pm6():
    criteria = [_crit("PS2", "met", "PS2_Strong"), _crit("PM6", "met", "PM6_Moderate")]
    cleaned, forced = apply_cross_criterion_exclusions(criteria)
    by = {c["code"]: c for c in cleaned}
    assert forced is None
    assert by["PS2"]["status"] == "met"
    assert by["PM6"]["status"] == "not_met"


def test_pp3_and_bp4_conflict_demotes_both():
    criteria = [_crit("PP3", "met", "PP3_Supporting"), _crit("BP4", "met", "BP4_Supporting")]
    cleaned, _ = apply_cross_criterion_exclusions(criteria)
    by = {c["code"]: c for c in cleaned}
    assert by["PP3"]["status"] == "not_met"
    assert by["BP4"]["status"] == "not_met"


def test_pvs1_demotes_pp3():
    criteria = [_crit("PVS1", "met", "PVS1"), _crit("PP3", "met", "PP3_Supporting")]
    cleaned, _ = apply_cross_criterion_exclusions(criteria)
    by = {c["code"]: c for c in cleaned}
    assert by["PVS1"]["status"] == "met"
    assert by["PP3"]["status"] == "not_met"


def test_ps1_precludes_pm5():
    criteria = [_crit("PS1", "met", "PS1"), _crit("PM5", "met", "PM5")]
    cleaned, _ = apply_cross_criterion_exclusions(criteria)
    by = {c["code"]: c for c in cleaned}
    assert by["PS1"]["status"] == "met"
    assert by["PM5"]["status"] == "not_met"


def test_pm4_precludes_bp3():
    criteria = [_crit("PM4", "met", "PM4_Moderate"), _crit("BP3", "met", "BP3_Supporting")]
    cleaned, _ = apply_cross_criterion_exclusions(criteria)
    by = {c["code"]: c for c in cleaned}
    assert by["PM4"]["status"] == "met"
    assert by["BP3"]["status"] == "not_met"


def test_pm1_pm5_kept_when_no_cspec_combine_clause():
    criteria = [_crit("PM1", "met", "PM1_Moderate"), _crit("PM5", "met", "PM5")]
    cleaned, _ = apply_cross_criterion_exclusions(criteria, gene=None)
    by = {c["code"]: c for c in cleaned}
    assert by["PM1"]["status"] == "met"
    assert by["PM5"]["status"] == "met"


def test_pm1_precludes_pm5_for_rasopathy_gene():
    criteria = [_crit("PM1", "met", "PM1_Moderate"), _crit("PM5", "met", "PM5")]
    cleaned, _ = apply_cross_criterion_exclusions(criteria, gene="KRAS")
    by = {c["code"]: c for c in cleaned}
    assert by["PM1"]["status"] == "met"
    assert by["PM5"]["status"] == "not_met", "RASopathy VCEP does not combine PM5 with PM1"


def test_cardiomyopathy_cspec_keeps_pm5_and_demotes_pm1():
    """THE test that would have caught the inverted resolution. The CM CSpec
    does not stop at "PM5 should not be combined with PM1" — it continues "If
    both are applicable at MODERATE weight, use of PM5 is most appropriate
    since it is variant specific". HeartVar kept PM1 and demoted PM5, i.e. it
    kept the region-level argument and threw away the variant-specific one:
    the opposite of what all 8 CM genes (GN002/GN095/GN098/GN099/GN100/GN101/
    GN102/GN103) instruct. Tier-neutral (both +2) but the provenance of the
    call was wrong, and provenance is what a curator audits."""
    for gene in ("MYH7", "MYBPC3", "TNNI3", "TNNT2"):
        criteria = [_crit("PM1", "met", "PM1_Moderate"),
                    _crit("PM5", "met", "PM5_Moderate")]
        cleaned, _ = apply_cross_criterion_exclusions(criteria, gene=gene)
        by = {c["code"]: c for c in cleaned}
        assert by["PM5"]["status"] == "met", (
            f"{gene}: CSpec says PM5 is the appropriate one to use"
        )
        assert by["PM1"]["status"] == "not_met", f"{gene}: PM1 must be demoted"
        assert "variant specific" in by["PM1"]["evidence"]
        assert compute_points_total(cleaned) == 2


def test_cm_tie_break_keeps_the_stronger_criterion_when_strengths_differ():
    """The CSpec rule is stated for the MODERATE/MODERATE case. When PM1 is
    stronger than PM5 — which happens after the PP3+PM1 correlated-evidence
    cap has already touched PM1, or if a future spec grades PM1 up — following
    the letter of the rule would silently discard a point. Keep the stronger
    one and label the tie-break as HeartVar's."""
    criteria = [_crit("PM1", "met", "PM1_Moderate"),
                _crit("PM5", "met", "PM5_Supporting")]
    cleaned, _ = apply_cross_criterion_exclusions(criteria, gene="MYH7")
    by = {c["code"]: c for c in cleaned}
    assert by["PM1"]["status"] == "met"
    assert by["PM5"]["status"] == "not_met"
    assert "HeartVar's tie-break" in by["PM5"]["evidence"]


def test_recessive_ps4_without_case_control_is_demoted():
    criteria = [_crit("PS4", "met", "PS4", "3 unrelated probands carried the variant")]
    cleaned, _ = apply_cross_criterion_exclusions(criteria, gene=None, inheritance="AR")
    assert {c["code"]: c for c in cleaned}["PS4"]["status"] == "not_met"


def test_recessive_ps4_with_case_control_is_kept():
    """The escape from the recessive demotion now needs the STRUCTURED
    facts.or_lower_ci, not odds-ratio language in the prose."""
    c = _crit("PS4", "met", "PS4", "case-control study, see PMID 12345678")
    c["facts"] = {"or_lower_ci": 2.0}
    cleaned, _ = apply_cross_criterion_exclusions([c], gene=None, inheritance="AR")
    assert {x["code"]: x for x in cleaned}["PS4"]["status"] == "met"


def test_recessive_ps4_prose_odds_ratio_no_longer_escapes_the_demotion():
    """The sibling of the bug A2 fixed on the main PS4 route. This regex was the
    ESCAPE from the demotion, so a false positive GRANTED weight — including on
    prose that says an odds ratio was NOT reported. A case-control PS4 needs an
    actual machine-readable statistic."""
    for prose in (
        "case-control odds ratio 4.1 (95% CI 2.0-8.3)",
        "no odds ratio was reported for this variant",
        "either the proband or: a sibling carried it",
        "95% CI reported for the allele frequency, not a case-control test",
    ):
        criteria = [_crit("PS4", "met", "PS4", prose)]
        cleaned, _ = apply_cross_criterion_exclusions(
            criteria, gene=None, inheritance="AR",
        )
        got = {c["code"]: c for c in cleaned}["PS4"]
        assert got["status"] == "not_met", prose


def test_recessive_ps4_ignores_a_non_numeric_or_lower_ci():
    """Fail closed: a truthy-but-unusable value must not buy the escape."""
    for bad in (True, "2.0", None, 0, -1):
        c = _crit("PS4", "met", "PS4", "case-control")
        c["facts"] = {"or_lower_ci": bad}
        cleaned, _ = apply_cross_criterion_exclusions(
            [c], gene=None, inheritance="AR",
        )
        assert {x["code"]: x for x in cleaned}["PS4"]["status"] == "not_met", bad


def test_dominant_ps4_proband_count_is_kept():
    criteria = [_crit("PS4", "met", "PS4", "5 unrelated probands")]
    cleaned, _ = apply_cross_criterion_exclusions(criteria, gene=None, inheritance="AD")
    assert {c["code"]: c for c in cleaned}["PS4"]["status"] == "met"


def test_bs3_without_assay_and_pmid_is_demoted():
    criteria = [_crit("BS3", "met", "BS3_Strong", "benign in ClinVar, high allele frequency")]
    cleaned, _ = apply_cross_criterion_exclusions(criteria)
    assert {c["code"]: c for c in cleaned}["BS3"]["status"] == "not_met"


def test_bs3_with_assay_and_pmid_is_kept():
    criteria = [_crit(
        "BS3", "met", "BS3_Strong",
        "patch-clamp functional assay showed normal channel activity (PMID 30000001)",
    )]
    cleaned, _ = apply_cross_criterion_exclusions(criteria)
    assert {c["code"]: c for c in cleaned}["BS3"]["status"] == "met"


def test_ps3_without_assay_and_pmid_is_demoted():
    for bad in (
        "pathogenic in ClinVar with 3-star review",
        "REVEL 0.92 and located in a mutational hotspot",
        "a different pathogenic missense at this residue (PM5)",
        "SpliceAI predicts a donor-loss splice effect",
    ):
        criteria = [_crit("PS3", "met", "PS3_Strong", bad)]
        cleaned, _ = apply_cross_criterion_exclusions(criteria)
        assert {c["code"]: c for c in cleaned}["PS3"]["status"] == "not_met", bad


def test_ps3_with_assay_and_pmid_is_kept():
    """The VERIFICATION gate: a real assay + PMID survives. Validation facts are
    attached because the Brnich strength ladder is a separate later step and
    would otherwise suppress these before the gate is reached."""
    for good in (
        "patch-clamp electrophysiology showed loss of current (PMID 31000002)",
        "in-vitro expression study showed loss of function (PMID:29000003)",
    ):
        c = _crit("PS3", "met", "PS3_Strong", good)
        c["facts"] = dict(VALIDATED)
        cleaned, _ = apply_cross_criterion_exclusions([c], gene="MYH7")
        assert {x["code"]: x for x in cleaned}["PS3"]["status"] == "met", good


VALIDATED = {"assay_lab_controls": True}


def _ps3_facts(evidence, facts, strength="PS3_Strong", gene="MYH7",
               validated=True):
    """gene defaults to MYH7 because GN002 publishes a PS3 Strong row, so the
    gene ceiling does not interfere. `validated` merges in VALIDATED so these
    tests exercise the verification gate rather than the strength ladder; set it
    False to test the ladder itself."""
    c = _crit("PS3", "met", strength, evidence)
    if facts is None and validated:
        facts = dict(VALIDATED)
    elif isinstance(facts, dict) and validated:
        facts = {**facts, **VALIDATED}
    c["facts"] = facts
    cleaned, _ = apply_cross_criterion_exclusions([c], gene=gene)
    return {x["code"]: x for x in cleaned}["PS3"]


def test_ps3_facts_keep_when_prose_keyword_absent():
    ps3 = _ps3_facts(
        "This variant abolished channel current in the referenced study (PMID 31000002).",
        {"assay_pmids": [31000002], "assay_type": "patch-clamp electrophysiology"},
    )
    assert ps3["status"] == "met"
    assert ps3["criteria_strength"] == "PS3_Supporting"


def test_ps3_facts_unrecognized_assay_type_falls_back_to_prose():
    ps3 = _ps3_facts(
        "Observed in a clinical case (PMID 31000002).",
        {"assay_pmids": [31000002], "assay_type": "clinical observation"},
    )
    assert ps3["status"] == "not_met"


def test_ps3_facts_non_numeric_pmid_falls_back():
    ps3 = _ps3_facts(
        "Damaging in a study (PMID n/a).",
        {"assay_pmids": ["n/a"], "assay_type": "patch-clamp electrophysiology"},
    )
    assert ps3["status"] == "not_met"


def test_ps3_facts_pmid_not_echoed_falls_back_and_demotes():
    ps3 = _ps3_facts(
        "Predicted deleterious by in-silico tools.",
        {"assay_pmids": [31000002], "assay_type": "patch-clamp electrophysiology"},
    )
    assert ps3["status"] == "not_met"


def test_ps3_facts_pmid_echoed_and_good_prose_still_kept():
    ps3 = _ps3_facts(
        "Patch-clamp electrophysiology showed loss of current (PMID 31000002).",
        {"assay_pmids": [31000002], "assay_type": "patch-clamp electrophysiology"},
    )
    assert ps3["status"] == "met"


def test_ps3_facts_none_is_byte_identical_to_prose_gate():
    good = _ps3_facts("patch-clamp electrophysiology showed loss of current (PMID 31000002)", None)
    bad = _ps3_facts("pathogenic in ClinVar with 3-star review", None)
    assert good["status"] == "met"
    assert bad["status"] == "not_met"


def test_ps3_is_not_applied_when_the_assay_validation_is_undescribed():
    """The Brnich 2019 starting point, and the single biggest behaviour change of
    2026-09-01: "evaluation of functional assays should start from the assumption
    of no evidence". A real assay with a real PMID earns NOTHING unless the paper
    describes its validation, because a default strength is exactly what the SVI
    framework replaced. Conservative by construction — a curator who has read the
    paper can apply PS3 by hand."""
    undescribed = _ps3_facts(
        "patch-clamp electrophysiology showed loss of current (PMID 31000002)",
        {"assay_pmids": [31000002], "assay_type": "patch-clamp electrophysiology"},
        validated=False,
    )
    assert undescribed["status"] == "not_met"
    assert "start from the assumption of no evidence" in undescribed["evidence"]


def _ps4_strength(strength, evidence, gene="MYH7", inh="AD"):
    criteria = [_crit("PS4", "met", strength, evidence)]
    cleaned, _ = apply_cross_criterion_exclusions(criteria, gene=gene, inheritance=inh)
    return {c["code"]: c for c in cleaned}["PS4"]


def test_ps4_capped_to_supporting_when_few_probands():
    for strength in ("PS4_Strong", "PS4", "PS4_Moderate"):
        ps4 = _ps4_strength(strength, "identified in 5 unrelated probands with HCM")
        assert ps4["status"] == "met"
        assert ps4["criteria_strength"] == "PS4_Supporting", strength


def test_ps4_uses_the_rasopathy_ladder_on_lztr1_and_mras():
    """LZTR1 and MRAS fell through BOTH PS4 cap branches and got no cap at all.

    The RASopathy arm is selected by `_vcep_freq(gene)["vcep"] == "RASopathy"`,
    true for all 16 RASopathy genes — but the per-gene ladder was read only from
    each gene's PS4 *Comments*, and LZTR1 and MRAS are the two RASopathy genes
    whose PS4 record carries none. The own-ladder branch was skipped for want of
    a ladder and the generic Kelly branch was skipped by `not _ps4_is_raso`, so
    a met PS4 kept whatever the model emitted and a bare "PS4" scored Strong
    (+4). LZTR1 c.848G>A scored +12 Pathogenic against the panel's LP.

    They DO publish a ladder — as POINTS, in the strength rows rather than in
    Comments:

        Supporting  ">=1 points."   Moderate  ">=3 points."   Strong  ">=5 points."

    and the conversion is published too. LZTR1's vcep_specifications: "full
    points (1) awarded with consistent RASopathy phenotype". One point per
    fully-phenotyped case, so those thresholds are the SAME >=1 / >=3 / >=5 the
    other 14 RASopathy genes publish as "independent occurrences" — not the
    Kelly >=2/>=6/>=15 cardiomyopathy adaptation, which is a different panel's
    rule and would be wrong here.

    ⚠ Deliberately NOT generalised to any gene with a points ladder. FBN1's PS4
    rungs are also points ("If >= 4 points.") but FBN1 publishes no
    points-per-proband mapping — that is a CLOSED gap in
    SPEC_criteria_completeness.md. The conversion is gated on the RASopathy
    VCEP, which publishes the per-case value.
    """
    for gene in ("LZTR1", "MRAS"):
        ps4 = _ps4_strength("PS4", "identified in 1 proband with Noonan syndrome",
                            gene=gene)
        assert ps4["status"] == "met", (gene, ps4)
        assert ps4["criteria_strength"] == "PS4_Supporting", (gene, ps4["criteria_strength"])

        ps4 = _ps4_strength("PS4", "identified in 3 unrelated probands", gene=gene)
        assert ps4["criteria_strength"] == "PS4_Moderate", (gene, ps4["criteria_strength"])

        ps4 = _ps4_strength("PS4", "identified in 5 unrelated probands", gene=gene)
        assert ps4["criteria_strength"] in (None, "PS4"), (gene, ps4["criteria_strength"])

    from backend.acmg.hard_coded import _ps4_spec_occurrence_ladder
    assert _ps4_spec_occurrence_ladder("FBN1") == {}


def test_ps4_graded_on_kcnq1_and_bmpr2_own_published_ladders():
    """Before this, KCNQ1 and BMPR2 were graded on the Kelly >=2/>=6/>=15
    ladder because their thresholds live in the STRENGTH ROWS rather than in
    Comments, and the reader only looked at Comments. Their own published
    rungs are stricter in places and looser in others, so the Kelly fallback
    was wrong in both directions.

      KCNQ1  >=2 Supporting / 3-5 Moderate / >=6 Strong
      BMPR2  >1 Supporting / >3 Moderate / >4 Strong   (strict bounds)
    """
    ps4 = _ps4_strength("PS4", "identified in 3 unrelated probands with LQTS",
                        gene="KCNQ1")
    assert ps4["status"] == "met"
    assert ps4["criteria_strength"] == "PS4_Moderate", ps4["criteria_strength"]

    ps4 = _ps4_strength("PS4", "identified in 6 unrelated probands with LQTS",
                        gene="KCNQ1")
    assert ps4["criteria_strength"] in (None, "PS4"), ps4["criteria_strength"]

    ps4 = _ps4_strength("PS4", "reported in 5 unrelated patients with PAH",
                        gene="BMPR2")
    assert ps4["criteria_strength"] in (None, "PS4"), ps4["criteria_strength"]

    for gene in ("KCNQ1", "BMPR2"):
        ps4 = _ps4_strength("PS4", "identified in 1 proband", gene=gene)
        assert ps4["status"] == "not_met", (gene, ps4)


def test_ps4_moderate_at_six_to_fourteen_probands():
    ps4 = _ps4_strength("PS4_Strong", "reported in 8 unrelated probands with DCM")
    assert ps4["criteria_strength"] == "PS4_Moderate"


def test_ps4_strong_kept_at_fifteen_plus_probands():
    ps4 = _ps4_strength("PS4_Strong", "observed in 22 unrelated probands with HCM")
    assert ps4["criteria_strength"] == "PS4_Strong"


def test_ps4_case_control_route_needs_a_structured_statistic_not_prose():
    """The case-control route skips the ENTIRE proband cap and bare PS4 then
    defaults to Strong (+4), so what opens that route matters. It used to be a
    prose regex (`odds[\\s-]*ratio|\\bor\\s*[=:]|95\\s*%?\\s*ci|...`) — which the
    phrase "odds ratio" satisfies wherever it appears, including in a sentence
    saying no odds ratio was reported. It now requires the structured
    `facts.or_lower_ci`. Prose alone falls through to the proband ladder, and
    with no count quoted that means Supporting."""
    ps4 = _ps4_strength("PS4_Strong",
                        "case-control odds ratio 8.2 (95% CI 3.1-21.7) in HCM cohort")
    assert ps4["criteria_strength"] == "PS4_Supporting"


def test_ps4_prose_denying_an_odds_ratio_no_longer_earns_strong():
    """The test that would have caught the bug. "no odds ratio" contains "odds
    ratio", so the old regex read a DENIAL of case-control evidence as proof of
    it and handed the variant an uncapped PS4_Strong (+4) — 4 points, on the
    pathogenic side, from a negative statement."""
    ps4 = _ps4_strength(
        "PS4_Strong",
        "Seen in 2 probands; no odds ratio or 95% CI was reported in the source.",
    )
    assert ps4["criteria_strength"] == "PS4_Supporting", (
        "prose denying case-control evidence must not open the case-control "
        "route"
    )


def test_ps4_unquantified_defaults_to_supporting():
    ps4 = _ps4_strength("PS4_Strong", "reported in individuals with cardiomyopathy")
    assert ps4["criteria_strength"] == "PS4_Supporting"


def test_ps4_never_raised():
    ps4 = _ps4_strength("PS4_Supporting", "seen in 30 unrelated probands")
    assert ps4["criteria_strength"] == "PS4_Supporting"


def _ps4_facts(strength, evidence, facts, gene="MYH7", inh="AD"):
    c = _crit("PS4", "met", strength, evidence)
    c["facts"] = facts
    cleaned, _ = apply_cross_criterion_exclusions([c], gene=gene, inheritance=inh)
    return {x["code"]: x for x in cleaned}["PS4"]


def test_ps4_facts_proband_count_reaches_strong_with_pmid():
    ps4 = _ps4_facts("PS4", "Case series across families; PMID 21622575",
                     {"proband_count": 22})
    assert compute_points_total([ps4]) == 4


def test_ps4_facts_proband_count_moderate_with_pmid():
    ps4 = _ps4_facts("PS4", "Reported cohort; PMID 12345678", {"proband_count": 8})
    assert ps4["criteria_strength"] == "PS4_Moderate"
    assert compute_points_total([ps4]) == 2


def test_ps4_facts_count_without_pmid_held_at_supporting():
    ps4 = _ps4_facts("PS4", "Reported in a large cohort", {"proband_count": 40})
    assert ps4["criteria_strength"] == "PS4_Supporting"
    assert compute_points_total([ps4]) == 1


def test_ps4_facts_or_lower_ci_marks_case_control_route():
    ps4 = _ps4_facts("PS4", "Case-control enrichment vs controls",
                     {"or_lower_ci": 8.2})
    assert compute_points_total([ps4]) == 4


def test_ps4_facts_absent_falls_back_to_prose():
    ps4 = _ps4_facts("PS4", "reported in 8 unrelated probands with DCM", None)
    assert ps4["criteria_strength"] == "PS4_Moderate"


from backend.acmg.hard_coded import _pp1_strength_from_lod_meioses  # noqa: E402


def _pp1_lit(strength, evidence, facts, gene="MYH7", inh="AD",
            has_curator_segregation=False):
    """Run a met PP1 (source='test' → literature-sourced) through the
    exclusions block and return the resulting PP1 dict."""
    c = _crit("PP1", "met", strength, evidence)
    c["facts"] = facts
    cleaned, _ = apply_cross_criterion_exclusions(
        [c], gene=gene, inheritance=inh,
        has_curator_segregation=has_curator_segregation,
    )
    return {x["code"]: x for x in cleaned}["PP1"]


def test_pp1_ladder_kelly_maps_supporting_moderate_strong():
    assert _pp1_strength_from_lod_meioses(None, 3, 2, 0) == "Supporting"
    assert _pp1_strength_from_lod_meioses(None, 5, 3, 0) == "Moderate"
    assert _pp1_strength_from_lod_meioses(None, 7, 4, 0) == "Strong"
    assert _pp1_strength_from_lod_meioses(2.1, 0, 1, 0) == "Strong"
    assert _pp1_strength_from_lod_meioses(1.5, 0, 1, 0) == "Moderate"
    assert _pp1_strength_from_lod_meioses(0.9, 0, 1, 0) == "Supporting"
    assert _pp1_strength_from_lod_meioses(None, 0, 2, 0) == "Supporting"
    assert _pp1_strength_from_lod_meioses(None, 0, 1, 0) is None
    assert _pp1_strength_from_lod_meioses(None, 0, 0, 0) is None


def test_pp1_helper_parity_with_no_ai_path():
    for carriers in range(0, 4):
        for meioses in range(0, 9):
            tier = _pp1_strength_from_lod_meioses(None, meioses, carriers, 0)
            if carriers >= 1:
                if meioses >= 7:
                    expect = "Strong"
                elif meioses >= 5:
                    expect = "Moderate"
                elif meioses >= 3:
                    expect = "Supporting"
                elif carriers >= 2:
                    expect = "Supporting"
                else:
                    expect = None
            else:
                expect = None
            assert tier == expect, (carriers, meioses, tier, expect)


def test_pp1_literature_no_pmid_is_demoted():
    pp1 = _pp1_lit("PP1", "Co-segregated in the family",
                   {"seg_affected_carriers": 4, "seg_meioses": 6})
    assert pp1["status"] == "not_met"
    assert "PMID" in pp1["evidence"]


def test_pp1_literature_with_pmid_retained():
    pp1 = _pp1_lit("PP1_Moderate", "Segregation reported; PMID 20301486",
                   {"seg_affected_carriers": 4, "seg_meioses": 5, "pmid": 20301486})
    assert pp1["status"] == "met"
    assert pp1["criteria_strength"] == "PP1_Moderate"


def test_pp1_literature_needs_one_genotyped_carrier():
    pp1 = _pp1_lit("PP1", "Family study; PMID 20301486",
                   {"seg_affected_carriers": 0, "seg_meioses": 8, "pmid": 20301486})
    assert pp1["status"] == "not_met"
    assert "genotyped affected carrier" in pp1["evidence"].lower()


def test_pp1_literature_single_carrier_no_meioses_demoted():
    pp1 = _pp1_lit("PP1", "One relative also carries it; PMID 20301486",
                   {"seg_affected_carriers": 1, "pmid": 20301486})
    assert pp1["status"] == "not_met"


def test_pp1_literature_cap_only_never_raises():
    pp1 = _pp1_lit("PP1", "Segregation; PMID 20301486",
                   {"seg_affected_carriers": 5, "seg_meioses": 7, "pmid": 20301486})
    assert pp1["status"] == "met"
    assert pp1["criteria_strength"] == "PP1_Supporting"


def test_pp1_literature_capped_at_moderate_even_with_high_lod():
    pp1 = _pp1_lit("PP1_Strong", "Strong segregation, LOD 3.0; PMID 20301486",
                   {"seg_affected_carriers": 8, "lod": 3.0, "pmid": 20301486})
    assert pp1["status"] == "met"
    assert pp1["criteria_strength"] == "PP1_Moderate"


def test_pp1_literature_llm_strong_facts_supporting_capped_down():
    pp1 = _pp1_lit("PP1_Strong", "Segregation; PMID 20301486",
                   {"seg_affected_carriers": 3, "seg_meioses": 3, "pmid": 20301486})
    assert pp1["status"] == "met"
    assert pp1["criteria_strength"] == "PP1_Supporting"


def test_pp1_unaffected_carrier_tempers_one_tier():
    pp1 = _pp1_lit("PP1_Strong", "Segregation with one unaffected carrier; PMID 20301486",
                   {"seg_affected_carriers": 5, "seg_meioses": 7,
                    "seg_unaffected_carriers": 1, "pmid": 20301486})
    assert pp1["status"] == "met"
    assert pp1["criteria_strength"] == "PP1_Moderate"


def test_pp1_curator_path_not_demoted_without_pmid():
    pp1 = _pp1_lit("PP1_Strong", "Curator segregation: 4 carriers across 7 meioses",
                   {"seg_affected_carriers": 4, "seg_meioses": 7},
                   has_curator_segregation=True)
    assert pp1["status"] == "met"
    assert pp1["criteria_strength"] == "PP1_Strong"


def test_pp1_inferred_source_immune_to_literature_hardening():
    c = _crit("PP1", "met", "PP1_Strong", "Co-segregation (curator-entered).")
    c["source"] = "inferred"
    cleaned, _ = apply_cross_criterion_exclusions([c], gene="MYH7", inheritance="AD")
    pp1 = {x["code"]: x for x in cleaned}["PP1"]
    assert pp1["status"] == "met" and pp1["criteria_strength"] == "PP1_Strong"


_HC_NOT_EMITTED_BY_THE_ENGINE = ("PM1", "PP1", "BS4")


def _hc(*entries):
    """Assemble a hard-coded-style set: every code defaults to not_met unless
    listed in ``entries`` as (code, strength).

    Mirrors what compute_hard_coded_criteria actually emits, so the three codes
    it owns in name only (see above) are left out and stay available to the
    supplementary layer."""
    from backend.app import HARD_CODED_CRITERIA_CODES
    met = {code: strength for code, strength in entries}
    out = []
    for code in HARD_CODED_CRITERIA_CODES:
        if code in _HC_NOT_EMITTED_BY_THE_ENGINE and code not in met:
            continue
        if code in met:
            out.append(_crit(code, "met", met[code]))
        else:
            out.append(_crit(code, "not_met", None))
    return out


def test_build_no_ai_pads_to_28_codes():
    full = build_no_ai_criteria(_hc(), [])
    assert [c["code"] for c in full] == list(CANONICAL_CRITERIA_ORDER)
    assert len(full) == 28


def test_build_no_ai_tags_ai_codes_not_assessed():
    full = build_no_ai_criteria(_hc(), [])
    not_assessed = {c["code"] for c in full if c["status"] == "not_assessed"}
    assert not_assessed == set(AI_EVALUATED_CRITERIA_CODES) | {"PM1", "PP1", "BS4"}
    for c in full:
        if c["status"] == "not_assessed":
            assert c["source"] == "unassessed"


def test_build_no_ai_not_assessed_scores_zero():
    full = build_no_ai_criteria(_hc(), [])
    assert compute_points_total(full) == 0
    assert classification_for(compute_points_total(full)) == "VUS"


def test_build_no_ai_supplementary_fills_a_code():
    supp = [_crit("PP1", "met", "PP1_Strong", "co-segregation in 5 meioses")]
    full = build_no_ai_criteria(_hc(), supp)
    by_code = {c["code"]: c for c in full}
    assert by_code["PP1"]["status"] == "met"
    assert by_code["PP1"]["criteria_strength"] == "PP1_Strong"


def test_build_no_ai_supplementary_never_overrides_hard_coded():
    supp = [_crit("PM2", "met", "PM2_Moderate", "should be ignored")]
    full = build_no_ai_criteria(_hc(("PM2", "PM2_Supporting")), supp)
    by_code = {c["code"]: c for c in full}
    assert by_code["PM2"]["criteria_strength"] == "PM2_Supporting"


def test_no_ai_pvs1_pm2_is_likely_pathogenic():
    full = build_no_ai_criteria(_hc(("PVS1", "PVS1_VeryStrong"), ("PM2", "PM2_Supporting")), [])
    cleaned, forced = apply_cross_criterion_exclusions(full, gene="MYBPC3", inheritance="AD")
    pts = compute_points_total(cleaned)
    assert pts == 9
    assert (forced or classification_for(pts)) == "Likely pathogenic"


def test_no_ai_ba1_forces_benign_even_with_pathogenic_met():
    full = build_no_ai_criteria(_hc(("BA1", "BA1"), ("PVS1", "PVS1_VeryStrong")), [])
    cleaned, forced = apply_cross_criterion_exclusions(full, gene="TTN", inheritance="AR")
    assert forced == "Benign"
    assert {c["code"]: c for c in cleaned}["PVS1"]["status"] == "not_met"


def test_no_ai_frequency_only_is_likely_benign():
    full = build_no_ai_criteria(_hc(("BS1", "BS1_Strong")), [])
    cleaned, forced = apply_cross_criterion_exclusions(full)
    pts = compute_points_total(cleaned)
    assert pts == -4
    assert (forced or classification_for(pts)) == "Likely benign"


def _supp_codes(out):
    return {c["code"]: c for c in out}


def test_infer_supplementary_empty_inputs_returns_nothing():
    assert infer_supplementary_criteria({}, {}, "MYH7", _hc()) == []


def test_infer_pp1_from_segregation_strength_ladder():
    def pp1(carriers, meioses):
        out = infer_supplementary_criteria(
            {}, {"seg_affected_carriers": carriers, "seg_meioses": meioses}, "MYH7", _hc())
        return _supp_codes(out).get("PP1", {}).get("criteria_strength")
    assert pp1(2, 7) == "PP1_Strong"
    assert pp1(2, 5) == "PP1_Moderate"
    assert pp1(2, 3) == "PP1_Supporting"
    assert pp1(2, 0) == "PP1_Supporting"
    assert pp1(1, 0) is None
    assert pp1(0, 9) is None


def test_infer_bs4_requires_two_affected_noncarriers():
    def bs4(n):
        out = infer_supplementary_criteria({}, {"seg_affected_noncarriers": n}, "MYH7", _hc())
        return _supp_codes(out).get("BS4", {}).get("criteria_strength")
    assert bs4(1) is None
    assert bs4(2) == "BS4"


def test_infer_bs4_strength_matches_the_ai_path():
    """The no-key path emitted BS4_Supporting (-1) while the AI path emitted bare
    BS4 (-4), so identical curator input scored 3 points apart depending only on
    whether an API key was configured.

    Biesecker et al. 2024 (ClinGen guidance for PP1/BS4 and PP4, AJHG
    111:24-38) settles it: "In autosomal-dominant inheritance,
    autosomal-recessive inheritance with homozygosity, and X-linked inheritance,
    non-segregations do provide evidence of benignity ... we recommend assigning
    -4.0 points to such variants for these observations." So the AI path was
    right and the cap on this path was the bug.
    """
    from backend.acmg.tiers import _points_for
    out = infer_supplementary_criteria(
        {}, {"seg_affected_noncarriers": 3}, "MYH7", _hc(),
    )
    bs4 = _supp_codes(out)["BS4"]
    assert bs4["criteria_strength"] == "BS4"
    assert _points_for(bs4["criteria_strength"], "BS4") == -4, (
        "Biesecker 2024 specifies -4.0 points for a non-segregation observation"
    )


def test_infer_pm1_from_domain_plp():
    met = {"ok": True, "domain_name": "kinase", "domain_start": 1, "domain_end": 9,
           "total_plp": 4, "has_two_star_plp": True}
    sup = {"ok": True, "domain_name": "kinase", "domain_start": 1, "domain_end": 9,
           "total_plp": 1, "has_two_star_plp": False}
    none = {"ok": True, "domain_name": "kinase", "domain_start": 1, "domain_end": 9,
            "total_plp": 0, "has_two_star_plp": False}
    def pm1(dp):
        out = infer_supplementary_criteria({"domain_plp": dp}, {}, "MYH7", _hc())
        return _supp_codes(out).get("PM1")
    assert pm1(met)["criteria_strength"] == "PM1_Moderate"
    assert pm1(sup)["criteria_strength"] == "PM1_Supporting"
    assert pm1(none)["status"] == "not_met"


def test_infer_ps1_wins_over_pm5_when_both_present():
    both = {"clinvar_pm5_candidates": {"ok": True, "ps1_count_two_star": 1, "count_two_star": 2}}
    out = _supp_codes(infer_supplementary_criteria(both, {}, "MYH7", _hc()))
    assert out.get("PS1", {}).get("criteria_strength") == "PS1_Strong"
    assert "PM5" not in out
    pm5only = {"clinvar_pm5_candidates": {"ok": True, "ps1_count_two_star": 0, "count_two_star": 2}}
    out2 = _supp_codes(infer_supplementary_criteria(pm5only, {}, "MYH7", _hc()))
    assert out2.get("PM5", {}).get("criteria_strength") == "PM5_Moderate"
    assert "PS1" not in out2


def test_infer_clinvar_below_two_star_does_not_fire():
    weak = {"clinvar_pm5_candidates": {"ok": True, "ps1_count_two_star": 0, "count_two_star": 0}}
    assert infer_supplementary_criteria(weak, {}, "MYH7", _hc()) == []


def test_supplementary_pp1_flows_into_no_ai_score():
    supp = infer_supplementary_criteria(
        {}, {"seg_affected_carriers": 3, "seg_meioses": 7}, "MYBPC3", _hc())
    full = build_no_ai_criteria(_hc(("PM2", "PM2_Supporting"), ("PVS1", "PVS1_Strong")), supp)
    cleaned, forced = apply_cross_criterion_exclusions(full, gene="MYBPC3", inheritance="AD")
    pts = compute_points_total(cleaned)
    assert pts == 9
    assert (forced or classification_for(pts)) == "Likely pathogenic"


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()


def test_ps3_only_bmpr2_publishes_a_no_controls_floor():
    """🔴 PIN. This is the guard that keeps the fix scoped.

    If a re-harvest makes another gene match, this fails and a human decides
    whether that gene really publishes a FLOOR or whether the regex has started
    matching a CEILING. The Cardiomyopathy row is the trap: "data from
    individual in vitro studies are unlikely to meet the criteria required to
    assign this rule at MORE THAN supporting level" is a maximum, and an earlier
    attempt at this fix quoted it as authority for a minimum.
    """
    import json as _json
    from pathlib import Path as _Path
    from backend.acmg.hard_coded import _ps3_spec_floor

    spec = _json.loads(
        (_Path(__file__).resolve().parents[1] / "data" / "vcep_criteria_spec.json")
        .read_text()
    )
    matched = {g: _ps3_spec_floor(g) for g in spec["genes"]}
    assert {g: v for g, v in matched.items() if v} == {"BMPR2": "Supporting"}, (
        "the PS3 no-controls floor is scoped to BMPR2; another gene now matches"
    )
    assert len(matched) >= 26


def test_ps3_the_cardiomyopathy_ceiling_is_not_read_as_a_floor():
    """MYH7's Supporting row caps one assay class; it licenses nothing."""
    from backend.acmg.hard_coded import _ps3_spec_floor, _ps3_spec_ceiling
    assert _ps3_spec_floor("MYH7") is None
    assert _ps3_spec_ceiling("MYH7") == "Strong"


def test_ps3_fbn1_delegates_to_brnich_so_it_keeps_the_no_evidence_start():
    """FBN1's three PS3 rows say only "Follow the 'Funtional Assay SVI
    Documentation'", so for FBN1 the gene spec IS Brnich 2019 and the current
    demote-to-zero behaviour is already right."""
    from backend.acmg.hard_coded import _ps3_spec_floor
    assert _ps3_spec_floor("FBN1") is None
    ps3 = _ps3_facts(
        "Reduced fibrillin-1 deposition in a functional assay (PMID 31000002).",
        {"assay_pmids": [31000002], "assay_type": "protein expression"},
        gene="FBN1", validated=False,
    )
    assert ps3["status"] == "not_met"


def test_ps3_bmpr2_holds_at_supporting_when_controls_are_undescribed():
    """The behaviour change. BMPR2 GN125 Supporting row, VERBATIM: "If no known
    variant validation controls (i.e. established pathogenic and benign
    variants) were used, then score at the supporting strength.\""""
    ps3 = _ps3_facts(
        "BMP signalling was abolished in a reporter assay (PMID 31000002).",
        {"assay_pmids": [31000002], "assay_type": "reporter assay"},
        strength="PS3_Strong", gene="BMPR2", validated=False,
    )
    assert ps3["status"] == "met"
    assert ps3["criteria_strength"] == "PS3_Supporting"
    assert "BMPR2" in ps3["evidence"]
    assert "no known variant validation controls" in ps3["evidence"]


def test_ps3_bmpr2_floor_does_not_raise_a_validated_assay():
    """The floor is a floor. An assay that earns Moderate keeps Moderate."""
    ps3 = _ps3_facts(
        "Abolished signalling (PMID 31000002).",
        {"assay_pmids": [31000002], "assay_type": "reporter assay",
         "assay_variant_controls": 12},
        strength="PS3_Strong", gene="BMPR2", validated=False,
    )
    assert ps3["criteria_strength"] == "PS3_Moderate"


def test_ps3_bmpr2_floor_is_still_bounded_by_the_gene_ceiling():
    """Floor and ceiling both apply; PS3 takes the lower. BMPR2's ceiling is
    Strong, so the floor is what binds here — this test exists so that a future
    gene with a Supporting ceiling and a Moderate floor cannot silently exceed
    its ceiling."""
    from backend.acmg.hard_coded import _ps3_spec_floor, _ps3_spec_ceiling
    floor, ceiling = _ps3_spec_floor("BMPR2"), _ps3_spec_ceiling("BMPR2")
    from backend.acmg.hard_coded import _PS3_RANK
    assert _PS3_RANK[f"PS3_{floor}"] <= _PS3_RANK[f"PS3_{ceiling}"]


def test_ps3_non_bmpr2_genes_are_unchanged_by_the_floor():
    """Blast radius. Every other spec'd gene still demotes on undescribed
    controls, which is the pre-2026-09-07 behaviour."""
    for gene in ("MYH7", "MYBPC3", "KCNQ1", "PTPN11", "TNNT2"):
        ps3 = _ps3_facts(
            "Damaging in a functional study (PMID 31000002).",
            {"assay_pmids": [31000002], "assay_type": "patch-clamp electrophysiology"},
            gene=gene, validated=False,
        )
        assert ps3["status"] == "not_met", f"{gene} changed behaviour"
