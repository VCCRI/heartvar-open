"""PS2 / PM6 / PS4 graded on the rules the registry publishes.

All three were previously scored by a hardcoded rule that overrode the gene's
own ladder, and all three were justified by the SAME wrong belief: that the
operational thresholds were not available. They were — in the registry's
"Comments" and "VCEP Specifications" fields, which the harvester did not read
until 2026-09-07.

  PS2/PM6  the SVI de novo points table gives points per occurrence, and the
           RASopathy panel says which row to use ("full points awarded with
           RASopathy phenotypes"). GN043's Comments state the occurrence rule
           outright: "PS2_Very Strong: >=2 independent occurrences of PS2".
  PS4      the code SKIPPED RASopathy genes because "their VCEP uses a
           points-based PS4 with phenotype categories, not the flat count
           ladder". GN043's Comments: "PS4: >=5 independent occurrences /
           PS4_Moderate: >=3 / PS4_Supporting: >=1". It is a flat count ladder.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from backend.acmg.hard_coded import (
    _DENOVO_PHENOTYPE_ROW,
    _DENOVO_TABLE1,
    _denovo_points_tier,
    _eval_pm6,
    _eval_ps2,
    _ps4_spec_occurrence_ladder,
)

SPEC = json.loads(
    (Path(__file__).resolve().parents[1] / "data" / "vcep_criteria_spec.json").read_text()
)["genes"]


def test_svi_table1_matches_the_published_document():
    """ClinGen SVI de novo v1.1 Table 1, points per occurrence."""
    assert _DENOVO_TABLE1["highly_specific"] == (2.0, 1.0)
    assert _DENOVO_TABLE1["consistent_not_specific"] == (1.0, 0.5)
    assert _DENOVO_TABLE1["consistent_heterogeneous"] == (0.5, 0.25)
    assert _DENOVO_TABLE1["not_consistent"] == (0.0, 0.0)


def test_only_vceps_that_publish_a_phenotype_row_are_listed():
    """🔴 The guard against the failure that blocked this twice. A VCEP whose
    row is unknown must NOT get a guessed one — the rows are 2x apart."""
    assert set(_DENOVO_PHENOTYPE_ROW) == {"RASopathy", "PotassiumChannel"}
    assert _denovo_points_tier(2, 0, "FBN1", "PS2") is None


def test_ps2_verystrong_at_two_confirmed_occurrences():
    """GN043 Comments, verbatim: "PS2_Very Strong: >=2 independent occurrences
    of PS2". Derived independently from Table 1 (2 points each = 4), so the two
    sources agree."""
    assert _denovo_points_tier(1, 0, "NRAS", "PS2") == "Strong"
    assert _denovo_points_tier(2, 0, "NRAS", "PS2") == "VeryStrong"


def test_ps2_verystrong_from_the_mixed_route():
    """Same Comments: "OR >=2 independent occurrences of PM6 and one occurrence
    of PS2" — 2 + 1 + 1 = 4 points."""
    assert _denovo_points_tier(1, 2, "SOS1", "PS2") == "VeryStrong"


def test_pm6_strong_at_two_assumed_occurrences():
    """GN043: "PM6_Strong: >=2 independent occurrences of PM6"."""
    assert _denovo_points_tier(0, 1, "RAF1", "PM6") == "Moderate"
    assert _denovo_points_tier(0, 2, "RAF1", "PM6") == "Strong"


def test_kcnq1_row_switches_on_pp4():
    """GN112 keys its phenotype row to PP4 (LQT1-specific features)."""
    assert _denovo_points_tier(1, 0, "KCNQ1", "PS2") == "Moderate"
    assert _denovo_points_tier(1, 0, "KCNQ1", "PS2", pp4_met=True) == "Strong"


@pytest.mark.parametrize("gene", ["MYH7", "ACTC1", "MYBPC3", "TNNT2", "BMPR2"])
def test_single_rung_genes_are_skipped_entirely(gene):
    """🔴 The regression guard. These publish ONE PS2 rung (Strong) and ONE PM6
    rung (Moderate). Their text defers to SVI guidance for whether the
    criterion is MET, not for how strong it is, so applying the points ladder
    would DOWNGRADE them: one occurrence is 1 point = Moderate under the
    conservative row, which is 2 points below their published Strong."""
    assert _denovo_points_tier(1, 0, gene, "PS2") is None
    assert _denovo_points_tier(5, 0, gene, "PS2") is None
    assert _eval_ps2(True, gene, 5, 0)["criteria_strength"] == "PS2_Strong"


def test_back_compat_with_no_gene_and_no_counts():
    """An older client that sends neither must behave exactly as before."""
    assert _eval_ps2(True)["criteria_strength"] == "PS2_Strong"
    assert _eval_pm6("trio", "unconfirmed", "AD", False, False
                     )["criteria_strength"] == "PM6_Moderate"


def test_ps2_still_requires_confirmed_parentage():
    assert _eval_ps2(False, "NRAS", 9, 0)["status"] == "not_met"


def test_ps4_ladder_is_read_for_the_rasopathy_genes():
    """The belief this disproves: "their VCEP uses a points-based PS4 with
    phenotype categories, not the flat count ladder"."""
    assert _ps4_spec_occurrence_ladder("PTPN11") == {
        "Supporting": 1, "Moderate": 3, "Strong": 5}


def test_ps4_ladder_is_absent_for_the_cardiomyopathy_genes():
    """Their PS4 is ODDS-RATIO ONLY, verified across all three rungs on
    2026-09-07: STRONG needs the lower bound of the 95% CI around the OR >=20,
    MODERATE >=10, SUPPORTING >=5. No rung offers a proband-count route, so a
    case cohort is required and HeartVar has none.

    FBN1 is excluded for a different reason: its rungs are a POINTS ladder
    ("If >= 4 points." / "If 2-3.5 points." / "If 1-1.5 points."), and the
    points-per-proband mapping is not published.
    """
    for gene in ("MYH7", "ACTC1", "MYBPC3", "TNNT2", "MYL2", "MYL3", "TPM1", "FBN1"):
        assert _ps4_spec_occurrence_ladder(gene) == {}, gene


def test_ps4_ladder_is_read_for_kcnq1_and_bmpr2():
    """These two publish a proband ladder in their STRENGTH ROWS rather than in
    Comments, so the comments-only reader missed them and they fell through to
    the generic Kelly ladder. Verbatim registry text:

      KCNQ1  Supporting "PS4_supporting is Met by 2 independent observations of
                         the variant in probands."
             Moderate   "PS4_moderate is Met by 3-5 probands."
             Strong     "PS4 is Met by 6 or more probands."
      BMPR2  Supporting "Prior observation of the variant in >1 unrelated
                         patients with the same phenotype..."
             Moderate   ">3 unrelated patients"
             Strong     ">4 unrelated patients"

    BMPR2's bounds are STRICT, so >1 is a minimum of 2, >3 of 4, >4 of 5.
    """
    assert _ps4_spec_occurrence_ladder("KCNQ1") == {
        "Supporting": 2, "Moderate": 3, "Strong": 6}
    assert _ps4_spec_occurrence_ladder("BMPR2") == {
        "Supporting": 2, "Moderate": 4, "Strong": 5}


def test_ps4_ladder_never_reads_an_odds_ratio_as_a_proband_count():
    """The Cardiomyopathy rungs carry ">=20" / ">=10" / ">=5" as OR bounds, and
    the row text also contains the worked example "variant detected in 5 out of
    3,500 cases and 1 out of 60,000 controls". A count reader that ignores
    context would turn either into a proband threshold and grade PS4 on it.
    """
    for gene in ("MYH7", "ACTC1", "TNNT2"):
        assert _ps4_spec_occurrence_ladder(gene) == {}, gene


def test_the_ps4_ladder_covers_eighteen_genes_and_no_more():
    """Pin, so a re-harvest that widens or narrows it fails here.

    14 RASopathy genes publish the ladder in Comments as independent
    occurrences; KCNQ1 and BMPR2 publish proband counts in their strength rows;
    LZTR1 and MRAS publish the same >=1/>=3/>=5 rungs as POINTS, at the
    RASopathy VCEP's published one-point-per-case value.

    MYH7 (odds-ratio only, all three rungs) and FBN1 must stay out. FBN1 DOES
    publish a points-per-proband mapping — 1 point for a proband meeting Ghent
    criteria or with ectopia lentis, 0.5 for thoracic aortic disease alone or an
    undescribed phenotype — so the reason is not absence. It is that the mapping
    is PHENOTYPE-CONDITIONAL and HeartVar collects one proband count with no
    per-proband phenotype, so there is no conversion to make. The RASopathy
    value is flat, which is why LZTR1 and MRAS convert and FBN1 does not.
    """
    hits = sorted(g for g in SPEC if _ps4_spec_occurrence_ladder(g))
    assert len(hits) == 18, hits
    assert "PTPN11" in hits and "KCNQ1" in hits and "BMPR2" in hits
    assert "LZTR1" in hits and "MRAS" in hits
    assert "MYH7" not in hits and "FBN1" not in hits
    for gene in ("LZTR1", "MRAS"):
        assert _ps4_spec_occurrence_ladder(gene) == {
            "Supporting": 1, "Moderate": 3, "Strong": 5}, gene


def test_case_control_only_genes_are_identified():
    """The eight Cardiomyopathy genes accept ONLY the odds-ratio route, so a
    proband-count PS4 is not earned on them at any strength. FBN1 (points) and
    the RASopathy genes (count ladder) must not be swept in.
    """
    from backend.acmg.hard_coded import _ps4_case_control_only
    for gene in ("MYH7", "ACTC1", "MYBPC3", "TNNT2", "MYL2", "MYL3", "TPM1"):
        assert _ps4_case_control_only(gene) is True, gene
    for gene in ("PTPN11", "BRAF", "KCNQ1", "BMPR2", "FBN1"):
        assert _ps4_case_control_only(gene) is False, gene
    for bad in ("", "NOTAGENE", None):
        assert _ps4_case_control_only(bad) is False


def test_junk_inputs_never_raise():
    for bad in ("", "x", None, True, -1, 3.7):
        _denovo_points_tier(bad, bad, "NRAS", "PS2")
        _ps4_spec_occurrence_ladder(bad if isinstance(bad, str) else None)


@pytest.mark.parametrize("gene,code", [("SOS1", "PM1"), ("SOS2", "PM1")])
def test_the_comment_level_withdrawals_are_enforced(gene, code):
    """Both are confirmed by the RASopathy specification paper: "No
    well-defined functional domains were observed for SOS1/SOS2, LZTR1, or
    PPP1CB; therefore, PM1 cannot be used for these genes"
    (Wilcox et al. 2025, PMID 40496714)."""
    assert SPEC[gene][code]["applicability"] == "not_applicable"


def test_ppp1cb_pp2_is_applicable_the_withdrawal_was_a_sign_typo():
    """PPP1CB PP2 was withdrawn on a STALE, HIDDEN comment. Reverted 2026-09-08.

    GN128 v1.0/v1.1 carried a sign typo — "Missense z score is <3.09 in gnomAD"
    — corrected to ">3.09" in v1.2 (2024-12-02). The registry hides withdrawn
    content with a `hide` class rather than deleting it, and the harvester read
    Comments without checking visibility, so it picked up a comment the page
    never renders and used it to withdraw a criterion whose Supporting row IS
    rendered: "Missense z score is >3.09 in gnomAD."

    The contradiction was actually noticed when PP2 was pinned as withdrawn on
    2026-09-07 ("against a Supporting row stating the condition as '>3.09'") and
    the wrong side was pinned. Two independent sources say PP2 applies:

      * Wilcox et al. 2025 (PMID 40496714): "Only BRAF, MAP2K1, PTPN11, and
        PPP1CB have missense z-scores greater than 3.09; therefore, PP2 can only
        be applied to these genes."
      * PPP1CB's own eRepo summary: "The missense Z-score is 4.33 which is above
        the threshold set by the Rasopathy VCEP (PP2)."
    """
    assert SPEC["PPP1CB"]["PP2"]["applicability"] == "applicable"
    assert "Supporting" in (SPEC["PPP1CB"]["PP2"].get("strengths") or {})
    assert "<3.09" in (SPEC["PPP1CB"]["PP2"].get("comments") or "")


@pytest.mark.parametrize("gene,code", [("BMPR2", "BS3"), ("LZTR1", "PM3"),
                                       ("LZTR1", "PVS1")])
def test_the_near_misses_are_left_applicable(gene, code):
    """🔴 The scope boundary. BMPR2's BS3 is "not applicable FOR SPLICING
    EFFECTS" — one assay class, not the criterion. LZTR1's two are bare "Not
    applicable." against specs saying the rule applies "when curating for AR
    disease only" — conditional on inheritance, not withdrawn."""
    assert SPEC[gene][code]["applicability"] == "applicable"
