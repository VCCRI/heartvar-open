"""PM3 phase weighting, and the SpliceAI PP3 delta.

Both paths were changed on 2026-08-25 and NEITHER had a test — the full suite
passed unchanged through a PM3_Moderate -> PM3_Supporting demotion and a
0.5 -> 0.2 PP3 threshold. That silence is the reason this file exists.

PM3 — ClinGen SVI "Recommendation for in trans Criterion (PM3)" v1.0 (approved
2019-05-02), Table 1, points per in-trans proband:

                                     confirmed in trans | phase unknown
    Pathogenic or Likely pathogenic          1.0        |     0.5
    Homozygous occurrence (max 1.0)     0.5 P / 0.25 LP |     N/A
    Uncertain significance (max 0.5)         0.25       |     0.0

Table 2, total points -> strength: PM3_Supporting 0.5, PM3 1.0, PM3_Strong 2.0,
PM3_VeryStrong 4.0. So a phase-unknown observation is HALF a confirmed one, and
"inherited from an affected parent" is phase-unknown by the branch's own
admission.

SpliceAI — KCNQ1 GN112 says ">= 0.2"; ACGS 2024 v1.2 calls ">0.2" the most
accurate single splice threshold. HeartVar sat at 0.5 while its OWN BP4 arm and
PVS1 splice routing already used 0.2, so one delta could be "possible splice
impact" for one purpose and not for another.
"""
from __future__ import annotations

from backend.acmg.hard_coded import _SPLICEAI_PP3_DELTA, _eval_pm3
from backend.acmg.tiers import _points_for

GENE = "SOMEGENE"


def _pm3(in_trans="", denovo=""):
    return _eval_pm3(GENE, "AR", "het", in_trans, denovo)


def test_confirmed_in_trans_is_moderate():
    out = _pm3(in_trans="yes")
    assert out["status"] == "met"
    assert out["criteria_strength"] == "PM3_Moderate", (
        "1.0 point (Table 1, confirmed in trans) = PM3 at Moderate (Table 2)"
    )


def test_phase_unknown_is_supporting_not_moderate():
    out = _pm3(denovo="inherited_affected")
    assert out["status"] == "met"
    assert out["criteria_strength"] == "PM3_Supporting", (
        "inherited-from-affected-parent is a POSSIBLE trans configuration, so "
        "it sits in Table 1's phase-unknown column at 0.5 points = "
        "PM3_Supporting; Moderate is reserved for confirmed phase"
    )


def test_phase_unknown_is_worth_exactly_half_of_confirmed():
    """The point ratio in Table 1, asserted through the scoring model."""
    confirmed = _points_for(_pm3(in_trans="yes")["criteria_strength"], "PM3")
    unknown = _points_for(_pm3(denovo="inherited_affected")["criteria_strength"], "PM3")
    assert (confirmed, unknown) == (2, 1)


def test_pm3_not_met_without_either_observation():
    assert _pm3()["status"] == "not_met"


def test_pm3_needs_recessive_and_heterozygous():
    assert _eval_pm3(GENE, "AD", "het", "yes", "")["status"] == "not_met"
    assert _eval_pm3(GENE, "AR", "hom", "yes", "")["status"] == "not_met"


def test_pm3_respects_a_cspec_not_applicable_gene():
    """The cardiomyopathy / RASopathy / FBN1 VCEPs mark PM3 Not Applicable."""
    out = _eval_pm3("MYH7", "AR", "het", "yes", "")
    assert out["status"] == "not_met"
    assert "Not Applicable" in out["evidence"]


def test_spliceai_pp3_delta_is_the_spec_value():
    assert _SPLICEAI_PP3_DELTA == 0.2, (
        "KCNQ1 GN112 says '>= 0.2' and ACGS 2024 says '>0.2'; nothing HeartVar "
        "cites supports 0.5"
    )


def test_spliceai_pp3_delta_matches_the_other_arms_that_already_used_0_2():
    """The BP4 arm and the PVS1 splice routing both gate on 0.2. A PP3 threshold
    above them let the same delta mean two different things in one file."""
    import inspect

    from backend.acmg import hard_coded

    src = inspect.getsource(hard_coded)
    assert "spliceai_max >= 0.2" in src
    assert "spliceai_max >= _SPLICEAI_PP3_DELTA" in src
