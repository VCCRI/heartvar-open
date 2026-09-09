"""ACMG/AMP-2015 combining-rule floor on the benign side.

THE DEFECT. The Tavtigian-2020 point bands the engine scores against place
"Likely benign" at -1 to -6. A single Supporting-strength benign criterion is
worth exactly -1, so ONE weak benign line of evidence — in practice almost
always BP4 off an in-silico predictor — was enough to move a variant out of
VUS and call it Likely benign.

ACMG/AMP 2015 (Richards et al., Genet Med 2015, Table 5) does not allow that.
Its benign combining rules are:

    Benign         (i) BA1 stand-alone, or (ii) >= 2 Strong (BS1-BS4)
    Likely benign  (i) 1 Strong + 1 Supporting, or (ii) >= 2 Supporting (BP1-BP7)

A lone BP4 satisfies neither, so under the 2015 rules it is a VUS. This is the
best-known divergence between the verbal rules and the point system, and it
runs in the false-benign direction, which is the one that matters clinically.

REAL CASE that prompted this (2026-08-17): CHD7 NM_017780.4:c.2098A>G
p.(Asn700Asp) (8-60794987-A-G, hg38). ClinVar VCV001359827 is Uncertain
significance at 2 stars, no conflicts. HeartVar returned Likely benign off
BP4_Supporting alone (REVEL 0.168). Frequency contributed nothing in either
direction: gnomAD v4 FAF95 popmax 6.96e-05 sits in the dead zone between CHD7's
constraint-aware PM2 ceiling (4e-05) and BS1 floor (1e-04). The same signature
is in the stored eRepo set-2 benchmark as NRAS c.250A>G, where the expert panel
also said VUS.

WHAT THE FLOOR DOES. A tier of "Likely benign" is pushed back to "VUS" when the
met benign evidence amounts to a single criterion weaker than Strong.

WHAT IT DELIBERATELY DOES NOT DO. A lone Strong benign criterion (BS1 alone,
-4) still reads Likely benign. Strict 2015 would call that a VUS too (it wants
1 Strong AND 1 Supporting), but a strong frequency argument is genuinely good
benign evidence, that case is not the reported defect, and demoting it is a
much larger calibration shift. Kept on the point-system behaviour on purpose.

No pytest dependency — runnable with pytest or directly
(``python -m backend.tests.test_benign_combining_floor``).
"""
from __future__ import annotations

from backend.app import (
    _CRITERION_NAMES,
    apply_benign_combining_floor,
    classification_for,
    classification_for_criteria,
    compute_points_total,
    tier_without_clinvar_assertion,
)


def _crit(code, status="met", strength=None):
    """A criterion dict in the shape the aggregation functions expect. A met
    criterion defaults its strength to the bare code; not_met → None."""
    if strength is None:
        strength = code if status == "met" else None
    return {
        "code": code,
        "name": _CRITERION_NAMES.get(code, code),
        "status": status,
        "direction": "benign" if code.startswith("B") else "pathogenic",
        "criteria_strength": strength,
        "evidence": "",
        "source": "hard_coded",
    }


def test_lone_supporting_benign_criterion_stays_vus():
    """CHD7 c.2098A>G: BP4_Supporting alone must not reach Likely benign."""
    criteria = [_crit("BP4", strength="BP4_Supporting")]
    assert compute_points_total(criteria) == -1
    assert classification_for(-1) == "Likely benign"
    assert classification_for_criteria(-1, criteria) == "VUS"


def test_lone_moderate_benign_criterion_stays_vus():
    """One criterion is one criterion — Moderate strength does not rescue it."""
    criteria = [_crit("BP4", strength="BP4_Moderate")]
    assert compute_points_total(criteria) == -2
    assert classification_for_criteria(-2, criteria) == "VUS"


def test_two_supporting_benign_criteria_reach_likely_benign():
    """>= 2 Supporting is the 2015 Table 5 route to Likely benign."""
    criteria = [
        _crit("BP4", strength="BP4_Supporting"),
        _crit("BP1", strength="BP1_Supporting"),
    ]
    assert compute_points_total(criteria) == -2
    assert classification_for_criteria(-2, criteria) == "Likely benign"


def test_lone_strong_benign_criterion_still_likely_benign():
    """Deliberate carve-out: a lone BS1 keeps its point-system tier."""
    criteria = [_crit("BS1", strength="BS1_Strong")]
    assert compute_points_total(criteria) == -4
    assert classification_for_criteria(-4, criteria) == "Likely benign"


def test_strong_plus_supporting_reaches_likely_benign():
    criteria = [
        _crit("BS1", strength="BS1_Strong"),
        _crit("BP4", strength="BP4_Supporting"),
    ]
    assert classification_for_criteria(-5, criteria) == "Likely benign"


def test_two_strong_benign_criteria_still_benign():
    """The floor only intercepts Likely benign — Benign is untouched."""
    criteria = [
        _crit("BS1", strength="BS1_Strong"),
        _crit("BS2", strength="BS2_Strong"),
    ]
    assert compute_points_total(criteria) == -8
    assert classification_for_criteria(-8, criteria) == "Benign"


def test_bp6_does_not_count_as_a_second_benign_criterion():
    """BP6 is retired (SVI 2018) and scores 0, so BP4+BP6 is still one line of
    scoring benign evidence — it must not buy Likely benign."""
    criteria = [
        _crit("BP4", strength="BP4_Supporting"),
        _crit("BP6", strength="BP6"),
    ]
    assert compute_points_total(criteria) == -1
    assert classification_for_criteria(-1, criteria) == "VUS"


def test_not_met_benign_criteria_do_not_count():
    criteria = [
        _crit("BP4", strength="BP4_Supporting"),
        _crit("BP1", status="not_met"),
        _crit("BP7", status="not_assessed"),
    ]
    assert classification_for_criteria(-1, criteria) == "VUS"


def test_pathogenic_tiers_are_untouched():
    criteria = [_crit("PP3", strength="PP3_Supporting")]
    assert classification_for_criteria(1, criteria) == "VUS"
    assert classification_for_criteria(6, criteria) == "Likely pathogenic"
    assert classification_for_criteria(10, criteria) == "Pathogenic"


def test_forced_classification_wins_over_the_floor():
    """BA1 forces Benign in apply_cross_criterion_exclusions; the floor must
    not be able to override a stand-alone call."""
    criteria = [_crit("BA1", strength="BA1")]
    assert classification_for_criteria(-8, criteria, "Benign") == "Benign"


def test_floor_applies_to_the_no_pp5_proxy_tier():
    """tier_without_clinvar_assertion is the novel-variant proxy reported next
    to the headline tier — it must not disagree with it."""
    criteria = [
        _crit("BP4", strength="BP4_Supporting"),
        _crit("BP6", strength="BP6"),
    ]
    pts, tier = tier_without_clinvar_assertion(criteria)
    assert pts == -1
    assert tier == "VUS"


def test_floor_helper_only_rewrites_likely_benign():
    lone = [_crit("BP4", strength="BP4_Supporting")]
    assert apply_benign_combining_floor("Likely benign", lone) == "VUS"
    assert apply_benign_combining_floor("VUS", lone) == "VUS"
    assert apply_benign_combining_floor("Benign", lone) == "Benign"
    assert apply_benign_combining_floor("Pathogenic", lone) == "Pathogenic"


def test_floor_helper_tolerates_an_empty_criteria_list():
    assert apply_benign_combining_floor("Likely benign", []) == "Likely benign"
    assert apply_benign_combining_floor("Likely benign", None) == "Likely benign"


if __name__ == "__main__":
    import sys

    mod = sys.modules[__name__]
    failures = 0
    for _name in [n for n in dir(mod) if n.startswith("test_")]:
        try:
            getattr(mod, _name)()
            print(f"PASS {_name}")
        except AssertionError as exc:
            failures += 1
            print(f"FAIL {_name}: {exc}")
    print(f"\n{failures} failure(s)")
    sys.exit(1 if failures else 0)
