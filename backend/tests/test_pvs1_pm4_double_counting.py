"""PM4 must not be applied alongside PVS1 at any strength.

Abou Tayoun 2018 states it as a requirement, not a suggestion:

    "at the Moderate strength level, there is potential overlap in usage of
    PVS1_Moderate and PM4 (protein length changing variant). To prevent
    double-counting of this evidence type, we recommend that PM4 should not be
    applied for any variant in which PVS1, at any strength level, is also
    applied."

THE GAP. `apply_cross_criterion_exclusions` had a PM4-vs-BP3 guard and a
PVS1-demotes-PP3 guard, but nothing stopped PM4 and PVS1 co-firing. A whole-exon
in-frame deletion can plausibly earn `exon_loss_variant` (→ bare PVS1, +8) and a
protein-length change (→ PM4_Moderate, +2) from the SAME observation, for +10 off
one fact. Note "at any strength level" — restricting the guard to PVS1_Moderate,
where the paper motivates it, would miss exactly that worst case.
"""
from __future__ import annotations

from backend.acmg.hard_coded import apply_cross_criterion_exclusions


def _met(code, strength):
    return {"code": code, "status": "met", "criteria_strength": strength,
            "evidence": f"{code} test fixture"}


def _run(*entries):
    criteria, _forced = apply_cross_criterion_exclusions(list(entries))
    return {c["code"]: c for c in criteria}


def test_pm4_is_demoted_when_pvs1_is_met_very_strong():
    out = _run(_met("PVS1", "PVS1"), _met("PM4", "PM4_Moderate"))
    assert out["PVS1"]["status"] == "met", "PVS1 is the more specific criterion"
    assert out["PM4"]["status"] != "met", (
        "PM4 alongside bare PVS1 double-counts one protein-length observation "
        "for +10"
    )


def test_pm4_is_demoted_at_every_pvs1_strength():
    """'any strength level', verbatim — including the tiers below Moderate."""
    for strength in ("PVS1", "PVS1_Strong", "PVS1_Moderate", "PVS1_Supporting"):
        out = _run(_met("PVS1", strength), _met("PM4", "PM4_Moderate"))
        assert out["PM4"]["status"] != "met", strength


def test_the_demotion_names_pvs1_and_its_strength():
    out = _run(_met("PVS1", "PVS1_Moderate"), _met("PM4", "PM4_Moderate"))
    ev = out["PM4"]["evidence"]
    assert "PVS1_Moderate" in ev
    assert "double-count" in ev.lower()


def test_pm4_survives_when_pvs1_is_not_met():
    out = _run(
        {"code": "PVS1", "status": "not_met", "criteria_strength": None,
         "evidence": "not a null variant"},
        _met("PM4", "PM4_Moderate"),
    )
    assert out["PM4"]["status"] == "met", (
        "the guard must key on PVS1 being MET, not merely present"
    )


def test_pm4_survives_on_its_own():
    out = _run(_met("PM4", "PM4_Moderate"))
    assert out["PM4"]["status"] == "met"


def test_pvs1_is_never_the_one_demoted():
    """Direction matters: yielding PVS1 to PM4 would drop 8 points to 2."""
    out = _run(_met("PVS1", "PVS1"), _met("PM4", "PM4_Moderate"))
    assert out["PVS1"]["criteria_strength"] == "PVS1"
