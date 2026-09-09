"""BA1 / BS1 thresholds are INCLUSIVE at the floor, per every spec sampled.

THE BUG. Both criteria used strict ``>``:

    if gnomad_af is not None and gnomad_af > ba1_threshold:   # BA1
    if gnomad_af > bs1_threshold and gnomad_af <= ba1_threshold:  # BS1

Every specification states these bounds with ">=" — MYH7 GN002 ">=0.001",
KCNQ1 GN112 "greater than or equal to 0.004", LDLR GN013 ">=0.005". So a variant
sitting EXACTLY on a published threshold escaped both criteria: a benign-direction
under-call, which in a tool with documented pathogenic over-call problems is the
wrong way to be wrong.

The BA1 ceiling on BS1 had to flip to strict "<" in the same change. With BA1
firing at ">= ba1_threshold", leaving BS1 at "<= ba1_threshold" would let a
variant exactly on the BA1 cutoff satisfy BOTH and double-count one frequency
observation. The window is [bs1, ba1) and the two bounds are one decision.
"""
from __future__ import annotations

from backend.acmg.hard_coded import _eval_ba1, _eval_bs1

BA1 = 0.001
BS1 = 0.0001
AN = 50_000


def _ba1(af, an=AN):
    return _eval_ba1(af, BA1, "MYH7 GN002", gnomad_an=an)


def _bs1(af, ba1_met=False):
    return _eval_bs1(af, BS1, BA1, "AD", "MYH7 GN002", ba1_met=ba1_met)


def test_ba1_fires_exactly_on_the_threshold():
    out = _ba1(BA1)
    assert out["status"] == "met", (
        "MYH7 GN002 says '>=0.001'; an AF of exactly 0.001 must fire BA1"
    )
    assert out["criteria_strength"] == "BA1"


def test_ba1_does_not_fire_just_below_the_threshold():
    assert _ba1(BA1 * 0.999)["status"] != "met"


def test_ba1_still_needs_the_two_thousand_allele_floor_on_the_boundary():
    """The inclusive bound must not smuggle a low-AN site past the SVI gate."""
    out = _ba1(BA1, an=1999)
    assert out["status"] == "not_met"
    assert "2,000" in out["evidence"]


def test_bs1_fires_exactly_on_its_floor():
    out = _bs1(BS1)
    assert out["status"] == "met"
    assert out["criteria_strength"] == "BS1_Strong"


def test_bs1_does_not_fire_just_below_its_floor():
    assert _bs1(BS1 * 0.999)["status"] != "met"


def test_bs1_does_not_fire_exactly_on_the_ba1_ceiling():
    """This is the double-count guard: BA1 owns that value now."""
    assert _bs1(BA1)["status"] != "met", (
        "an AF exactly on the BA1 cutoff fires BA1; BS1 must not also fire or "
        "one frequency observation is counted twice"
    )


def test_bs1_fires_just_below_the_ba1_ceiling():
    out = _bs1(BA1 * 0.999)
    assert out["status"] == "met", "the window is [bs1, ba1), so just under counts"


def test_bs1_and_ba1_never_both_fire_across_the_range():
    """Sweep the interesting neighbourhood — the two criteria must partition it."""
    for af in (BS1 * 0.5, BS1, BS1 * 2, BA1 * 0.5, BA1 * 0.999, BA1, BA1 * 2):
        ba1_met = _ba1(af)["status"] == "met"
        bs1_met = _bs1(af, ba1_met=ba1_met)["status"] == "met"
        assert not (ba1_met and bs1_met), f"both fired at af={af}"
