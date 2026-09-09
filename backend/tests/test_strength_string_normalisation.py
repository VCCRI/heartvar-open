"""_points_for tolerates non-canonical criteria_strength spellings.

Everything HeartVar generates is canonical ("PS2_VeryStrong"), but the AI
layer's criteria_strength is free text — the prompt asks for a bare code and
nothing enforces it. An unrecognised strength used to score ZERO, silently, in
the under-call direction. These pin the normalisation and, just as important,
pin that it does NOT invent points for a tier that genuinely isn't one.
"""
import pytest

from backend.acmg.tiers import _points_for, compute_points_total


@pytest.mark.parametrize("strength", [
    "PS2_VeryStrong",
    "PS2_Very Strong",
    "PS2_Very_Strong",
    "PS2_very strong",
    "PS2_VERYSTRONG",
])
def test_very_strong_spellings_all_score_eight(strength):
    assert _points_for(strength, "PS2") == 8


def test_lowercase_tier_scores():
    assert _points_for("PS3_strong", "PS3") == 4
    assert _points_for("PM1_moderate", "PM1") == 2
    assert _points_for("PP1_supporting", "PP1") == 1


def test_benign_direction_survives_normalisation():
    assert _points_for("BS1_Strong", "BS1") == -4
    assert _points_for("BS3_very strong", "BS3") == -8
    assert _points_for("bp4_supporting", "BP4") == -1


def test_bare_codes_are_case_insensitive():
    assert _points_for("PS3", "PS3") == 4
    assert _points_for("ps3", "ps3") == 4
    assert _points_for("BA1", "BA1") == -8


def test_unrecognised_tier_still_scores_zero():
    """The normaliser must not manufacture points from a token that is not a
    strength. "PS4_Nonsense" is not "PS4 at Strong"."""
    assert _points_for("PS4_Nonsense", "PS4") == 0
    assert _points_for("PS4_Definitive", "PS4") == 0
    assert _points_for("Wingding", "PS4") == 0


def test_pp5_bp6_zeroing_is_case_insensitive():
    """The PP5/BP6 circularity guard keys off `code`. A lowercase code must not
    slip a reputable-source assertion into the point total."""
    assert _points_for("PP5", "PP5") == 0
    assert _points_for("PP5", "pp5") == 0
    assert _points_for("BP6", "bp6") == 0
    assert _points_for("PP5", "pp5", True) == 1


def test_total_over_mixed_spellings():
    crit = [
        {"code": "PS2", "status": "met", "criteria_strength": "PS2_Very Strong"},
        {"code": "PM2", "status": "met", "criteria_strength": "PM2_Supporting"},
        {"code": "BP4", "status": "not_met", "criteria_strength": "BP4"},
    ]
    assert compute_points_total(crit) == 9
