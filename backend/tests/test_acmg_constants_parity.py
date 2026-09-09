"""Cross-language parity guard for the ACMG scoring constants.

The Tavtigian-2020 point values, criterion names, hybrid-classification split
and classification thresholds are maintained in ONE place —
``static/acmg_constants.json`` — and consumed by BOTH the Python backend
(``backend.app``, loaded at import) and the browser (``static/heartvar.js``,
fetched at runtime). This test fails CI if the two sides ever drift:

  1. The Python module-level constants equal the JSON (and the JSON equals the
     known-correct Tavtigian values, including ``VeryStrong: 8``).
  2. Completeness — every criterion code has a name AND a point value.
  3. ``classification_for`` matches the JSON threshold table across a battery
     of boundary point totals.
  4. A static guard over ``static/heartvar.js``: it must NOT re-introduce a
     hardcoded tier-points object literal (the old
     ``_TIER_POINTS = { Strong: 4, ... }`` shape) — it must derive the tables
     from the fetched JSON instead.

No pytest dependency — runnable with pytest or directly
(``python -m backend.tests.test_acmg_constants_parity``).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from backend.app import (
    _TIER_POINTS,
    _BARE_CODE_POINTS,
    _CRITERION_NAMES,
    classification_for,
    HARD_CODED_CRITERIA_CODES,
    AI_EVALUATED_CRITERIA_CODES,
    CANONICAL_CRITERIA_ORDER,
)

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
_JSON_PATH = _PROJECT_ROOT / "static" / "acmg_constants.json"
_JS_PATH = _PROJECT_ROOT / "static" / "heartvar.js"


def _load_json() -> dict:
    return json.loads(_JSON_PATH.read_text())


def test_tier_points_match_json_and_known_correct():
    d = _load_json()
    assert _TIER_POINTS == d["tier_points"]
    assert _TIER_POINTS == {
        "VeryStrong": 8,
        "Strong": 4,
        "Moderate": 2,
        "Supporting": 1,
    }


def test_bare_code_points_match_json():
    d = _load_json()
    assert _BARE_CODE_POINTS == d["bare_code_points"]


def test_criterion_names_match_json():
    d = _load_json()
    assert _CRITERION_NAMES == d["criterion_names"]


def test_hard_coded_codes_match_json_with_order():
    d = _load_json()
    json_codes = d["hard_coded_criteria_codes"]
    assert set(HARD_CODED_CRITERIA_CODES) == set(json_codes)
    assert list(HARD_CODED_CRITERIA_CODES) == list(json_codes)


def test_ai_codes_match_json_with_order():
    d = _load_json()
    json_codes = d["ai_evaluated_criteria_codes"]
    assert set(AI_EVALUATED_CRITERIA_CODES) == set(json_codes)
    assert list(AI_EVALUATED_CRITERIA_CODES) == list(json_codes)


def test_canonical_order_matches_json():
    d = _load_json()
    assert list(CANONICAL_CRITERIA_ORDER) == list(d["criteria_display_order"])


def test_every_code_has_a_name_and_a_point_value():
    all_codes = set(HARD_CODED_CRITERIA_CODES) | set(AI_EVALUATED_CRITERIA_CODES)
    assert len(all_codes) == 28
    assert all_codes == set(CANONICAL_CRITERIA_ORDER)
    for code in all_codes:
        assert code in _CRITERION_NAMES, f"{code} missing from criterion_names"
        assert code in _BARE_CODE_POINTS, f"{code} missing from bare_code_points"


def test_hard_coded_and_ai_are_disjoint():
    assert not (set(HARD_CODED_CRITERIA_CODES) & set(AI_EVALUATED_CRITERIA_CODES))


def test_the_hybrid_split_is_21_rule_based_and_7_llm():
    """PINS the 2026-09-08 split: 18 + 10 became 21 + 7.

    PM1, PP1 and BS4 moved from the LLM side to the Python side, joining PS1 and
    PM5, which the merge had already been discarding. The LLM's verdicts on the
    three matched the Python derivation in substance, differing mostly by
    emitting a bare code where Python emits an equal-valued strength suffix (PM1
    against PM1_Moderate, both +2). PM1 and PP1 were also unstable across two
    identical AI runs, so the model was adding run-to-run variance without
    adding information.

    Handing any of the five back to the LLM has to change these two numbers,
    which is why they are pinned rather than derived.
    """
    assert len(HARD_CODED_CRITERIA_CODES) == 21, (
        "the rule-based half is no longer 21 codes"
    )
    assert len(AI_EVALUATED_CRITERIA_CODES) == 7, (
        "the LLM half is no longer 7 codes"
    )
    assert set(AI_EVALUATED_CRITERIA_CODES) == {
        "PS3", "PS4", "PP2", "PP4", "BS3", "BP1", "BP5",
    }
    for code in ("PM1", "PP1", "BS4", "PS1", "PM5"):
        assert code in HARD_CODED_CRITERIA_CODES, (
            f"{code} became Python-authoritative on 2026-09-08 and must stay "
            "on the hard-coded side"
        )


def _expected_label(points: int, thresholds: list[dict]) -> str:
    for band in thresholds:
        if points >= band["min"]:
            return band["label"]
    raise AssertionError("threshold table has no catch-all band")


def test_classification_for_matches_threshold_table():
    thresholds = _load_json()["classification_thresholds"]
    for p in (-20, -7, -6, -1, 0, 5, 6, 9, 10, 15):
        assert classification_for(p) == _expected_label(p, thresholds), (
            f"classification_for({p}) drifted from the JSON threshold table"
        )


def test_classification_for_boundaries_unchanged():
    assert classification_for(10) == "Pathogenic"
    assert classification_for(9) == "Likely pathogenic"
    assert classification_for(6) == "Likely pathogenic"
    assert classification_for(5) == "VUS"
    assert classification_for(0) == "VUS"
    assert classification_for(-1) == "Likely benign"
    assert classification_for(-6) == "Likely benign"
    assert classification_for(-7) == "Benign"


def test_js_does_not_hardcode_tier_points_literal():
    """The JS must derive tier points from the fetched JSON, not a literal.

    Matches the OLD shape ``{ ... Strong: 4 ... }`` / ``{ ... VeryStrong: 8
    ... }`` — an object literal pairing a tier name directly with its integer.
    The populate-from-fetch assignment ``_TIER_POINTS = d.tier_points`` has no
    such inline number, so it does not trip this guard.
    """
    src = _JS_PATH.read_text()
    bad = re.compile(
        r"['\"]?(?:VeryStrong|Strong|Moderate|Supporting)['\"]?\s*:\s*-?\d+"
    )
    hits = [m.group(0) for m in bad.finditer(src)]
    assert not hits, (
        "static/heartvar.js appears to hardcode a tier-points literal "
        f"(found {hits!r}); tier points must come from the fetched "
        "static/acmg_constants.json instead."
    )


def test_js_does_not_hardcode_bare_code_points_literal():
    """No re-introduced bare-code points table (the old ``PVS1: 8`` shape)."""
    src = _JS_PATH.read_text()
    assert not re.search(r"PVS1['\"]?\s*:\s*8", src), (
        "static/heartvar.js appears to hardcode the bare-code points table "
        "(PVS1: 8); it must come from the fetched static/acmg_constants.json."
    )


def test_js_has_no_second_classifier():
    """The browser must NOT carry its own scorer or tier classifier.

    The server owns both (``backend.acmg.tiers`` — compute_points_total and
    classification_for_criteria, including the ACMG-2015 benign combining
    floor) and the page renders the ``points_total``/``classification`` it is
    handed. The old JS twins were remnants of the removed BYOK client-AI path
    and were deleted 2026-08-17: an unreferenced second classifier cannot be
    kept honest, and this one had already gone stale (it still returned
    "Likely benign" for a lone BP4 after the backend stopped doing so).

    ``_pointsFor`` is gone too; ``CRIT_POINTS`` remains for DISPLAY-only
    per-criterion point chips and is not a scorer."""
    src = _JS_PATH.read_text()
    for banned in (
        "function classificationFor",
        "function computePointsTotal",
        "function _pointsFor",
    ):
        assert banned not in src, (
            f"static/heartvar.js re-introduces `{banned}` — the tier and the "
            "point total are the server's (backend/acmg/tiers.py). A browser "
            "twin will drift; call the server instead."
        )


if __name__ == "__main__":  # pragma: no cover - manual run convenience
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("All ACMG-constants parity checks passed.")
