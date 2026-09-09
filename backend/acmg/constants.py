"""ACMG scoring constants — SINGLE SOURCE OF TRUTH.

The Tavtigian 2020 point values, criterion names, hybrid-classification split,
classification thresholds and the canonical 28-code display order live in
``static/acmg_constants.json`` so the same numbers drive BOTH the server and the
browser (``static/heartvar.js`` fetches the same file). Editing the values in
code is a mistake — edit the JSON. A cross-language parity test
(``backend/tests/test_acmg_constants_parity.py``) fails CI if the two sides
drift.

This module computes the JSON path from its own location (parents[2] = repo
root) so it never imports ``backend.app``.
"""
from __future__ import annotations

import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

_ACMG_CONSTANTS_FILE = PROJECT_ROOT / "static" / "acmg_constants.json"
try:
    _ACMG_CONSTANTS = json.loads(_ACMG_CONSTANTS_FILE.read_text())
except (OSError, ValueError) as _e:
    raise RuntimeError(
        f"ACMG constants file is missing or malformed ({_ACMG_CONSTANTS_FILE}): {_e}. "
        "This file is required — it is the single source of truth for ACMG "
        "scoring shared with static/heartvar.js."
    ) from _e


_TIER_POINTS: dict[str, int] = dict(_ACMG_CONSTANTS["tier_points"])

_BARE_CODE_POINTS: dict[str, int] = dict(_ACMG_CONSTANTS["bare_code_points"])

_CLASSIFICATION_THRESHOLDS: list[dict] = list(_ACMG_CONSTANTS["classification_thresholds"])


HARD_CODED_CRITERIA_CODES: tuple[str, ...] = tuple(
    _ACMG_CONSTANTS["hard_coded_criteria_codes"]
)

AI_EVALUATED_CRITERIA_CODES: tuple[str, ...] = tuple(
    _ACMG_CONSTANTS["ai_evaluated_criteria_codes"]
)


_CRITERION_NAMES: dict[str, str] = dict(_ACMG_CONSTANTS["criterion_names"])


CANONICAL_CRITERIA_ORDER: tuple[str, ...] = tuple(
    _ACMG_CONSTANTS["criteria_display_order"]
)


_CLINVAR_STAR_STRENGTH = {4: "Strong", 3: "Strong", 2: "Moderate", 1: "Supporting"}
