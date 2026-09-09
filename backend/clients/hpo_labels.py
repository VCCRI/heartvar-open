"""HPO term-id → primary label lookup.

Loads ``backend/data/hpo_labels.json`` (built by ``build_hpo_labels.py``
from the official Human Phenotype Ontology release) on first call and
caches the dict for the process lifetime. Callers use
:func:`resolve_hpo_label` to turn ``HP:NNNNNNN`` into a human-readable
phenotype name; :func:`resolve_hpo_tokens` is the convenience wrapper
that parses a comma-separated curator input string and resolves each
token (HPO IDs map to labels, free-text tokens pass through unchanged).
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

_LABELS_PATH = Path(__file__).resolve().parents[1] / "data" / "hpo_labels.json"
_HPO_ID_RE = re.compile(r"^HP:\d{7}$")


@lru_cache(maxsize=1)
def _load_labels() -> dict[str, str]:
    """Return the HP-id → label map, or an empty dict if the file is
    missing (build script not yet run). Callers degrade gracefully —
    missing labels just mean the ID passes through unchanged."""
    if not _LABELS_PATH.exists():
        return {}
    try:
        with _LABELS_PATH.open("r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return {}


def resolve_hpo_label(hpo_id: str) -> str | None:
    """Return the canonical label for ``hpo_id`` or ``None`` if not found."""
    if not hpo_id:
        return None
    key = hpo_id.strip().upper()
    if not _HPO_ID_RE.match(key):
        return None
    return _load_labels().get(key)


def resolved_token_map(raw: str) -> dict[str, str]:
    """Return ``{original_token: resolved_label}`` for every comma-separated
    token in ``raw``. Used by the curate endpoint to ship resolved labels
    alongside the original input so the frontend can build PubMed search
    URLs without needing its own HPO table.
    """
    out: dict[str, str] = {}
    for tok in (raw or "").split(","):
        tok = tok.strip()
        if not tok:
            continue
        if _HPO_ID_RE.match(tok.upper()):
            label = resolve_hpo_label(tok)
            if label:
                out[tok] = label
        else:
            out[tok] = tok
    return out
