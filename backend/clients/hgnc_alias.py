"""HGNC gene-symbol alias resolver (LOCAL-FIRST, synchronous).

Maps an outdated / alias / previous gene symbol a curator types to the
current HGNC-approved symbol, so a stale symbol (e.g. ``CSX`` for
``NKX2-5``) doesn't silently produce empty evidence by sending the wrong
gene to the downstream evidence clients.

PREFERS the local map ``backend/data/hgnc_alias_map.json`` (built by
``scripts/build_hgnc_alias_db.py`` from HGNC's public complete set) and is
PURE / LOCAL — it never touches the network on the hot path. The map is
loaded once and cached for the life of the process; the ``HGNC_ALIAS_MAP_PATH``
env var overrides its location.

The map JSON has the shape::

    {"approved": ["A1BG", "A2M", ...],            # sorted approved symbols (UPPER)
     "aliases": {"ALIASUPPER": "ApprovedSymbol"}}  # alias/prev (UPPER) -> approved

``canonicalise_gene_symbol(gene)`` returns::

    {"input": <original string>,
     "approved": <approved symbol or None>,
     "is_alias": <bool>,        # True only when input matched via the alias map
     "recognized": <bool>,      # True if input is a known approved symbol or alias
     "ambiguous": <bool>}       # reserved; always False on the local path

FAIL-OPEN: if the local map is absent or unreadable, this returns a NEUTRAL
result (``recognized=True``, ``approved=gene`` unchanged, ``is_alias=False``)
and logs a warning ONCE — a missing optional alias map must never block a
curation. (Ambiguous aliases that map to >1 approved symbol are dropped at
BUILD time and so never appear here; ``ambiguous`` is kept in the return
shape for a possible future live HGNC fallback that could surface them.)

A live HGNC fallback (rest.genenames.org) is intentionally NOT added here to
keep this module pure/local; it can be a future addition layered above this.

Rebuild the map periodically::

    python3 scripts/build_hgnc_alias_db.py
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any

log = logging.getLogger("heartvar.hgnc_alias")

_BACKEND_ROOT = Path(__file__).resolve().parent.parent
MAP_PATH = Path(
    os.environ.get("HGNC_ALIAS_MAP_PATH")
    or (_BACKEND_ROOT / "data" / "hgnc_alias_map.json")
)

_APPROVED_UPPER_TO_CANONICAL: dict[str, str] | None = None
_ALIASES: dict[str, str] | None = None
_LOADED = False
_WARNED_MISSING = False


def _load_map() -> bool:
    """Lazily load the local alias map into the process caches. Returns True
    on success, False when the map is absent/unreadable (fail-open). Loads at
    most once; subsequent calls are no-ops once ``_LOADED`` is set."""
    global _APPROVED_UPPER_TO_CANONICAL, _ALIASES, _LOADED, _WARNED_MISSING
    if _LOADED:
        return _APPROVED_UPPER_TO_CANONICAL is not None
    _LOADED = True
    path = MAP_PATH
    if not path.exists():
        if not _WARNED_MISSING:
            log.warning(
                "HGNC alias map not found at %s — gene-symbol canonicalisation "
                "disabled (fail-open). Build it with "
                "scripts/build_hgnc_alias_db.py.",
                path,
            )
            _WARNED_MISSING = True
        return False
    try:
        with path.open("r", encoding="utf-8") as f:
            data: dict[str, Any] = json.load(f)
        approved = data.get("approved") or []
        aliases = data.get("aliases") or {}
        _APPROVED_UPPER_TO_CANONICAL = {str(s).upper(): str(s) for s in approved}
        _ALIASES = {str(k).upper(): str(v) for k, v in aliases.items()}
        return True
    except (OSError, ValueError):
        if not _WARNED_MISSING:
            log.warning(
                "HGNC alias map at %s is unreadable — gene-symbol "
                "canonicalisation disabled (fail-open).",
                path,
                exc_info=True,
            )
            _WARNED_MISSING = True
        _APPROVED_UPPER_TO_CANONICAL = None
        _ALIASES = None
        return False


def canonicalise_gene_symbol(gene: str) -> dict:
    """Resolve ``gene`` to its current HGNC-approved symbol via the local map.

    Returns a dict with keys ``input``, ``approved``, ``is_alias``,
    ``recognized`` and ``ambiguous``:

      - input already an approved symbol (case-insensitive) →
        ``approved`` = HGNC canonical casing, ``is_alias=False``,
        ``recognized=True``.
      - input is a known alias / previous symbol →
        ``approved`` = mapped approved symbol, ``is_alias=True``,
        ``recognized=True``.
      - input is unknown → ``approved=None``, ``is_alias=False``,
        ``recognized=False``.
      - the local map is absent/unreadable → NEUTRAL fail-open:
        ``approved`` = the input unchanged, ``is_alias=False``,
        ``recognized=True`` (never block a curation).

    ``ambiguous`` is always False on the local path (ambiguous aliases are
    dropped at build time); the key is reserved for a future live fallback.
    """
    original = gene
    lookup = (gene or "").strip().upper()

    if not _load_map():
        return {
            "input": original,
            "approved": gene,
            "is_alias": False,
            "recognized": True,
            "ambiguous": False,
        }

    assert _APPROVED_UPPER_TO_CANONICAL is not None
    assert _ALIASES is not None

    if not lookup:
        return {
            "input": original,
            "approved": None,
            "is_alias": False,
            "recognized": False,
            "ambiguous": False,
        }

    canonical = _APPROVED_UPPER_TO_CANONICAL.get(lookup)
    if canonical is not None:
        return {
            "input": original,
            "approved": canonical,
            "is_alias": False,
            "recognized": True,
            "ambiguous": False,
        }

    mapped = _ALIASES.get(lookup)
    if mapped is not None:
        return {
            "input": original,
            "approved": mapped,
            "is_alias": True,
            "recognized": True,
            "ambiguous": False,
        }

    return {
        "input": original,
        "approved": None,
        "is_alias": False,
        "recognized": False,
        "ambiguous": False,
    }
