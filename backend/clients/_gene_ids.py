"""Local gene-ID cross-reference — Ensembl gene id from a symbol or Entrez id.

VEP is the only place a curation learns a gene's Ensembl id, and it does NOT
always supply one. RefSeq-prefixed input (NM_/NR_/XM_/XR_) is queried against
VEP's RefSeq transcript set (``refseq=1`` — see ensembl_vep._REFSEQ_PARAMS),
whose per-transcript ``gene_id`` is an **EntrezGene** id, not an Ensembl one:

    NM_000501.3:c.1537G>A (ELN)  ->  gene_id = "2006"      (Entrez)
    ENST00000355349      (MYH7)  ->  gene_id = "ENSG00000092054"

Any consumer keyed on an ``ENSG`` therefore has to resolve one itself, or it
silently reports "unavailable" for every RefSeq-transcript curation — which is
essentially all clinical input.

Resolution is OFFLINE, from one of two local sources tried in this order, so it
works under ``HEARTVAR_OFFLINE_STRICT`` with no network call:

1. ``backend/data/gene_id_map.json.gz`` — the compact map built by
   ``scripts/build_gene_id_map.py``. TRACKED and gzipped (~0.5 MB), so it is
   COPYd into the runtime image with the rest of ``backend/data/``.
2. ``backend/data/hgnc_complete_set.txt`` — the full 17 MB HGNC set, when
   present. Dev machines and the data-builder image have it.

Source 1 exists because source 2 alone is a production trap: the HGNC set is
both gitignored AND listed in ``.dockerignore``, so it is a build-time cache
that never reaches the deployed container. A runtime-only dependency on it works
perfectly on a dev machine and silently fails-open in production, leaving Open
Targets "unavailable" for every RefSeq curation while looking correct locally.

Fail-open: with neither source present, resolution is disabled and warns once
rather than blocking a curation.
"""

from __future__ import annotations

import gzip
import json
import logging
import os
from pathlib import Path

log = logging.getLogger("heartvar.gene_ids")

_BACKEND_ROOT = Path(__file__).resolve().parent.parent
GENE_ID_MAP_PATH = Path(
    os.environ.get("GENE_ID_MAP_PATH")
    or (_BACKEND_ROOT / "data" / "gene_id_map.json.gz")
)
HGNC_TSV_PATH = Path(
    os.environ.get("HGNC_COMPLETE_SET_PATH")
    or (_BACKEND_ROOT / "data" / "hgnc_complete_set.txt")
)

_SYMBOL_TO_ENSG: dict[str, str] | None = None
_ENTREZ_TO_ENSG: dict[str, str] | None = None
_LOADED = False
_WARNED_MISSING = False


def _load_compact_map() -> tuple[dict[str, str], dict[str, str]] | None:
    """Read the shipped ``gene_id_map.json.gz`` (see scripts/build_gene_id_map.py).

    Returns (symbol→ENSG, entrez→ENSG) or None when the file is absent or
    unusable, so the caller can fall back to the full HGNC set. Values are
    stored as the ENSG integer suffix and rehydrated here.
    """
    path = GENE_ID_MAP_PATH
    if not path.exists():
        return None
    try:
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            data = json.load(fh)
        raw_sym = data["symbol_to_ensg"]
        raw_ent = data["entrez_to_ensg"]
    except (OSError, ValueError, KeyError, TypeError):
        log.warning(
            "gene id map at %s is unreadable — falling back to the HGNC set.",
            path, exc_info=True,
        )
        return None
    try:
        sym = {k: f"ENSG{int(v):011d}" for k, v in raw_sym.items()}
        ent = {k: f"ENSG{int(v):011d}" for k, v in raw_ent.items()}
    except (TypeError, ValueError):
        log.warning("gene id map at %s has non-integer values — ignoring.", path)
        return None
    if not sym and not ent:
        return None
    return sym, ent


def _load() -> bool:
    """Lazily populate the process caches. Returns True on success, False when
    no source is available (fail-open). Loads at most once; subsequent calls are
    no-ops once ``_LOADED`` is set.

    Prefers the small shipped map, because the 17 MB HGNC set it was built from
    is excluded from the runtime image (see the module docstring). Falls back to
    that set when present — dev machines and the builder image have it, and it is
    authoritative if someone has refreshed it without rebuilding the map.
    """
    global _SYMBOL_TO_ENSG, _ENTREZ_TO_ENSG, _LOADED, _WARNED_MISSING
    if _LOADED:
        return _SYMBOL_TO_ENSG is not None
    _LOADED = True

    compact = _load_compact_map()
    if compact is not None:
        _SYMBOL_TO_ENSG, _ENTREZ_TO_ENSG = compact
        return True

    path = HGNC_TSV_PATH
    if not path.exists():
        if not _WARNED_MISSING:
            log.warning(
                "No gene-ID source: neither %s nor %s exists — Ensembl gene-id "
                "resolution disabled (fail-open), so Open Targets will be "
                "unavailable for RefSeq-transcript curations. Rebuild the map "
                "with scripts/build_gene_id_map.py. See scripts/README.md.",
                GENE_ID_MAP_PATH, path,
            )
            _WARNED_MISSING = True
        return False
    sym: dict[str, str] = {}
    ent: dict[str, str] = {}
    try:
        with path.open(encoding="utf-8") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                i_sym = header.index("symbol")
                i_ensg = header.index("ensembl_gene_id")
                i_entrez = header.index("entrez_id")
            except ValueError:
                if not _WARNED_MISSING:
                    log.warning(
                        "hgnc_complete_set.txt at %s is missing symbol / "
                        "ensembl_gene_id / entrez_id columns — Ensembl gene-id "
                        "resolution disabled (fail-open).",
                        path,
                    )
                    _WARNED_MISSING = True
                return False
            width = max(i_sym, i_ensg, i_entrez) + 1
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < width:
                    continue
                ensg = f[i_ensg].strip()
                if not ensg.startswith("ENSG"):
                    continue
                symbol = f[i_sym].strip().upper()
                if symbol:
                    sym.setdefault(symbol, ensg)
                entrez = f[i_entrez].strip()
                if entrez:
                    ent.setdefault(entrez, ensg)
    except OSError:
        if not _WARNED_MISSING:
            log.warning(
                "hgnc_complete_set.txt at %s is unreadable — Ensembl gene-id "
                "resolution disabled (fail-open).",
                path,
                exc_info=True,
            )
            _WARNED_MISSING = True
        return False
    _SYMBOL_TO_ENSG = sym
    _ENTREZ_TO_ENSG = ent
    return True


def resolve_ensembl_gene_id(
    gene_id: str | None = None,
    gene_symbol: str | None = None,
) -> str | None:
    """Best-effort Ensembl gene id (``ENSG…``) for a VEP ``gene_id`` / symbol.

    Tried in order, most to least authoritative:

      1. ``gene_id`` is already an ``ENSG`` — returned as-is (version suffix
         stripped, so ``ENSG00000049540.16`` → ``ENSG00000049540``).
      2. ``gene_id`` is all digits — treated as an EntrezGene id (what VEP's
         RefSeq set returns) and mapped via the HGNC ``entrez_id`` column.
      3. ``gene_symbol`` — mapped via the HGNC ``symbol`` column, after
         canonicalising aliases/previous symbols where that map is available.

    Returns None when nothing resolves (unknown gene, or the HGNC set is not
    provisioned). Never raises.
    """
    raw = (gene_id or "").strip()
    if raw.startswith("ENSG"):
        return raw.split(".", 1)[0]

    symbol = (gene_symbol or "").strip().upper()
    if not raw.isdigit() and not symbol:
        return None

    if not _load():
        return None
    assert _SYMBOL_TO_ENSG is not None
    assert _ENTREZ_TO_ENSG is not None

    if raw.isdigit():
        hit = _ENTREZ_TO_ENSG.get(raw)
        if hit:
            return hit

    if not symbol:
        return None
    hit = _SYMBOL_TO_ENSG.get(symbol)
    if hit:
        return hit

    try:
        from .hgnc_alias import canonicalise_gene_symbol

        approved = (canonicalise_gene_symbol(symbol) or {}).get("approved")
    except Exception:  # noqa: BLE001 — resolution is best-effort, never fatal
        return None
    if approved:
        approved_upper = str(approved).strip().upper()
        if approved_upper != symbol:
            return _SYMBOL_TO_ENSG.get(approved_upper)
    return None
