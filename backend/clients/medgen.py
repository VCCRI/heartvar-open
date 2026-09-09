"""MedGen (NCBI) — gene-associated conditions.

PREFERS a local SQLite cache (``data/medgen.db``, built from NCBI
ClinVar's public ``gene_condition_source_id`` table by
``scripts/build_medgen_db.py``) and FALLS BACK to the live NCBI
E-utilities path (esearch + elink + esummary) when the DB file is
absent or the gene isn't present locally. The local read costs zero of
the shared 3 req/s NCBI budget that the literature clients depend on.

The live fallback resolves a gene symbol to its NCBI Gene ID with
``esearch``, walks the ``gene → medgen`` ``elink`` to get associated
MedGen UIDs, then batch-fetches their summaries with ``esummary``.
Records are filtered to keep only ``Disease or Syndrome`` / ``Finding``
SemanticTypes, then capped at ``MAX_CONDITIONS``.

OMIM URLs are constructed locally from MIM numbers (surfaced by the
local table, or by MedGen's ``conceptmeta`` XML fragment on the live
path). No OMIM API is contacted and no OMIM data is cached or logged —
this client is deliberately MedGen-only so the project stays viable as
an open-source / no-paid-API derivative.

Rebuild the local cache monthly::

    python3 scripts/build_medgen_db.py
"""

from __future__ import annotations

import logging
import os
import re
import sqlite3
from pathlib import Path
from typing import Any

import httpx

from ..localio import run_local
from ._http_retry import make_async_client
from ._ncbi_throttle import NCBI_API_KEY, NCBIError, eutils_get_json
from ._offline import offline_strict
from ._paths import PROJECT_ROOT

log = logging.getLogger("heartvar.medgen")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

_PROJECT_ROOT = PROJECT_ROOT


def _db_path() -> Path:
    """Resolve the local MedGen DB path at call time so a MEDGEN_DB_PATH
    set after import (e.g. by tests) is honoured."""
    override = os.environ.get("MEDGEN_DB_PATH", "").strip()
    if override:
        return Path(override)
    return _PROJECT_ROOT / "data" / "medgen.db"


DB_PATH = _PROJECT_ROOT / "data" / "medgen.db"

NCBI_CONTACT_EMAIL = (
    os.environ.get("NCBI_EMAIL", "").strip() or "heartvar@victorchang.edu.au"
)
NCBI_TOOL = "heartvar"

MAX_CONDITIONS = 5

_KEEP_SEMANTIC_TYPES = {
    "Disease or Syndrome",
    "Finding",
}

_HEART_KEYWORDS = (
    "cardio", "cardiac", "heart",
    "cardiomyopath", "myocardi",
    "arrhythmia", "tachycardia", "fibrillation",
    "atrial septal", "ventricular septal",
    "septal defect", "septation",
    "noncompaction", "non-compaction",
    "long qt", "qt syndrome", "brugada",
    "hypoplastic left heart", "hlhs",
    "tetralogy", "aortic stenosis", "aortic regurgitation",
    "patent ductus", "endocardial cushion",
    "channelopath",
    "noonan",
)


def _heart_related(title: str | None) -> bool:
    """Heuristic: does this condition involve the heart?

    Keyword fragments match case-insensitively as substrings against
    the MedGen ``title``. Used to bubble cardiac conditions to the top
    of the conditions list — HeartVar is a cardiac-curation tool so the
    user's attention should land on heart-relevant entries first when
    a gene has both cardiac and non-cardiac associations.
    """
    if not isinstance(title, str) or not title:
        return False
    lower = title.lower()
    return any(kw in lower for kw in _HEART_KEYWORDS)

# XML blob. Attributes within the tag are then searched separately so
# we don't rely on attribute order.
_NAME_TAG_RE = re.compile(r"<Name\b[^>]*>")
_SDUI_RE = re.compile(r'SDUI="(\d+)"')


def _params(**extra) -> dict:
    """Build E-utilities query params with NCBI-policy fields attached."""
    params = dict(extra)
    params["tool"] = NCBI_TOOL
    params["email"] = NCBI_CONTACT_EMAIL
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    return params


def _gene_search_url(gene: str) -> str:
    """Fallback MedGen browser URL — surfaced as the section-level
    "View on MedGen" link and as the error-state fallback when the
    structured lookup fails or returns nothing."""
    return f"https://www.ncbi.nlm.nih.gov/medgen/?term={gene}[gene]"


async def _resolve_gene_id(client: httpx.AsyncClient, gene: str) -> str | None:
    """Resolve a HGNC gene symbol to an NCBI Gene ID via esearch.

    Returns the first hit's UID, or ``None`` when esearch succeeded but
    returned no matches. Propagates :class:`NCBIError` from the shared
    helper so the caller can distinguish "transient failure" from
    "gene legitimately not in NCBI Gene".
    """
    params = _params(
        db="gene",
        term=f"{gene}[sym] AND Homo sapiens[orgn]",
        retmode="json",
        retmax="1",
    )
    payload = await eutils_get_json(
        client, f"{EUTILS}/esearch.fcgi", params, timeout=20.0
    )
    idlist = (payload.get("esearchresult") or {}).get("idlist") or []
    return idlist[0] if idlist else None


async def _gene_to_medgen_uids(
    client: httpx.AsyncClient, ncbi_gene_id: str
) -> list[str]:
    """Walk the gene → medgen elink and return the associated MedGen UIDs."""
    params = _params(
        dbfrom="gene",
        db="medgen",
        id=ncbi_gene_id,
        retmode="json",
    )
    payload = await eutils_get_json(
        client, f"{EUTILS}/elink.fcgi", params, timeout=20.0
    )
    uids: list[str] = []
    for linkset in payload.get("linksets") or []:
        for db in linkset.get("linksetdbs") or []:
            if db.get("dbto") != "medgen":
                continue
            for uid in db.get("links") or []:
                uids.append(str(uid))
    return list(dict.fromkeys(uids))


async def _fetch_summaries(
    client: httpx.AsyncClient, uids: list[str]
) -> list[dict]:
    """ESummary the MedGen UIDs in a single batched call. Empty ``uids``
    short-circuits to ``[]`` without contacting NCBI."""
    if not uids:
        return []
    params = _params(
        db="medgen",
        id=",".join(uids),
        retmode="json",
    )
    payload = await eutils_get_json(
        client, f"{EUTILS}/esummary.fcgi", params, timeout=25.0
    )
    result = payload.get("result") or {}
    out: list[dict] = []
    for uid in result.get("uids") or []:
        rec = result.get(uid)
        if isinstance(rec, dict):
            out.append(rec)
    return out


def _semantic_type_value(record: dict) -> str | None:
    """MedGen esummary returns semantictype as a {"value": "..."} dict —
    flatten to the underlying string."""
    st = record.get("semantictype")
    if isinstance(st, dict):
        return st.get("value")
    if isinstance(st, str):
        return st
    return None


def _record_passes_filter(record: dict) -> bool:
    """Keep only Disease/Syndrome and Finding records. Drop everything
    else (and drop records missing a SemanticType entirely)."""
    st = _semantic_type_value(record)
    return bool(st) and st in _KEEP_SEMANTIC_TYPES


def _extract_mim(record: dict) -> str | None:
    """Pull the canonical OMIM MIM number out of the conceptmeta XML.

    MedGen embeds OMIM cross-references inside ``conceptmeta`` as
    ``<Name SAB="OMIM" SDUI="123456" TTY="…" …>`` elements. Prefer the
    ``TTY="PT"`` (preferred-term) entry — the canonical disease MIM —
    and fall back to the first OMIM entry of any TTY when no PT is
    present. Returns ``None`` when the record has no OMIM cross-reference.
    """
    cm = record.get("conceptmeta")
    if not isinstance(cm, str) or not cm:
        return None
    omim_tags = [m.group(0) for m in _NAME_TAG_RE.finditer(cm) if 'SAB="OMIM"' in m.group(0)]
    for tag in omim_tags:
        if 'TTY="PT"' in tag:
            sdui = _SDUI_RE.search(tag)
            if sdui:
                return sdui.group(1)
    for tag in omim_tags:
        sdui = _SDUI_RE.search(tag)
        if sdui:
            return sdui.group(1)
    return None


def _build_result(
    gene: str,
    fallback_url: str,
    raw_conditions: list[dict],
) -> dict:
    """Turn a list of ``{"name", "cui", "mim"}`` dicts (in source order)
    into the canonical MedGen success payload.

    Shared by the local-DB path and the live E-utilities path so both
    return a byte-identical dict shape. ``raw_conditions`` must already
    be filtered/deduped — ``total_found`` is ``len(raw_conditions)``
    BEFORE the ``MAX_CONDITIONS`` cap, matching the legacy semantics.
    Heart-related conditions are bubbled to the top with a stable sort
    that preserves the source order within each group, then the list is
    capped — so cardiac conditions never get pushed off the visible
    list by a source-earlier non-cardiac entry.
    """
    all_conditions: list[dict] = []
    for idx, rec in enumerate(raw_conditions):
        cui = (rec.get("cui") or "").strip() or None
        name = rec.get("name") or "(unnamed)"
        mim = rec.get("mim")
        all_conditions.append({
            "name": name,
            "cui": cui,
            "mim": mim,
            "moi": None,
            "heart_related": _heart_related(name),
            "omim_url": f"https://www.omim.org/entry/{mim}" if mim else None,
            "medgen_url": f"https://www.ncbi.nlm.nih.gov/medgen/{cui}" if cui else None,
            "_idx": idx,
        })
    all_conditions.sort(key=lambda c: (not c["heart_related"], c["_idx"]))
    conditions = []
    for c in all_conditions[:MAX_CONDITIONS]:
        c.pop("_idx", None)
        conditions.append(c)

    return {
        "ok": True,
        "gene": gene,
        "conditions": conditions,
        "total_found": len(raw_conditions),
        "url": fallback_url,
    }


def _query_local_sync(gene: str) -> list[dict] | None:
    """Return the gene's conditions from the local DB as a list of
    ``{"name", "cui", "mim"}`` dicts in stored (source) order, or
    ``None`` to signal a MISS that must fall back to the live path.

    A miss is: the DB file is absent, the query errors, or the gene has
    no rows. An empty-but-present gene cannot occur (rows only exist for
    genes that have ≥1 association), so ``None`` unambiguously means
    "not resolvable locally — try the network".
    """
    db_path = _db_path()
    if not db_path.exists():
        return None
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    except sqlite3.Error as e:
        log.warning("MedGen local DB open failed (%s); falling back to live", e)
        return None
    try:
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        cur.execute(
            """
            SELECT name, cui, mim
            FROM conditions
            WHERE UPPER(gene_symbol) = UPPER(?)
            ORDER BY rowid
            """,
            (gene,),
        )
        rows = cur.fetchall()
    except sqlite3.Error as e:
        log.warning(
            "MedGen local DB query failed for %s (%s); falling back to live",
            gene, e,
        )
        return None
    finally:
        conn.close()

    if not rows:
        return None
    return [
        {"name": r["name"], "cui": r["cui"], "mim": r["mim"]}
        for r in rows
    ]


async def fetch_medgen(gene: str, hpo_terms: Any = None) -> dict:
    """Look up gene-associated conditions, preferring the local MedGen DB.

    Reads ``data/medgen.db`` (built by ``scripts/build_medgen_db.py``)
    first; on a miss — DB absent, query error, or gene not present
    locally — transparently falls back to the live NCBI E-utilities
    path. The return shape is identical on both paths.

    Parameters
    ----------
    gene : str
        HGNC gene symbol.
    hpo_terms : optional
        Accepted for API symmetry with other gene-keyed clients but
        currently unused — neither the local table nor MedGen
        ``esummary`` surfaces HPO cross-references in a structured form,
        so the optional HPO-ranking step in the spec is a no-op here.

    Returns
    -------
    dict
        On success::

            {
              "ok": True,
              "gene": "MYH7",
              "conditions": [
                {"name": ..., "cui": ..., "mim": "192600",
                 "moi": None,
                 "heart_related": True,
                 "omim_url": "https://www.omim.org/entry/192600",
                 "medgen_url": "https://www.ncbi.nlm.nih.gov/medgen/C…"},
                …
              ],
              "total_found": 12,        # before capping at MAX_CONDITIONS
              "url": "https://www.ncbi.nlm.nih.gov/medgen/?term=MYH7[gene]"
            }

        On any error or empty result::

            {"ok": False, "error": "…", "url": <fallback gene-search URL>}
    """
    gene = (gene or "").strip()
    if not gene:
        return {
            "ok": False,
            "error": "missing gene symbol",
            "url": _gene_search_url(""),
        }
    fallback_url = _gene_search_url(gene)

    raw = await run_local(_query_local_sync, gene)
    if raw is not None:
        return _build_result(gene, fallback_url, raw)

    if offline_strict():
        return {
            "ok": False,
            "error": "MedGen local DB miss (offline-strict mode; live NCBI fallback disabled)",
            "url": fallback_url,
        }

    return await _fetch_medgen_live(gene, hpo_terms)


async def _fetch_medgen_live(gene: str, hpo_terms: Any = None) -> dict:
    """Live NCBI MedGen E-utilities lookup — the graceful fallback used
    when the local DB is unavailable or doesn't carry the gene.

    Preserves the original three-call behaviour (esearch → elink →
    esummary) and the exact legacy return shape. ``gene`` is assumed
    already stripped and non-empty (the public ``fetch_medgen`` guards
    that before delegating here).
    """
    fallback_url = _gene_search_url(gene)

    async with make_async_client() as client:
        try:
            ncbi_gene_id = await _resolve_gene_id(client, gene)
        except NCBIError as e:
            return {
                "ok": False,
                "error": f"NCBI E-utilities request failed ({e})",
                "url": fallback_url,
            }
        if not ncbi_gene_id:
            return {
                "ok": False,
                "error": "gene symbol not found in NCBI Gene",
                "url": fallback_url,
            }
        try:
            uids = await _gene_to_medgen_uids(client, ncbi_gene_id)
        except NCBIError as e:
            return {
                "ok": False,
                "error": f"NCBI E-utilities request failed ({e})",
                "url": fallback_url,
            }
        if not uids:
            return {
                "ok": True,
                "gene": gene,
                "conditions": [],
                "total_found": 0,
                "url": fallback_url,
            }
        try:
            records = await _fetch_summaries(client, uids)
        except NCBIError as e:
            return {
                "ok": False,
                "error": f"NCBI E-utilities request failed ({e})",
                "url": fallback_url,
            }

    filtered = [r for r in records if _record_passes_filter(r)]

    raw_conditions = [
        {
            "name": rec.get("title") or "(unnamed)",
            "cui": (rec.get("conceptid") or "").strip() or None,
            "mim": _extract_mim(rec),
        }
        for rec in filtered
    ]
    return _build_result(gene, fallback_url, raw_conditions)
