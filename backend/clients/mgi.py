"""MGI mouse orthologue + phenotype lookup.

PREFERS a local SQLite cache (``data/mgi.db``, built by
``scripts/build_mgi_db.py`` from MGI's public flat reports + the MP
ontology) and FALLS BACK to the live three-host lookup chain
(genenames.org HGNC id → Alliance of Genome Resources orthologs →
Alliance phenotypes) when the DB is absent OR the gene is not found
locally. The local path eliminates the three per-variant HTTP calls that
risk rate-limiting on public deployment.

The public ``fetch_mgi(gene)`` signature and its return-dict shape are
identical on both paths — ``app.py``, the prompt builder, and the
frontend all consume specific keys (ok, found, gene, hgnc_id,
mouse_symbol, mgi_id, ortholog_confidence, is_best_score,
phenotype_count, cardiac_phenotypes, phenotype_sample, pubmed_count,
url, alliance_url). Note that ``ortholog_confidence`` and
``is_best_score`` are NOT present in the bulk reports, so they are
``None`` on the local path (the keys still exist in the dict).

Rebuild the database monthly::

    python3 scripts/build_mgi_db.py
"""

from __future__ import annotations

import asyncio
import logging
import os
import sqlite3
import urllib.parse
from pathlib import Path
from typing import Any

import httpx

from ..localio import run_local
from ._http_retry import make_async_client, request_with_retry
from ._offline import offline_strict
from ._paths import PROJECT_ROOT

log = logging.getLogger("heartvar.mgi")

HGNC_API = "https://rest.genenames.org/fetch/symbol"
ALLIANCE_API = "https://www.alliancegenome.org/api"

_PROJECT_ROOT = PROJECT_ROOT
DB_PATH = Path(os.environ.get("MGI_DB_PATH") or (_PROJECT_ROOT / "data" / "mgi.db"))

_HGNC_CACHE: dict[str, str | None] = {}
_ORTHO_CACHE: dict[str, dict[str, Any] | None] = {}
_HGNC_LOCK = asyncio.Lock()
_ORTHO_LOCK = asyncio.Lock()

_CARDIAC_KEYWORDS = (
    "heart",
    "cardiac",
    "cardiovascular",
    "cardiomyopath",
    "myocard",
    "ventric",
    "atrial",
    "atrium",
    "aort",
    "valve",
)


async def _hgnc_id(client: httpx.AsyncClient, gene: str) -> str | None:
    if gene in _HGNC_CACHE:
        return _HGNC_CACHE[gene]
    async with _HGNC_LOCK:
        if gene in _HGNC_CACHE:
            return _HGNC_CACHE[gene]
        r = await request_with_retry(
            client, "GET", f"{HGNC_API}/{urllib.parse.quote(gene, safe='')}",
            headers={"Accept": "application/json"},
            timeout=30.0, name="MGI/hgnc",
        )
        if r is None or r.status_code != 200:
            return None
        docs = ((r.json() or {}).get("response") or {}).get("docs") or []
        hgnc_id = docs[0].get("hgnc_id") if docs else None
        _HGNC_CACHE[gene] = hgnc_id
        return hgnc_id


async def _mouse_orthologue(client: httpx.AsyncClient, hgnc_id: str) -> dict[str, Any] | None:
    if hgnc_id in _ORTHO_CACHE:
        return _ORTHO_CACHE[hgnc_id]
    async with _ORTHO_LOCK:
        if hgnc_id in _ORTHO_CACHE:
            return _ORTHO_CACHE[hgnc_id]
        r = await request_with_retry(
            client, "GET", f"{ALLIANCE_API}/gene/{hgnc_id}/orthologs",
            timeout=20.0, name="MGI/orthologs",
        )
        if r is None or r.status_code != 200:
            return None
        payload = r.json() or {}
        for result in payload.get("results") or []:
            gen = result.get("geneToGeneOrthologyGenerated") or {}
            obj = gen.get("objectGene") or {}
            taxon = (obj.get("taxon") or {}).get("curie") or ""
            if taxon == "NCBITaxon:10090":
                ortho = {
                    "mgi_id": obj.get("primaryExternalId"),
                    "symbol": (obj.get("geneSymbol") or {}).get("displayText"),
                    "confidence": (gen.get("confidence") or {}).get("name"),
                    "is_best_score": (gen.get("isBestScore") or {}).get("name") == "Yes",
                }
                _ORTHO_CACHE[hgnc_id] = ortho
                return ortho
        _ORTHO_CACHE[hgnc_id] = None
        return None


async def _phenotypes(client: httpx.AsyncClient, mgi_id: str) -> list[dict[str, Any]]:
    r = await request_with_retry(
        client, "GET", f"{ALLIANCE_API}/gene/{mgi_id}/phenotypes",
        params={"limit": "200"},
        timeout=20.0, name="MGI/phenotypes",
    )
    if r is None or r.status_code != 200:
        return []
    payload = r.json() or {}
    return payload.get("results") or []


def _is_cardiac(statement: str) -> bool:
    s = statement.lower()
    return any(k in s for k in _CARDIAC_KEYWORDS)


def _build_result(
    gene: str,
    hgnc_id: str | None,
    mouse_symbol: str | None,
    mgi_id: str,
    ortholog_confidence: str | None,
    is_best_score: bool | None,
    statements: list[str],
    pubmed_ids: set[str],
) -> dict:
    """Assemble the canonical fetch_mgi return dict. Shared by the live
    and local paths so both produce a byte-identical shape."""
    cardiac = [st for st in statements if _is_cardiac(st)]
    return {
        "ok": True,
        "found": True,
        "gene": gene,
        "hgnc_id": hgnc_id,
        "mouse_symbol": mouse_symbol,
        "mgi_id": mgi_id,
        "ortholog_confidence": ortholog_confidence,
        "is_best_score": is_best_score,
        "phenotype_count": len(statements),
        "cardiac_phenotypes": sorted(set(cardiac)),
        "phenotype_sample": statements[:8],
        "pubmed_count": len(pubmed_ids),
        "url": f"https://www.informatics.jax.org/marker/{mgi_id}",
        "alliance_url": f"https://www.alliancegenome.org/gene/{mgi_id}",
    }


def _query_local_sync(gene: str) -> dict | None:
    """Synchronous local-DB lookup. Returns the canonical fetch_mgi dict
    on a hit, or ``None`` to signal a miss (DB absent, query error, or the
    gene not present in ``gene_ortholog``) so the caller falls back to the
    live chain."""
    if not DB_PATH.exists():
        return None
    try:
        conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    except sqlite3.Error as e:
        log.warning("MGI local DB open failed (%s) — falling back to live", e)
        return None
    try:
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        ortho = cur.execute(
            "SELECT hgnc_id, mgi_id, mouse_symbol "
            "FROM gene_ortholog WHERE human_symbol = ?",
            (gene,),
        ).fetchone()
        if ortho is None or not ortho["mgi_id"]:
            return None
        mgi_id = ortho["mgi_id"]
        rows = cur.execute(
            "SELECT mp_id, phenotype_statement, pubmed_id "
            "FROM gene_phenotype WHERE mgi_id = ?",
            (mgi_id,),
        ).fetchall()
    except sqlite3.Error:
        log.exception("MGI local DB query failed for %s", gene)
        return None
    finally:
        conn.close()

    statements: list[str] = []
    seen_mp: set[str] = set()
    pubmed_ids: set[str] = set()
    for row in rows:
        mp_id = row["mp_id"]
        st = (row["phenotype_statement"] or "").strip()
        if st and mp_id not in seen_mp:
            seen_mp.add(mp_id)
            statements.append(st)
        pmid = (row["pubmed_id"] or "").strip()
        if pmid:
            pubmed_ids.add(pmid)

    return _build_result(
        gene=gene,
        hgnc_id=ortho["hgnc_id"],
        mouse_symbol=ortho["mouse_symbol"],
        mgi_id=mgi_id,
        ortholog_confidence=None,
        is_best_score=None,
        statements=statements,
        pubmed_ids=pubmed_ids,
    )


async def _fetch_mgi_live(gene: str) -> dict:
    """Original live implementation: genenames.org → Alliance orthologs →
    Alliance phenotypes. Used as the fallback when the local DB is absent
    or the gene is not found locally."""
    try:
        async with make_async_client() as c:
            hgnc_id = await _hgnc_id(c, gene)
            if not hgnc_id:
                return {"ok": True, "found": False, "error_step": "hgnc", "gene": gene}
            ortho = await _mouse_orthologue(c, hgnc_id)
            if not ortho or not ortho.get("mgi_id"):
                return {
                    "ok": True,
                    "found": False,
                    "gene": gene,
                    "hgnc_id": hgnc_id,
                    "reason": "no mouse orthologue indexed",
                }
            phenos = await _phenotypes(c, ortho["mgi_id"])
    except httpx.HTTPError as e:
        return {"ok": False, "error": f"transport error: {e!r}"}

    statements: list[str] = []
    pubmed_ids: set[str] = set()
    for p in phenos:
        st = (p.get("phenotypeStatement") or "").strip()
        if not st:
            continue
        statements.append(st)
        for ref in p.get("pubmedPubModIDs") or []:
            ref_str = str(ref).strip()
            if ref_str.startswith("PMID:"):
                pubmed_ids.add(ref_str.split(":", 1)[1])
            elif ref_str.isdigit():
                pubmed_ids.add(ref_str)

    return _build_result(
        gene=gene,
        hgnc_id=hgnc_id,
        mouse_symbol=ortho.get("symbol"),
        mgi_id=ortho["mgi_id"],
        ortholog_confidence=ortho.get("confidence"),
        is_best_score=ortho.get("is_best_score"),
        statements=statements,
        pubmed_ids=pubmed_ids,
    )


async def fetch_mgi(gene: str) -> dict:
    """Look up a human gene's mouse orthologue + phenotype context.

    Prefers the local ``data/mgi.db`` cache; falls back to the live
    three-host lookup chain when the DB is absent or the gene is not
    found locally. Return shape is identical on both paths.
    """
    gene = (gene or "").strip()
    if not gene:
        return {"ok": False, "error": "gene is required"}
    local = await run_local(_query_local_sync, gene)
    if local is not None:
        return local
    if offline_strict():
        return {"ok": True, "found": False, "error_step": "hgnc", "gene": gene}
    return await _fetch_mgi_live(gene)
