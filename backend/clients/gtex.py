from __future__ import annotations

import asyncio
import logging
import os
import sqlite3
from pathlib import Path

import httpx

from ..localio import run_local
from ._http_retry import make_async_client, request_with_retry
from ._offline import offline_strict
from ._paths import PROJECT_ROOT

log = logging.getLogger("heartvar.gtex")

GTEX_API = "https://gtexportal.org/api/v2"
DATASET = "gtex_v10"
GENCODE_VERSION = "v39"
CARDIOVASCULAR_TISSUES = (
    "Heart_Left_Ventricle",
    "Heart_Atrial_Appendage",
    "Artery_Aorta",
    "Artery_Coronary",
)

_PROJECT_ROOT = PROJECT_ROOT
DB_PATH = Path(os.environ.get("GTEX_DB_PATH") or (_PROJECT_ROOT / "data" / "gtex.db"))

_GENCODE_CACHE: dict[str, str | None] = {}
_GENCODE_LOCK = asyncio.Lock()


def _db_path() -> Path:
    """Resolve the DB path at call time so the GTEX_DB_PATH env override (or
    a test monkeypatching the module attribute) is always honoured."""
    return DB_PATH


def _query_local_sync(gene: str) -> dict | None:
    """Synchronous local-DB lookup. Returns the byte-identical fetch_gtex
    success dict on a hit, or ``None`` when the DB file is absent or the
    gene has no GTEx entry — both of which mean the caller should fall back
    to the live API."""
    db_path = _db_path()
    if not db_path.exists():
        return None
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    except sqlite3.Error:
        log.exception("GTEx local DB open failed at %s", db_path)
        return None
    try:
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        cur.execute(
            """
            SELECT gene_symbol, gencode_id, tissue_id, median_tpm
            FROM gtex_cardiac_expression
            WHERE UPPER(gene_symbol) = UPPER(?)
            """,
            (gene,),
        )
        rows = cur.fetchall()
    except sqlite3.Error:
        log.exception("GTEx local DB query failed for %s", gene)
        return None
    finally:
        conn.close()

    if not rows:
        return None

    gencode_id = rows[0]["gencode_id"]
    by_tissue = {r["tissue_id"]: r["median_tpm"] for r in rows}
    return {
        "ok": True,
        "found": True,
        "gene": gene,
        "gencode_id": gencode_id,
        "dataset": DATASET,
        "tissues": [
            {
                "tissue_id": t,
                "tissue_label": t.replace("_", " "),
                "median_tpm": by_tissue.get(t),
            }
            for t in CARDIOVASCULAR_TISSUES
        ],
        "url": f"https://gtexportal.org/home/gene/{gencode_id}",
    }


async def _resolve_gencode_id(c: httpx.AsyncClient, gene: str) -> str | None:
    """Return the GENCODE ID for ``gene`` under GENCODE_VERSION, caching
    the symbol→id mapping for the life of the process. Returns ``None``
    when the gene is not in GTEx."""
    if gene in _GENCODE_CACHE:
        return _GENCODE_CACHE[gene]
    async with _GENCODE_LOCK:
        if gene in _GENCODE_CACHE:
            return _GENCODE_CACHE[gene]
        lookup = await request_with_retry(
            c, "GET", f"{GTEX_API}/reference/gene",
            params={"geneId": gene, "gencodeVersion": GENCODE_VERSION,
                    "page": 0, "itemsPerPage": 5},
            timeout=30.0, name="GTEx/lookup",
        )
        if lookup is None:
            raise httpx.RequestError("GTEx gene lookup transient failure after retries")
        if lookup.status_code != 200:
            raise httpx.HTTPStatusError(
                f"gene lookup {lookup.status_code}: {lookup.text[:200]}",
                request=lookup.request, response=lookup,
            )
        matches = lookup.json().get("data", []) or []
        if not matches:
            _GENCODE_CACHE[gene] = None
            return None
        chosen = next((m for m in matches if m.get("geneSymbol") == gene), matches[0])
        gencode_id = chosen.get("gencodeId")
        _GENCODE_CACHE[gene] = gencode_id
        return gencode_id


async def _fetch_gtex_live(gene: str) -> dict:
    """Live GTEx Portal lookup of the gene's median expression in
    cardiovascular tissues (heart LV / atrial appendage + aorta / coronary
    artery) from GTEx v10. Two-step: gene symbol → GENCODE ID → median
    expression. The first step is cached in-process (see _resolve_gencode_id).
    This is the FALLBACK path — fetch_gtex prefers the local DB and calls here
    only when the DB is absent or the gene is not found locally. Cardiac
    expression supports PP2 for missense variants in established cardiac-disease
    genes; aortic / coronary expression extends the evidence to aortopathy and
    vascular genes.
    """
    async with make_async_client() as c:
        try:
            gencode_id = await _resolve_gencode_id(c, gene)
        except httpx.HTTPStatusError as e:
            return {"ok": False, "error": str(e)}
        if gencode_id is None:
            return {"ok": True, "found": False, "gene": gene,
                    "error": f"GTEx has no entry for {gene} in GENCODE {GENCODE_VERSION}"}

        params = [
            ("gencodeId", gencode_id),
            ("datasetId", DATASET),
        ]
        for t in CARDIOVASCULAR_TISSUES:
            params.append(("tissueSiteDetailId", t))
        expr = await request_with_retry(
            c, "GET", f"{GTEX_API}/expression/medianGeneExpression",
            params=params, timeout=30.0, name="GTEx/expression",
        )
        if expr is None:
            return {"ok": False, "error": "GTEx transient failure after retries",
                    "gencode_id": gencode_id}
        if expr.status_code != 200:
            return {"ok": False, "error": f"expression {expr.status_code}: {expr.text[:200]}",
                    "gencode_id": gencode_id}

        rows = expr.json().get("data", []) or []
        by_tissue = {row["tissueSiteDetailId"]: row["median"] for row in rows}
        return {
            "ok": True,
            "found": True,
            "gene": gene,
            "gencode_id": gencode_id,
            "dataset": DATASET,
            "tissues": [
                {
                    "tissue_id": t,
                    "tissue_label": t.replace("_", " "),
                    "median_tpm": by_tissue.get(t),
                }
                for t in CARDIOVASCULAR_TISSUES
            ],
            "url": f"https://gtexportal.org/home/gene/{gencode_id}",
        }


async def fetch_gtex(gene: str) -> dict:
    """Look up the gene's median expression in heart tissues from GTEx v10.

    PREFERS the local SQLite cache built by ``scripts/build_gtex_db.py`` (no
    network, no rate limit). FALLS BACK to the live GTEx Portal API when the
    DB file is absent or the gene is not present locally — preserving the
    original two-step behaviour. Returns the median TPM per heart tissue;
    cardiac expression supports PP2 for missense variants in established
    cardiac-disease genes.

    The return dict shape is identical on both paths:
    ``{ok, found, gene, gencode_id, dataset, tissues:[{tissue_id,
    tissue_label, median_tpm}], url}`` on a hit.
    """
    local = await run_local(_query_local_sync, gene)
    if local is not None:
        return local
    if offline_strict():
        return {"ok": True, "found": False, "gene": gene,
                "error": "GTEx not in local cache (offline-strict mode; live fallback disabled)"}
    return await _fetch_gtex_live(gene)
