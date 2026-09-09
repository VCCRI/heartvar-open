"""Local-DB-backed BioGRID protein-protein interactions client.

Reads from ``data/biogrid.db`` built by ``build_biogrid_db.py``. Replaces
the previous webservice.thebiogrid.org REST client: no network calls,
no API key required (BioGRID's bulk download is freely available for
academic use). Return shape is the same as the legacy client so the
curate pipeline and the frontend evidence renderer don't need to change.

Rebuild the database quarterly:

    python3 scripts/build_biogrid_db.py
"""

from __future__ import annotations

import logging
import sqlite3
from collections import defaultdict

from ..localio import run_local
from ._paths import db_path
from .chdgene import _load as _load_chdgene

log = logging.getLogger("heartvar.biogrid")

DB_PATH = db_path("biogrid.db", "BIOGRID_DB_PATH")


def _chd_symbols() -> set[str]:
    index, _meta = _load_chdgene()
    return {s.upper() for s in (index or {}).keys()}


def _query_sync(gene: str) -> dict:
    """Synchronous DB lookup. Called from the async wrapper via to_thread."""
    target = gene.strip().upper()
    if not DB_PATH.exists():
        return {
            "ok": False,
            "error": (
                f"BioGRID local DB not found at {DB_PATH}. "
                "Run `python3 scripts/build_biogrid_db.py` to build it."
            ),
        }

    conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT symbol_a, symbol_b, pmid, throughput FROM interactions
            WHERE symbol_a = ?
            UNION ALL
            SELECT symbol_a, symbol_b, pmid, throughput FROM interactions
            WHERE symbol_b = ?
            """,
            (target, target),
        )
        rows = cur.fetchall()
    except sqlite3.Error as e:
        log.exception("BioGRID local DB query failed for %s", gene)
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    finally:
        conn.close()

    partner_pubs: dict[str, set[str]] = defaultdict(set)
    partner_lt: dict[str, set[str]] = defaultdict(set)
    total_rows = 0
    low_throughput_rows = 0
    for symbol_a, symbol_b, pmid, throughput in rows:
        a = (symbol_a or "").upper()
        b = (symbol_b or "").upper()
        if a == target and b and b != target:
            partner = b
        elif b == target and a and a != target:
            partner = a
        else:
            continue
        total_rows += 1
        is_low = (throughput or "").strip() == "Low Throughput"
        if is_low:
            low_throughput_rows += 1
        if pmid:
            partner_pubs[partner].add(pmid)
            if is_low:
                partner_lt[partner].add(pmid)

    chd_set = _chd_symbols()
    ranked = sorted(
        partner_pubs.items(),
        key=lambda kv: (
            0 if kv[0] in chd_set else 1,
            -len(partner_lt.get(kv[0], set())),
            -len(kv[1]),
            kv[0],
        ),
    )
    top = []
    for sym, pmids in ranked[:5]:
        top.append({
            "symbol": sym,
            "publication_count": len(pmids),
            "low_throughput_count": len(partner_lt.get(sym, set())),
            "is_chdgene": sym in chd_set,
        })

    chd_interactors = [t["symbol"] for t in top if t["is_chdgene"]]
    all_chd = sorted(p for p in partner_pubs if p in chd_set)

    partners20 = []
    for sym, pmids in ranked[:20]:
        partners20.append({
            "symbol": sym,
            "evidence_count": len(pmids),
            "low_throughput_evidence_count": len(partner_lt.get(sym, set())),
            "is_chdgene": sym in chd_set,
        })

    return {
        "ok": True,
        "gene": gene,
        "total_interactions": total_rows,
        "low_throughput_interactions": low_throughput_rows,
        "unique_partners": len(partner_pubs),
        "top_partners": top,
        "chd_interactors_top": chd_interactors,
        "chd_interactors_all": all_chd,
        "partners": partners20,
        "url": f"https://thebiogrid.org/search.php?search={gene}&organism=9606",
    }


async def fetch_biogrid(gene: str) -> dict:
    """Look up curated physical PPIs for ``gene`` in the local cache.

    Returns the same dict shape as the previous REST client (plus a few
    new keys consumed by the ``/api/biogrid`` route). SQLite calls are
    wrapped in ``run_local`` so the parallel DB gather isn't
    blocked.
    """
    gene = (gene or "").strip()
    if not gene:
        return {"ok": False, "error": "gene is required"}
    return await run_local(_query_sync, gene)
