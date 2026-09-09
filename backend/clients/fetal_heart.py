"""Local-DB-backed fetal cardiac expression lookup.

Reads from ``data/fetal_heart.db`` built by ``build_fetal_heart_db.py``.
Returns per-(cell_type, stage) pseudobulk mean expression (log1p CPM)
and the fraction of cells in each group expressing the gene, derived
from Farah et al. 2024 ("Heart of Cells", Nature 627:854) via the
UCSC Cell Browser.

The qualitative bands (low / moderate / high / very high) are
recalibrated against this dataset's distribution at module-load time
using the non-zero mean_expr quantiles — fetal scRNA-seq mean log1p(CPM)
is on a different scale than GTEx bulk TPM, so reusing GTEx cut-points
would misclassify everything as "low".

Rebuild the database monthly:

    python3 scripts/build_fetal_heart_db.py
"""

from __future__ import annotations

import asyncio
import logging
import sqlite3

from ..localio import run_local
from ._paths import db_path

log = logging.getLogger("heartvar.fetal_heart")

DB_PATH = db_path("fetal_heart.db", "FETAL_HEART_DB_PATH")

DATASET_LABEL = "Farah2024_HOC"
PORTAL_URL = "https://cells.ucsc.edu/?ds=hoc"

CELL_TYPE_LABELS: dict[str, str] = {
    "vCM":         "Ventricular cardiomyocyte",
    "aCM":         "Atrial cardiomyocyte",
    "ncCM":        "Non-chamber cardiomyocyte",
    "Fibro":       "Fibroblast",
    "SMC":         "Smooth muscle cell",
    "BEC":         "Blood endothelial cell",
    "LEC":         "Lymphatic endothelial cell",
    "Endocardial": "Endocardial cell",
    "Epicardial":  "Epicardial cell",
    "WBC":         "Immune (WBC)",
    "Neuronal":    "Neural",
    "P-RBC":       "Erythroid precursor",
}

CELL_TYPE_ORDER: list[str] = [
    "vCM", "aCM", "ncCM",
    "Endocardial", "BEC", "LEC",
    "Fibro", "SMC", "Epicardial",
    "WBC", "Neuronal", "P-RBC",
]

_BANDS: tuple[float, float, float] | None = None
_BAND_LOCK = asyncio.Lock()


def _load_bands_sync() -> tuple[float, float, float] | None:
    """Read the dataset-relative p50/p80/p95 thresholds from the
    ``fetal_heart_meta`` companion table written at build time. Falls
    back to ``None`` if the DB is missing or hasn't been rebuilt with
    the metadata table (older builds)."""
    if not DB_PATH.exists():
        return None
    conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    try:
        cur = conn.cursor()
        cur.execute(
            "SELECT key, value FROM fetal_heart_meta "
            "WHERE key IN ('band_p50','band_p80','band_p95')"
        )
        rows = dict(cur.fetchall())
    except sqlite3.Error:
        log.exception("fetal_heart band-threshold lookup failed")
        return None
    finally:
        conn.close()
    try:
        return (float(rows["band_p50"]), float(rows["band_p80"]), float(rows["band_p95"]))
    except (KeyError, ValueError):
        return None


async def _get_bands() -> tuple[float, float, float] | None:
    global _BANDS
    if _BANDS is not None:
        return _BANDS
    async with _BAND_LOCK:
        if _BANDS is not None:
            return _BANDS
        _BANDS = await run_local(_load_bands_sync)
        return _BANDS


def _classify_band(
    mean_expr: float, bands: tuple[float, float, float] | None
) -> str:
    """Map a mean log1p(CPM) value to one of low / moderate / high /
    very high. Zero / near-zero rows always classify as ``absent`` so
    the UI can distinguish "gene silent here" from "gene present at
    low level"."""
    if mean_expr <= 0:
        return "absent"
    if bands is None:
        return "low"
    p50, p80, p95 = bands
    if mean_expr >= p95:
        return "very_high"
    if mean_expr >= p80:
        return "high"
    if mean_expr >= p50:
        return "moderate"
    return "low"


def _query_sync(gene: str) -> dict:
    """Synchronous DB pull for ``fetch_fetal_heart``."""
    if not DB_PATH.exists():
        return {
            "ok": False,
            "error": (
                f"Fetal heart local DB not found at {DB_PATH}. "
                "Run `python3 scripts/build_fetal_heart_db.py` to build it."
            ),
        }
    conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    try:
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        cur.execute(
            """
            SELECT cell_type, stage, mean_expr, pct_expressing, n_cells
            FROM fetal_heart_expression
            WHERE gene_symbol = ?
            """,
            (gene,),
        )
        rows = cur.fetchall()
    except sqlite3.Error as e:
        log.exception("fetal_heart DB query failed for %s", gene)
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    finally:
        conn.close()
    return {"ok": True, "gene": gene, "rows": [dict(r) for r in rows]}


async def fetch_fetal_heart(gene: str) -> dict:
    """Look up per-(cell_type, stage) fetal cardiac expression for ``gene``.

    Returns the mean log1p(CPM) and the fraction of cells expressing the
    gene in each (cell_type, stage) group from the Farah 2024 dataset.
    Cell-type-keyed entries are ordered by the panel's display order and
    sorted by stage within each cell type. The qualitative band is
    recalibrated against the dataset's own non-zero distribution.
    """
    gene = (gene or "").strip()
    if not gene:
        return {"ok": False, "error": "gene is required"}

    bands = await _get_bands()
    raw = await run_local(_query_sync, gene)
    if not raw.get("ok"):
        return raw

    rows = raw.get("rows") or []
    if not rows:
        return {
            "ok": True,
            "found": False,
            "gene": gene,
            "dataset": DATASET_LABEL,
            "url": PORTAL_URL,
        }

    stages_seen = sorted(
        {r["stage"] for r in rows},
        key=lambda s: int(s.split()[0]) if s.split()[0].isdigit() else 99,
    )

    order_index = {ct: i for i, ct in enumerate(CELL_TYPE_ORDER)}
    rows.sort(
        key=lambda r: (
            order_index.get(r["cell_type"], 999),
            stages_seen.index(r["stage"]) if r["stage"] in stages_seen else 99,
        )
    )

    cell_types_out = []
    for r in rows:
        cell_types_out.append({
            "cell_type": r["cell_type"],
            "cell_type_label": CELL_TYPE_LABELS.get(r["cell_type"], r["cell_type"]),
            "stage": r["stage"],
            "mean_expr": round(float(r["mean_expr"]), 3),
            "pct_expressing": round(float(r["pct_expressing"]), 1),
            "n_cells": int(r["n_cells"]),
            "band": _classify_band(float(r["mean_expr"]), bands),
        })

    return {
        "ok": True,
        "found": True,
        "gene": gene,
        "dataset": DATASET_LABEL,
        "stages": stages_seen,
        "cell_types": cell_types_out,
        "band_thresholds": (
            {"moderate": round(bands[0], 3),
             "high": round(bands[1], 3),
             "very_high": round(bands[2], 3)}
            if bands else None
        ),
        "url": PORTAL_URL,
    }
