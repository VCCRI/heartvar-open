#!/usr/bin/env python3
"""build_fetal_heart_db.py — build the local fetal heart scRNA-seq cache.

Downloads the public Farah 2024 "Heart of Cells" single-cell RNA-seq matrix
and metadata from the UCSC Cell Browser, computes a per-(gene, cell_type,
stage) pseudobulk summary, and writes it to ``data/fetal_heart.db``. The
backend fetal-heart client reads from this database to surface fetal
cardiac expression alongside the adult GTEx panel in HeartVar's Gene
Context tab.

UCSC serves the matrix as RAW integer counts at ``exprMatrix.tsv.gz`` (stored
Uint32 per the dataset's ``matrixArrType``; there is no separate counts file)
plus a wide TSV metadata file. This script normalises counts per cell to
log1p(CPM) before averaging within each (cell_type, stage) group, and
also records the fraction of cells in the group with non-zero expression
(canonical dotplot semantics).

Cell-type grouping uses ``major_cell_class`` (12 broad classes: vCM,
aCM, ncCM, Fibro, SMC, BEC, LEC, Endocardial, Epicardial, WBC, Neuronal,
P-RBC). The public Cell Browser metadata does NOT include the fine
named subtypes (vCM-Early / vCM-Late / etc.) that appear in the paper's
figures — those labels are not in the released annotations.

Build-time dependency: ``numpy`` (declared in ``requirements.txt`` so the
deploy host can run this builder; the runtime itself only reads the resulting
SQLite and does not import numpy).

Usage::

    python3 build_fetal_heart_db.py                 # use cached downloads if present
    python3 build_fetal_heart_db.py --force-download

The cached gzip / TSV live under ``data/``; delete them to force a fresh
download on the next run.
"""

from __future__ import annotations

import csv
import gzip
import shutil
import ssl
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import _dbbuild

import numpy as np

EXPR_URL = "https://cells.ucsc.edu/hoc/all-heart/exprMatrix.tsv.gz"
META_URL = "https://cells.ucsc.edu/hoc/all-heart/meta.tsv"

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DB_PATH = DATA_DIR / "fetal_heart.db"
EXPR_PATH = DATA_DIR / "fetal_heart_exprMatrix.tsv.gz"
META_PATH = DATA_DIR / "fetal_heart_meta.tsv"

DATASET_LABEL = "Farah2024_HOC"

CPM_SCALE = 1_000_000.0

CELL_TYPE_LABELS = {
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

DDL = """
CREATE TABLE fetal_heart_expression (
    gene_symbol     TEXT NOT NULL,
    cell_type       TEXT NOT NULL,
    stage           TEXT NOT NULL,
    mean_expr       REAL NOT NULL,
    pct_expressing  REAL NOT NULL,
    n_cells         INTEGER NOT NULL
)
"""

INSERT_SQL = (
    "INSERT INTO fetal_heart_expression VALUES (?,?,?,?,?,?)"
)

META_DDL = """
CREATE TABLE fetal_heart_meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
)
"""


def _download_urllib(url: str, dest: Path) -> None:
    """Stream the URL into ``dest`` using Python's stdlib."""
    t0 = time.perf_counter()
    bytes_seen = 0
    with urllib.request.urlopen(url) as r, open(dest, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            bytes_seen += len(chunk)
            if bytes_seen % (1 << 25) < (1 << 20):
                print(
                    f"  …{bytes_seen / (1024 * 1024):.1f} MB",
                    flush=True,
                )
    size_mb = dest.stat().st_size / (1024 * 1024)
    print(f"  ✓ {size_mb:.1f} MB in {time.perf_counter() - t0:.1f} s")


def _download_curl(url: str, dest: Path) -> None:
    """Fallback path via system curl — trusts macOS keychain."""
    curl = shutil.which("curl")
    if not curl:
        raise RuntimeError(
            "Python's TLS chain rejected the server certificate (likely an "
            "enterprise inspection proxy) and `curl` isn't on PATH."
        )
    print(f"  …falling back to {curl} (uses system keychain)", flush=True)
    t0 = time.perf_counter()
    cmd = [
        curl, "--location", "--fail", "--retry", "3",
        "--silent", "--show-error",
        "--output", str(dest), url,
    ]
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise RuntimeError(f"curl exited with status {proc.returncode}")
    size_mb = dest.stat().st_size / (1024 * 1024)
    print(f"  ✓ {size_mb:.1f} MB in {time.perf_counter() - t0:.1f} s")


def download(url: str, dest: Path) -> None:
    print(f"Downloading {url}\n  → {dest}", flush=True)
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        _download_urllib(url, dest)
    except (ssl.SSLError, urllib.error.URLError) as e:
        msg = str(e)
        if "CERTIFICATE_VERIFY_FAILED" in msg or "SSL" in msg:
            if dest.exists():
                dest.unlink()
            _download_curl(url, dest)
        else:
            raise


def _format_stage(age_value: str) -> str:
    """UCSC stores PCW as "09" / "11" / "13" / "15". Render as "9 PCW"
    (no leading zero) to match the plan's expected stage labels."""
    s = (age_value or "").strip().lstrip("0") or "0"
    return f"{s} PCW"


def load_meta(meta_path: Path) -> dict:
    """Read meta.tsv into a structured form keyed by cell_id.

    Returns a dict with:
      - ``cell_to_meta``: cell_id → (cell_type, stage, n_count)
      - ``group_keys``  : ordered list of (cell_type, stage) tuples
      - ``group_index`` : (cell_type, stage) → int (column in the per-row aggregates)
    """
    t0 = time.perf_counter()
    cell_to_meta: dict[str, tuple[str, str, float]] = {}
    group_keys: list[tuple[str, str]] = []
    group_index: dict[tuple[str, str], int] = {}

    with open(meta_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if header and header[0] == "":
            header[0] = "cell_id"
        required = {"cell_id", "major_cell_class", "Age", "nCount_RNA"}
        missing = required - set(header)
        if missing:
            raise RuntimeError(
                f"meta.tsv missing expected columns: {sorted(missing)}. "
                f"Header was: {header[:20]}"
            )
        idx_cell = header.index("cell_id")
        idx_class = header.index("major_cell_class")
        idx_age = header.index("Age")
        idx_count = header.index("nCount_RNA")

        for row in reader:
            if len(row) <= max(idx_cell, idx_class, idx_age, idx_count):
                continue
            cell_id = row[idx_cell]
            cell_type = row[idx_class].strip()
            stage = _format_stage(row[idx_age])
            try:
                n_count = float(row[idx_count])
            except ValueError:
                continue
            if not cell_id or not cell_type or n_count <= 0:
                continue
            key = (cell_type, stage)
            if key not in group_index:
                group_index[key] = len(group_keys)
                group_keys.append(key)
            cell_to_meta[cell_id] = (cell_type, stage, n_count)

    print(
        f"  ✓ meta: {len(cell_to_meta):,} cells across {len(group_keys)} "
        f"(cell_type × stage) groups in {time.perf_counter() - t0:.1f} s"
    )
    return {
        "cell_to_meta": cell_to_meta,
        "group_keys": group_keys,
        "group_index": group_index,
    }


def align_meta_to_matrix(
    matrix_header: list[str], meta: dict
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Given the expression-matrix column order, build:

    - ``cell_groups``  : int array, group index for each matrix column
                        (-1 marks columns absent from meta — skipped at agg time)
    - ``inv_total``    : float array, 1/nCount for each matrix column
    - ``n_per_group``  : int array, cells per group (the divisor for means)
    """
    cell_to_meta = meta["cell_to_meta"]
    group_index = meta["group_index"]
    n_groups = len(meta["group_keys"])
    n_cells = len(matrix_header)
    cell_groups = np.full(n_cells, -1, dtype=np.int32)
    inv_total = np.zeros(n_cells, dtype=np.float64)
    n_per_group = np.zeros(n_groups, dtype=np.int64)

    missing = 0
    for j, cell_id in enumerate(matrix_header):
        m = cell_to_meta.get(cell_id)
        if m is None:
            missing += 1
            continue
        cell_type, stage, n_count = m
        gi = group_index[(cell_type, stage)]
        cell_groups[j] = gi
        inv_total[j] = 1.0 / n_count
        n_per_group[gi] += 1

    if missing:
        print(f"  ⚠ {missing} matrix columns absent from meta — skipped at agg time")
    return cell_groups, inv_total, n_per_group


def build_db(expr_path: Path, meta_path: Path, db_path: Path) -> int:
    """Stream the expression matrix, compute pseudobulk, write to SQLite."""
    print(f"Loading metadata from {meta_path}", flush=True)
    meta = load_meta(meta_path)
    group_keys = meta["group_keys"]
    n_groups = len(group_keys)

    print(f"Opening matrix at {expr_path}", flush=True)
    t_total = time.perf_counter()

    tmp_path = _dbbuild.staging_db_path(db_path)
    conn = _dbbuild.connect(tmp_path)
    conn.execute("PRAGMA journal_mode = OFF")
    conn.execute("PRAGMA synchronous = OFF")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute(DDL)
    cur = conn.cursor()

    BATCH = 50_000
    batch: list[tuple] = []
    inserted = 0
    genes_seen = 0
    genes_skipped_empty = 0
    duplicates_skipped = 0
    seen_genes: set[str] = set()

    with gzip.open(expr_path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if header[0].lower() not in ("gene", "geneid", ""):
            raise RuntimeError(
                f"Unexpected exprMatrix header start: {header[:3]!r}; "
                "expected first column to be 'gene'."
            )
        cell_ids = header[1:]
        print(f"  matrix header: 1 gene column + {len(cell_ids):,} cells")

        cell_groups, inv_total, n_per_group = align_meta_to_matrix(cell_ids, meta)
        unused_slot = n_groups
        cell_groups_safe = np.where(cell_groups < 0, unused_slot, cell_groups).astype(np.int64)

        t_stream = time.perf_counter()
        for row in reader:
            if len(row) != 1 + len(cell_ids):
                print(
                    f"  ⚠ malformed row at gene {row[0]!r}: "
                    f"{len(row)} cols vs expected {1 + len(cell_ids)}; skipped",
                    flush=True,
                )
                continue
            gene = row[0].split("|", 1)[0]
            if gene in seen_genes:
                duplicates_skipped += 1
                continue
            seen_genes.add(gene)
            counts = np.asarray(row[1:], dtype=np.float64)
            total = float(counts.sum())
            if total <= 0:
                genes_skipped_empty += 1
                continue
            genes_seen += 1

            log_cpm = np.log1p(counts * inv_total * CPM_SCALE)
            sum_per_group = np.bincount(
                cell_groups_safe, weights=log_cpm, minlength=n_groups + 1
            )[:n_groups]
            pos = (counts > 0).astype(np.float64)
            pos_per_group = np.bincount(
                cell_groups_safe, weights=pos, minlength=n_groups + 1
            )[:n_groups]

            for gi in range(n_groups):
                n = int(n_per_group[gi])
                if n == 0:
                    continue
                if pos_per_group[gi] == 0:
                    continue
                mean_expr = float(sum_per_group[gi] / n)
                pct = float(100.0 * pos_per_group[gi] / n)
                cell_type, stage = group_keys[gi]
                batch.append(
                    (gene, cell_type, stage, mean_expr, pct, n)
                )

            if len(batch) >= BATCH:
                cur.executemany(INSERT_SQL, batch)
                inserted += len(batch)
                batch.clear()
                if genes_seen % 1000 == 0:
                    rate = genes_seen / (time.perf_counter() - t_stream)
                    print(
                        f"  …{genes_seen:,} genes processed "
                        f"({rate:.1f} genes/s, {inserted:,} rows)",
                        flush=True,
                    )

        if batch:
            cur.executemany(INSERT_SQL, batch)
            inserted += len(batch)

    print(
        f"  ✓ {genes_seen:,} genes processed "
        f"(skipped {genes_skipped_empty:,} empty, "
        f"{duplicates_skipped:,} duplicate symbols) "
        f"→ {inserted:,} rows in {time.perf_counter() - t_stream:.1f} s"
    )

    print("Creating indexes")
    t_idx = time.perf_counter()
    cur.execute("CREATE INDEX idx_fetal_gene ON fetal_heart_expression (gene_symbol)")
    cur.execute(
        "CREATE INDEX idx_fetal_gene_celltype "
        "ON fetal_heart_expression (gene_symbol, cell_type)"
    )
    conn.commit()
    print(f"  ✓ indexes built in {time.perf_counter() - t_idx:.1f} s")

    cur.execute(
        "CREATE TABLE fetal_heart_groups ("
        "cell_type TEXT NOT NULL, "
        "stage TEXT NOT NULL, "
        "n_cells INTEGER NOT NULL, "
        "PRIMARY KEY (cell_type, stage))"
    )
    for gi, (cell_type, stage) in enumerate(group_keys):
        n = int(n_per_group[gi])
        if n == 0:
            continue
        cur.execute(
            "INSERT INTO fetal_heart_groups VALUES (?,?,?)",
            (cell_type, stage, n),
        )

    cur.execute(META_DDL)
    cur.execute(
        "SELECT mean_expr FROM fetal_heart_expression "
        "WHERE mean_expr > 0 ORDER BY mean_expr"
    )
    values = [row[0] for row in cur.fetchall()]
    if values:
        n_val = len(values)
        p50 = values[int(0.50 * n_val)]
        p80 = values[int(0.80 * n_val)]
        p95 = values[int(0.95 * n_val)]
    else:
        p50 = p80 = p95 = 0.0
    cur.executemany(
        "INSERT INTO fetal_heart_meta VALUES (?,?)",
        [
            ("dataset_source", DATASET_LABEL),
            ("portal_url",     "https://cells.ucsc.edu/?ds=hoc"),
            ("band_p50",       f"{p50:.6f}"),
            ("band_p80",       f"{p80:.6f}"),
            ("band_p95",       f"{p95:.6f}"),
        ],
    )
    conn.commit()
    print("Vacuuming")
    t_vac = time.perf_counter()
    conn.execute("VACUUM")
    print(f"  ✓ vacuumed in {time.perf_counter() - t_vac:.1f} s")
    conn.close()

    _dbbuild.publish(tmp_path, db_path)

    print(f"\nBuilt {db_path} in {time.perf_counter() - t_total:.1f} s total")
    return inserted


def main(argv: list[str]) -> int:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    force = "--force-download" in argv

    for url, path in [(META_URL, META_PATH), (EXPR_URL, EXPR_PATH)]:
        if not path.exists() or force:
            download(url, path)
        else:
            size_mb = path.stat().st_size / (1024 * 1024)
            print(
                f"Using cached {path} ({size_mb:.1f} MB) — "
                "pass --force-download to refresh."
            )

    count = build_db(EXPR_PATH, META_PATH, DB_PATH)
    db_size_mb = DB_PATH.stat().st_size / (1024 * 1024)
    print(f"Done. {count:,} pseudobulk rows written to {DB_PATH} ({db_size_mb:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
