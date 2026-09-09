#!/usr/bin/env python3
"""build_gnomad_constraint_db.py — build the local gnomAD gene-constraint cache.

Downloads gnomAD v4.1's public per-transcript constraint metrics TSV and
writes one row per gene (the MANE Select transcript, falling back to the
canonical transcript) into ``data/gnomad_constraint.db``. The backend gnomAD
client serves the gene-level ``gnomad_constraint`` block (pLI / LOEUF / oe_mis
/ mis_z / syn_z) from this database in place of the live gnomAD GraphQL
``gene`` query.

The file is tens of MB and needs no auth. gnomAD release constraint metrics
are static per release, so a one-off build is enough; re-run on a new gnomAD
release. Re-running is safe — the table is dropped and rebuilt every time.

Usage::

    python3 build_gnomad_constraint_db.py                 # use cached download
    python3 build_gnomad_constraint_db.py --force-download

The cached TSV lives at ``data/gnomad.v4.1.constraint_metrics.tsv``; delete it
to force a fresh download on the next run.

Column mapping (TSV header -> GraphQL gnomad_constraint field), located BY
HEADER NAME at build time (never by position):

    pLI            <- lof.pLI
    oe_lof         <- lof.oe
    oe_lof_upper   <- lof.oe_ci.upper    (LOEUF)
    oe_mis         <- mis.oe
    oe_mis_upper   <- mis.oe_ci.upper
    mis_z          <- mis.z_score
    syn_z          <- syn.z_score

Transcript selection per gene: the gnomAD GraphQL ``gene`` query keys on the
Ensembl gene/transcript model and returns ``gene_id`` as ``ENSG…``. The
constraint TSV interleaves Ensembl rows (``gene_id`` = ``ENSG…``, transcript
``ENST…``) with RefSeq rows (``gene_id`` = Entrez int, transcript ``NM_…``) for
the same symbol, so we (1) prefer the Ensembl ``ENSG…`` row, then (2) within
that, the ``mane_select`` row, then ``canonical``, then first-seen. This
reproduces the GraphQL ``gene.gene_id`` (ENSG) and the MANE/canonical
transcript's metrics.

Data source: gnomAD v4.1 (https://gnomad.broadinstitute.org), released under
CC0. Please attribute "gnomAD v4.1".
"""

from __future__ import annotations

import argparse
import csv
import shutil
import ssl
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import _dbbuild

URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/"
    "constraint/gnomad.v4.1.constraint_metrics.tsv"
)

USER_AGENT = "HeartVar/1.0"

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DB_PATH = DATA_DIR / "gnomad_constraint.db"
DOWNLOAD_PATH = DATA_DIR / "gnomad.v4.1.constraint_metrics.tsv"

COL_MAP = {
    "lof.pLI": "pLI",
    "lof.oe": "oe_lof",
    "lof.oe_ci.upper": "oe_lof_upper",
    "mis.oe": "oe_mis",
    "mis.oe_ci.upper": "oe_mis_upper",
    "mis.z_score": "mis_z",
    "syn.z_score": "syn_z",
}
GENE_COL = "gene"
GENE_ID_COL = "gene_id"
CANONICAL_COL = "canonical"
MANE_COL = "mane_select"

DDL = """
CREATE TABLE gnomad_constraint (
    gene_symbol     TEXT,
    gene_id         TEXT,
    pLI             REAL,
    oe_lof          REAL,
    oe_lof_upper    REAL,
    oe_mis          REAL,
    oe_mis_upper    REAL,
    mis_z           REAL,
    syn_z           REAL
)
"""

INSERT_SQL = "INSERT INTO gnomad_constraint VALUES (?,?,?,?,?,?,?,?,?)"


def _maybe_float(s: str | None) -> float | None:
    if s is None:
        return None
    s = s.strip()
    if not s or s in ("-", "NA", "NaN", "nan"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _is_true(s: str | None) -> bool:
    """gnomAD constraint TSV writes booleans as 'true'/'false' (and the row
    may be blank for non-canonical/non-MANE transcripts)."""
    return (s or "").strip().lower() == "true"


def _download_urllib(url: str, dest: Path) -> None:
    """Stream the URL into ``dest`` using Python's stdlib (certifi CA bundle).
    Fails on enterprise SSL-inspection networks; ``download()`` falls back to
    system curl in that case."""
    t0 = time.perf_counter()
    bytes_seen = 0
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req) as r, open(dest, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            bytes_seen += len(chunk)
            if bytes_seen % (1 << 22) < (1 << 20):
                print(f"  …{bytes_seen / (1024 * 1024):.1f} MB", flush=True)
    size_mb = dest.stat().st_size / (1024 * 1024)
    print(f"  ✓ {size_mb:.1f} MB in {time.perf_counter() - t0:.1f} s")


def _download_curl(url: str, dest: Path) -> None:
    """Fallback path: shell out to system curl (uses the macOS keychain on
    darwin, so enterprise SSL-inspection roots are trusted automatically)."""
    curl = shutil.which("curl")
    if not curl:
        raise RuntimeError(
            "Python's TLS chain rejected the server certificate (likely an "
            "enterprise inspection proxy) and `curl` isn't on PATH. Install "
            "curl, or add the corporate root CA to certifi's bundle."
        )
    print(f"  …falling back to {curl} (uses system keychain)", flush=True)
    t0 = time.perf_counter()
    cmd = [
        curl, "--location", "--fail", "--retry", "3",
        "--silent", "--show-error",
        "--user-agent", USER_AGENT,
        "--output", str(dest), url,
    ]
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise RuntimeError(f"curl exited with status {proc.returncode}")
    size_mb = dest.stat().st_size / (1024 * 1024)
    print(f"  ✓ {size_mb:.1f} MB in {time.perf_counter() - t0:.1f} s")


def download(url: str, dest: Path) -> None:
    """Download ``url`` to ``dest``. Tries stdlib first (progress output); on
    TLS verification failure — common on enterprise SSL-inspection networks —
    falls back to system curl."""
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


def _transcript_rank(gene_id: str, mane: bool, canonical: bool) -> int:
    """Selection priority for a gene's transcript rows (lower wins).

    The gnomAD GraphQL ``gene`` query returns the Ensembl gene model
    (``gene_id`` = ``ENSG…``), so an Ensembl row always outranks a RefSeq row
    for the same symbol. Within each model, MANE Select > canonical > other.
    Ranks: ENSG+MANE 0, ENSG+canonical 1, ENSG other 2, RefSeq+MANE 3,
    RefSeq+canonical 4, RefSeq other 5.
    """
    ensembl = (gene_id or "").strip().upper().startswith("ENSG")
    base = 0 if ensembl else 3
    if mane:
        return base + 0
    if canonical:
        return base + 1
    return base + 2


def build_db(src: Path, db_path: Path) -> int:
    """Stream the constraint TSV into SQLite, one row per gene (MANE Select /
    canonical transcript). Returns gene rows inserted."""
    print(f"Building {db_path} from {src}", flush=True)
    tmp_path = _dbbuild.staging_db_path(db_path)

    best: dict[str, tuple[int, tuple]] = {}
    t0 = time.perf_counter()
    seen_rows = 0

    with src.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        required = [GENE_COL, GENE_ID_COL, CANONICAL_COL, MANE_COL, *COL_MAP]
        missing = [c for c in required if c not in idx]
        if missing:
            raise RuntimeError(
                f"gnomAD constraint TSV header missing expected columns: "
                f"{missing}. gnomAD may have changed the schema. "
                f"First 10 headers: {header[:10]}"
            )
        max_col = max(idx[c] for c in required)

        for row in reader:
            if len(row) <= max_col:
                continue
            gene_symbol = (row[idx[GENE_COL]] or "").strip()
            if not gene_symbol:
                continue
            seen_rows += 1
            gene_id = (row[idx[GENE_ID_COL]] or "").strip()
            mane = _is_true(row[idx[MANE_COL]])
            canonical = _is_true(row[idx[CANONICAL_COL]])
            rank = _transcript_rank(gene_id, mane, canonical)
            prior = best.get(gene_symbol)
            if prior is not None and prior[0] <= rank:
                continue
            rec = (
                gene_symbol,
                gene_id or None,
                _maybe_float(row[idx["lof.pLI"]]),
                _maybe_float(row[idx["lof.oe"]]),
                _maybe_float(row[idx["lof.oe_ci.upper"]]),
                _maybe_float(row[idx["mis.oe"]]),
                _maybe_float(row[idx["mis.oe_ci.upper"]]),
                _maybe_float(row[idx["mis.z_score"]]),
                _maybe_float(row[idx["syn.z_score"]]),
            )
            best[gene_symbol] = (rank, rec)

    conn = _dbbuild.connect(tmp_path)
    conn.execute("PRAGMA journal_mode = OFF")
    conn.execute("PRAGMA synchronous = OFF")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute(DDL)
    cur = conn.cursor()
    rows = [rec for (_rank, rec) in best.values()]
    cur.executemany(INSERT_SQL, rows)
    inserted = len(rows)
    print(
        f"  ✓ {inserted:,} gene rows inserted "
        f"(from {seen_rows:,} transcript rows) in {time.perf_counter() - t0:.1f} s"
    )

    print("Creating index")
    t_idx = time.perf_counter()
    cur.execute(
        "CREATE INDEX idx_gnomad_constraint_gene_upper "
        "ON gnomad_constraint (UPPER(gene_symbol))"
    )
    conn.commit()
    print(f"  ✓ index built in {time.perf_counter() - t_idx:.1f} s")

    for sym in ("MYH7", "TTN"):
        r = cur.execute(
            "SELECT gene_symbol, gene_id, pLI, oe_lof_upper, mis_z "
            "FROM gnomad_constraint WHERE UPPER(gene_symbol) = ?",
            (sym,),
        ).fetchone()
        print(f"  spot-check {sym}: {r}")

    conn.close()

    _dbbuild.publish(tmp_path, db_path)
    return inserted


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="re-download the constraint TSV even if a cached copy exists",
    )
    args = parser.parse_args(argv)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    if not DOWNLOAD_PATH.exists() or args.force_download:
        download(URL, DOWNLOAD_PATH)
    else:
        size_mb = DOWNLOAD_PATH.stat().st_size / (1024 * 1024)
        print(
            f"Using cached {DOWNLOAD_PATH} ({size_mb:.1f} MB) — "
            "pass --force-download to refresh."
        )
    count = build_db(DOWNLOAD_PATH, DB_PATH)
    db_size_mb = DB_PATH.stat().st_size / (1024 * 1024)
    print(f"\nDone. {count:,} gene rows written to {DB_PATH} ({db_size_mb:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
