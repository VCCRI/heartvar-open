#!/usr/bin/env python3
"""build_gtex_db.py — build the local GTEx cardiac-expression SQLite cache.

Downloads GTEx's public v10 gene-median-TPM GCT
(``GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz``) from the
adult-gtex Google Cloud Storage bucket, extracts the four cardiovascular
tissue columns (``Heart_Left_Ventricle``, ``Heart_Atrial_Appendage``,
``Artery_Aorta`` and ``Artery_Coronary``), and writes one row per
(gene, tissue) into ``data/gtex.db``. The backend GTEx
client reads from this database in place of the live GTEx Portal API
(two HTTP calls per gene: symbol→GENCODE id then median expression).

GTEx releases are static per version, so a one-off build is enough; re-run
on a new GTEx version. Re-running is safe — the table is dropped and
rebuilt every time.

Usage::

    python3 build_gtex_db.py                 # use cached download if present
    python3 build_gtex_db.py --force-download

The cached GCT lives at
``data/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz`` (~8.85 MB);
delete it to force a fresh download on the next run.

Data source: GTEx Portal, dbGaP phs000424 (license NIH-GDS). Please
attribute "GTEx Portal, dbGaP phs000424".
"""

from __future__ import annotations

import argparse
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

URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
    "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
)

USER_AGENT = "HeartVar/1.0"

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DB_PATH = DATA_DIR / "gtex.db"
DOWNLOAD_PATH = DATA_DIR / "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"

TISSUE_HEADER_ALIASES: dict[str, str] = {
    "Heart_Left_Ventricle": "Heart_Left_Ventricle",
    "Heart_Atrial_Appendage": "Heart_Atrial_Appendage",
    "Artery_Aorta": "Artery_Aorta",
    "Artery_Coronary": "Artery_Coronary",
    "Heart - Left Ventricle": "Heart_Left_Ventricle",
    "Heart - Atrial Appendage": "Heart_Atrial_Appendage",
    "Artery - Aorta": "Artery_Aorta",
    "Artery - Coronary": "Artery_Coronary",
}

WANTED_TISSUE_IDS = (
    "Heart_Left_Ventricle",
    "Heart_Atrial_Appendage",
    "Artery_Aorta",
    "Artery_Coronary",
)

DDL = """
CREATE TABLE gtex_cardiac_expression (
    gene_symbol   TEXT,
    gencode_id    TEXT,
    tissue_id     TEXT,
    median_tpm    REAL
)
"""

INSERT_SQL = "INSERT INTO gtex_cardiac_expression VALUES (?,?,?,?)"


def _maybe_float(s: str) -> float | None:
    s = s.strip()
    if not s or s == "-":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _download_urllib(url: str, dest: Path) -> None:
    """Stream the URL into ``dest`` using Python's stdlib (certifi CA
    bundle). Fails on enterprise SSL-inspection networks; ``download()``
    falls back to system curl in that case."""
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
    """Download ``url`` to ``dest``. Tries stdlib first (progress output);
    on TLS verification failure — common on enterprise SSL-inspection
    networks — falls back to system curl."""
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


def build_db(src: Path, db_path: Path) -> int:
    """Stream the gzipped GCT into SQLite. Returns rows inserted.

    GCT layout:
      line 1: "#1.2"
      line 2: "<nrows>\\t<ncols>"
      line 3: header — "Name", "Description", then one column per tissue
      line 4+: data rows (Name = versioned GENCODE id, Description = symbol)
    """
    print(f"Building {db_path} from {src}", flush=True)
    tmp_path = _dbbuild.staging_db_path(db_path)
    conn = _dbbuild.connect(tmp_path)
    conn.execute("PRAGMA journal_mode = OFF")
    conn.execute("PRAGMA synchronous = OFF")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute(DDL)
    cur = conn.cursor()

    BATCH = 20_000
    batch: list[tuple] = []
    inserted = 0
    genes_seen = 0
    t0 = time.perf_counter()

    with gzip.open(src, "rt", encoding="utf-8", newline="") as fh:
        magic = fh.readline().rstrip("\n")
        if not magic.startswith("#1."):
            raise RuntimeError(
                f"unexpected GCT magic line {magic!r}; expected '#1.2'. "
                "GTEx may have changed the file format."
            )
        fh.readline()
        header_line = fh.readline().rstrip("\n")
        header = header_line.split("\t")
        idx = {name: i for i, name in enumerate(header)}

        if "Name" not in idx or "Description" not in idx:
            raise RuntimeError(
                "GCT header missing 'Name'/'Description' columns: "
                f"{header[:5]}. GTEx may have changed the schema."
            )
        name_i = idx["Name"]
        desc_i = idx["Description"]

        tissue_cols: dict[str, int] = {}
        for hdr, col_i in idx.items():
            canonical = TISSUE_HEADER_ALIASES.get(hdr)
            if canonical and canonical not in tissue_cols:
                tissue_cols[canonical] = col_i
        missing = [t for t in WANTED_TISSUE_IDS if t not in tissue_cols]
        if missing:
            raise RuntimeError(
                f"GCT header missing expected cardiac tissue column(s): {missing}. "
                f"Available tissue headers: {sorted(idx)[:10]}…"
            )
        print(
            "  ✓ located cardiac columns: "
            + ", ".join(f"{t}@col{tissue_cols[t]}" for t in WANTED_TISSUE_IDS)
        )

        max_col = max(name_i, desc_i, *tissue_cols.values())
        for line in fh:
            row = line.rstrip("\n").split("\t")
            if len(row) <= max_col:
                continue
            gencode_id = row[name_i].strip()
            gene_symbol = row[desc_i].strip()
            if not gencode_id or not gene_symbol:
                continue
            genes_seen += 1
            for tissue_id in WANTED_TISSUE_IDS:
                tpm = _maybe_float(row[tissue_cols[tissue_id]])
                batch.append((gene_symbol, gencode_id, tissue_id, tpm))
            if len(batch) >= BATCH:
                cur.executemany(INSERT_SQL, batch)
                inserted += len(batch)
                batch.clear()
        if batch:
            cur.executemany(INSERT_SQL, batch)
            inserted += len(batch)

    print(
        f"  ✓ {inserted:,} rows inserted across {genes_seen:,} genes "
        f"({len(WANTED_TISSUE_IDS)} tissues/gene) in {time.perf_counter() - t0:.1f} s"
    )

    print("Creating index")
    t_idx = time.perf_counter()
    cur.execute(
        "CREATE INDEX idx_gtex_gene_upper "
        "ON gtex_cardiac_expression (UPPER(gene_symbol))"
    )
    conn.commit()
    print(f"  ✓ index built in {time.perf_counter() - t_idx:.1f} s")

    conn.close()

    _dbbuild.publish(tmp_path, db_path)
    return inserted


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="re-download the GCT even if a cached copy exists",
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
    print(f"\nDone. {count:,} rows written to {DB_PATH} ({db_size_mb:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
