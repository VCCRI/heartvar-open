#!/usr/bin/env python3
"""build_clinvar_db.py — build the local ClinVar SQLite cache.

Downloads ClinVar's public ``variant_summary.txt.gz`` from NCBI's FTP,
parses the GRCh38 rows, and writes them into ``data/clinvar.db``. The
backend ClinVar client and the ``/api/clinvar`` route read from this
database in place of the live E-utilities calls.

ClinVar publishes ``variant_summary`` weekly; running this script monthly
keeps the local cache reasonably fresh for a curation aid. Re-running is
safe — the variants table is dropped and rebuilt every time.

Usage::

    python3 build_clinvar_db.py                 # use cached download if present
    python3 build_clinvar_db.py --force-download

The cached gzip lives at ``data/variant_summary.txt.gz``; delete it to
force a fresh download on the next run.
"""

from __future__ import annotations

import csv
import gzip
import shutil
import ssl
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import _dbbuild

URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DB_PATH = DATA_DIR / "clinvar.db"
DOWNLOAD_PATH = DATA_DIR / "variant_summary.txt.gz"

WANTED_COLS = [
    "VariationID",
    "Type",
    "Name",
    "GeneSymbol",
    "ClinicalSignificance",
    "ReviewStatus",
    "NumberSubmitters",
    "PhenotypeList",
    "ChromosomeAccession",
    "Start",
    "Stop",
    "ReferenceAllele",
    "AlternateAllele",
    "Assembly",
    "LastEvaluated",
]

DDL = """
CREATE TABLE variants (
    variation_id            INTEGER,
    obj_type                TEXT,
    name                    TEXT,
    gene_symbol             TEXT,
    clinical_significance   TEXT,
    review_status           TEXT,
    number_submitters       INTEGER,
    phenotype_list          TEXT,
    chromosome_accession    TEXT,
    start                   INTEGER,
    stop                    INTEGER,
    reference_allele        TEXT,
    alternate_allele        TEXT,
    assembly                TEXT,
    last_evaluated          TEXT
)
"""

INSERT_SQL = (
    "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
)


def _maybe_int(s: str) -> int | None:
    s = s.strip()
    if not s or s == "-":
        return None
    try:
        return int(s)
    except ValueError:
        return None


def _download_urllib(url: str, dest: Path) -> None:
    """Stream the URL into ``dest`` using Python's stdlib. Uses certifi's
    CA bundle, so it fails on networks with enterprise SSL inspection
    where the intercept root isn't in certifi. ``download()`` falls back
    to system ``curl`` (which trusts the macOS keychain) in that case.
    """
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
    """Fallback path: shell out to system curl. On macOS curl uses Secure
    Transport / the system keychain, so enterprise SSL-inspection roots
    are trusted automatically without extra configuration."""
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
        "--output", str(dest), url,
    ]
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise RuntimeError(f"curl exited with status {proc.returncode}")
    size_mb = dest.stat().st_size / (1024 * 1024)
    print(f"  ✓ {size_mb:.1f} MB in {time.perf_counter() - t0:.1f} s")


def download(url: str, dest: Path) -> None:
    """Download ``url`` to ``dest``. Tries stdlib first (better progress
    output); on TLS verification failure — common on networks with
    enterprise SSL inspection — falls back to system curl."""
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
    """Stream the gzipped TSV into SQLite. Returns rows inserted."""
    print(f"Building {db_path} from {src}", flush=True)
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
    skipped_assembly = 0
    skipped_malformed = 0
    t0 = time.perf_counter()

    with gzip.open(src, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if header and header[0].startswith("#"):
            header[0] = header[0].lstrip("#")
        idx = {name: i for i, name in enumerate(header)}
        missing = [c for c in WANTED_COLS if c not in idx]
        if missing:
            raise RuntimeError(
                f"variant_summary header missing expected columns: {missing}. "
                "ClinVar may have changed the schema."
            )

        for row in reader:
            try:
                assembly = row[idx["Assembly"]]
            except IndexError:
                skipped_malformed += 1
                continue
            if assembly != "GRCh38":
                skipped_assembly += 1
                continue
            try:
                rec = (
                    _maybe_int(row[idx["VariationID"]]),
                    row[idx["Type"]],
                    row[idx["Name"]],
                    row[idx["GeneSymbol"]],
                    row[idx["ClinicalSignificance"]],
                    row[idx["ReviewStatus"]],
                    _maybe_int(row[idx["NumberSubmitters"]]) or 0,
                    row[idx["PhenotypeList"]],
                    row[idx["ChromosomeAccession"]],
                    _maybe_int(row[idx["Start"]]),
                    _maybe_int(row[idx["Stop"]]),
                    row[idx["ReferenceAllele"]],
                    row[idx["AlternateAllele"]],
                    assembly,
                    row[idx["LastEvaluated"]],
                )
            except IndexError:
                skipped_malformed += 1
                continue
            batch.append(rec)
            if len(batch) >= BATCH:
                cur.executemany(INSERT_SQL, batch)
                inserted += len(batch)
                batch.clear()
                if inserted % 500_000 == 0:
                    print(
                        f"  …{inserted:,} rows ({time.perf_counter() - t0:.1f} s)",
                        flush=True,
                    )
        if batch:
            cur.executemany(INSERT_SQL, batch)
            inserted += len(batch)

    print(
        f"  ✓ {inserted:,} GRCh38 rows inserted "
        f"(skipped {skipped_assembly:,} non-GRCh38, "
        f"{skipped_malformed:,} malformed) "
        f"in {time.perf_counter() - t0:.1f} s"
    )

    print("Creating indexes")
    t_idx = time.perf_counter()
    cur.execute("PRAGMA temp_store=FILE")
    cur.execute(f"PRAGMA temp_store_directory='{tmp_path.parent}'")
    cur.execute("CREATE INDEX idx_variants_gene ON variants (gene_symbol)")
    cur.execute("CREATE INDEX idx_variants_name ON variants (name)")
    cur.execute("CREATE INDEX idx_variants_variation_id ON variants (variation_id)")
    cur.execute(
        "CREATE INDEX idx_variants_gene_covering ON variants ("
        "gene_symbol, number_submitters DESC, variation_id ASC, "
        "name, clinical_significance, review_status, phenotype_list, "
        "start, stop, obj_type, last_evaluated)"
    )
    conn.commit()
    built = [r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='index' "
        "AND tbl_name='variants' AND name NOT LIKE 'sqlite_%' ORDER BY name")]
    print(f"  ✓ indexes built in {time.perf_counter() - t_idx:.1f} s: "
          f"{', '.join(built)}")

    required = {
        "idx_variants_gene",
        "idx_variants_name",
        "idx_variants_variation_id",
        "idx_variants_gene_covering",
    }
    missing = sorted(required - set(built))
    if missing:
        raise SystemExit(
            f"clinvar build is missing required index/indexes: {missing}. "
            f"Built: {built}. Refusing to publish — the app would be slow with "
            "no way to tell from the data."
        )

    conn.close()
    _dbbuild.publish(tmp_path, db_path)
    return inserted


def main(argv: list[str]) -> int:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    if not DOWNLOAD_PATH.exists() or "--force-download" in argv:
        download(URL, DOWNLOAD_PATH)
    else:
        size_mb = DOWNLOAD_PATH.stat().st_size / (1024 * 1024)
        print(
            f"Using cached {DOWNLOAD_PATH} ({size_mb:.1f} MB) — "
            "pass --force-download to refresh."
        )
    count = build_db(DOWNLOAD_PATH, DB_PATH)
    db_size_mb = DB_PATH.stat().st_size / (1024 * 1024)
    print(f"\nDone. {count:,} records written to {DB_PATH} ({db_size_mb:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
