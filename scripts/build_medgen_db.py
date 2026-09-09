#!/usr/bin/env python3
"""build_medgen_db.py — build the local MedGen gene→condition cache.

Downloads NCBI ClinVar's public ``gene_condition_source_id`` table and
writes its gene→condition associations into ``data/medgen.db``. The
backend MedGen client (``backend/clients/medgen.py``) reads from this
database in place of the three live NCBI E-utilities calls (esearch +
elink + esummary), freeing the shared 3 req/s NCBI budget the literature
clients also depend on. The live E-utilities path is kept as a graceful
fallback for genes absent from the local table.

``gene_condition_source_id`` is one row per gene–disease–source
association (a single condition can appear from several sources, e.g.
MONDO + OMIM + GeneReviews). NCBI refreshes the file DAILY; running this
script monthly keeps the local cache reasonably fresh for a curation
aid. Re-running is safe — the conditions table is dropped and rebuilt
every time.

Usage::

    python3 scripts/build_medgen_db.py                 # use cached download if present
    python3 scripts/build_medgen_db.py --force-download

The cached TSV lives at ``data/gene_condition_source_id``; delete it (or
pass ``--force-download``) to force a fresh download on the next run.

License: NCBI public domain. Only the integer OMIM MIM cross-reference is
used to build an omim.org URL — no OMIM text is fetched, cached, or
logged. NLM acknowledgment requested.
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

URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id"

USER_AGENT = "HeartVar/1.0"

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DB_PATH = DATA_DIR / "medgen.db"
DOWNLOAD_PATH = DATA_DIR / "gene_condition_source_id"

COL_GENE = "AssociatedGenes"
COL_CUI = "ConceptID"
COL_MIM = "DiseaseMIM"
COL_NAME = "DiseaseName"
COL_SOURCE = "SourceName"
COL_LASTUPDATED = "LastUpdated"
WANTED_COLS = [COL_GENE, COL_CUI, COL_MIM, COL_NAME, COL_SOURCE, COL_LASTUPDATED]

DDL = """
CREATE TABLE conditions (
    gene_symbol   TEXT,
    cui           TEXT,
    mim           TEXT,
    name          TEXT,
    source        TEXT,
    last_updated  TEXT
)
"""

INSERT_SQL = "INSERT INTO conditions VALUES (?,?,?,?,?,?)"


def _clean(s: str | None) -> str | None:
    """Strip whitespace; treat empty / NCBI '-' placeholder as None."""
    if s is None:
        return None
    s = s.strip()
    if not s or s == "-":
        return None
    return s


def _download_urllib(url: str, dest: Path) -> None:
    """Stream the URL into ``dest`` using Python's stdlib. Uses certifi's
    CA bundle, so it fails on networks with enterprise SSL inspection
    where the intercept root isn't in certifi. ``download()`` falls back
    to system ``curl`` (which trusts the macOS keychain) in that case."""
    t0 = time.perf_counter()
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    bytes_seen = 0
    with urllib.request.urlopen(req) as r, open(dest, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            bytes_seen += len(chunk)
    size_mb = dest.stat().st_size / (1024 * 1024)
    print(f"  ✓ {size_mb:.2f} MB in {time.perf_counter() - t0:.1f} s")


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
        "--user-agent", USER_AGENT,
        "--output", str(dest), url,
    ]
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise RuntimeError(f"curl exited with status {proc.returncode}")
    size_mb = dest.stat().st_size / (1024 * 1024)
    print(f"  ✓ {size_mb:.2f} MB in {time.perf_counter() - t0:.1f} s")


def download(url: str, dest: Path) -> None:
    """Download ``url`` to ``dest``. Tries stdlib first; on TLS
    verification failure — common on networks with enterprise SSL
    inspection — falls back to system curl."""
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
    """Stream the TSV into SQLite. Returns rows inserted.

    Rows without an AssociatedGenes symbol are skipped (the gene→condition
    key is undefined there — those rows carry only RelatedGenes). Rows
    without a ConceptID are skipped too (no MedGen CUI to surface).

    The source file lists a single condition from multiple sources
    (MONDO, OMIM, GeneReviews, HPO, …). We dedupe to ONE row per
    ``(UPPER(gene), cui)`` so a gene's condition isn't listed several
    times — mirroring the live client, where each MedGen UID is one
    condition. The first row seen (file order) wins, but if a later
    duplicate carries a MIM cross-reference and the kept row doesn't, the
    MIM is back-filled so the OMIM URL stays populated. File order is
    preserved on insert so the client's cardiac-first stable sort matches
    the live ordering as closely as the source allows.
    """
    print(f"Building {db_path} from {src}", flush=True)
    tmp_path = _dbbuild.staging_db_path(db_path)
    conn = _dbbuild.connect(tmp_path)
    conn.execute("PRAGMA journal_mode = OFF")
    conn.execute("PRAGMA synchronous = OFF")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute(DDL)
    cur = conn.cursor()

    t0 = time.perf_counter()
    skipped_no_gene = 0
    skipped_no_cui = 0
    skipped_malformed = 0

    seen: dict[tuple[str, str], list] = {}

    with open(src, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if header and header[0].startswith("#"):
            header[0] = header[0].lstrip("#")
        idx = {name: i for i, name in enumerate(header)}
        missing = [c for c in WANTED_COLS if c not in idx]
        if missing:
            raise RuntimeError(
                f"gene_condition_source_id header missing expected columns: "
                f"{missing}. NCBI may have changed the schema."
            )

        for row in reader:
            try:
                gene = _clean(row[idx[COL_GENE]])
                cui = _clean(row[idx[COL_CUI]])
                mim = _clean(row[idx[COL_MIM]])
                name = _clean(row[idx[COL_NAME]])
                source = _clean(row[idx[COL_SOURCE]])
                last_updated = _clean(row[idx[COL_LASTUPDATED]])
            except IndexError:
                skipped_malformed += 1
                continue
            if not gene:
                skipped_no_gene += 1
                continue
            if not cui:
                skipped_no_cui += 1
                continue
            key = (gene.upper(), cui)
            existing = seen.get(key)
            if existing is None:
                seen[key] = [gene, cui, mim, name, source, last_updated]
            elif existing[2] is None and mim is not None:
                existing[2] = mim

    rows = [tuple(v) for v in seen.values()]
    cur.executemany(INSERT_SQL, rows)
    inserted = len(rows)
    print(
        f"  ✓ {inserted:,} unique gene→condition rows inserted "
        f"(skipped {skipped_no_gene:,} no-AssociatedGenes, "
        f"{skipped_no_cui:,} no-ConceptID, "
        f"{skipped_malformed:,} malformed) "
        f"in {time.perf_counter() - t0:.1f} s"
    )

    print("Creating index")
    t_idx = time.perf_counter()
    cur.execute(
        "CREATE INDEX idx_conditions_gene ON conditions (UPPER(gene_symbol))"
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
        help="re-download the bulk file even if a cached copy exists",
    )
    args = parser.parse_args(argv)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    if not DOWNLOAD_PATH.exists() or args.force_download:
        download(URL, DOWNLOAD_PATH)
    else:
        size_mb = DOWNLOAD_PATH.stat().st_size / (1024 * 1024)
        print(
            f"Using cached {DOWNLOAD_PATH} ({size_mb:.2f} MB) — "
            "pass --force-download to refresh."
        )
    count = build_db(DOWNLOAD_PATH, DB_PATH)
    db_size_mb = DB_PATH.stat().st_size / (1024 * 1024)
    print(f"\nDone. {count:,} records written to {DB_PATH} ({db_size_mb:.2f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
