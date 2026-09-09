#!/usr/bin/env python3
"""build_uniprot_db.py — build the local UniProt SQLite cache.

Streams UniProt's reviewed human proteome (SwissProt, ~20k entries) from
the public REST ``/uniprotkb/stream`` endpoint with the exact field set
the live HeartVar client requests, runs every entry through the SAME
parsing path the live client uses (``backend.clients.uniprot.parse_entry``),
and stores the resulting structured dict as JSON in ``data/uniprot.db``.

The backend UniProt client (``backend/clients/uniprot.py``) reads from this
database in place of a live ``/uniprotkb/search`` call, falling back to the
live call when the DB is absent or a gene isn't found locally.

UniProt releases roughly every eight weeks; re-running this monthly keeps
the local cache fresh for a curation aid. Re-running is safe — the
``uniprot_entry`` table is dropped and rebuilt every time.

Each entry is indexed under its primary gene symbol AND every synonym
(uppercased), mirroring the live client's ``gene_exact`` matching. When a
symbol maps to several entries the FIRST one streamed wins (mirroring the
live client's ``size=1`` "first match" behaviour); later collisions are
skipped.

Data: UniProt Knowledgebase, licensed CC-BY 4.0 (attribute UniProt).

Usage::

    python3 scripts/build_uniprot_db.py                 # use cached download if present
    python3 scripts/build_uniprot_db.py --force-download

The cached JSON download lives at ``data/uniprot_human_reviewed.json``;
delete it (or pass ``--force-download``) to refresh on the next run.
"""

from __future__ import annotations

import argparse
import json
import shutil
import ssl
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import _dbbuild

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from backend.clients.uniprot import (  # noqa: E402
    UNIPROT_FIELDS,
    gene_symbols,
    parse_entry,
)

CONTACT_UA = "HeartVar/1.0"

STREAM_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=organism_id:9606+AND+reviewed:true"
    "&format=json"
    "&fields=" + ",".join(UNIPROT_FIELDS)
)

DATA_DIR = PROJECT_ROOT / "data"
DB_PATH = DATA_DIR / "uniprot.db"
DOWNLOAD_PATH = DATA_DIR / "uniprot_human_reviewed.json"

DDL = """
CREATE TABLE uniprot_entry (
    gene_symbol   TEXT,
    accession     TEXT,
    payload       TEXT
)
"""

INSERT_SQL = "INSERT INTO uniprot_entry (gene_symbol, accession, payload) VALUES (?,?,?)"


def _download_urllib(url: str, dest: Path) -> None:
    """Stream the URL into ``dest`` using Python's stdlib. Uses certifi's
    CA bundle, so it fails on networks with enterprise SSL inspection where
    the intercept root isn't in certifi. ``download()`` falls back to system
    ``curl`` (which trusts the macOS keychain) in that case."""
    t0 = time.perf_counter()
    bytes_seen = 0
    req = urllib.request.Request(url, headers={"User-Agent": CONTACT_UA})
    with urllib.request.urlopen(req) as r, open(dest, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            bytes_seen += len(chunk)
            if bytes_seen % (1 << 24) < (1 << 20):
                print(f"  …{bytes_seen / (1024 * 1024):.1f} MB", flush=True)
    size_mb = dest.stat().st_size / (1024 * 1024)
    print(f"  ✓ {size_mb:.1f} MB in {time.perf_counter() - t0:.1f} s")


def _download_curl(url: str, dest: Path) -> None:
    """Fallback path: shell out to system curl. On macOS curl uses Secure
    Transport / the system keychain, so enterprise SSL-inspection roots are
    trusted automatically without extra configuration."""
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
        "--user-agent", CONTACT_UA,
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
    print(f"Downloading reviewed human proteome\n  → {dest}", flush=True)
    print("  (this can take a minute or more)", flush=True)
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


def _iter_entries(src: Path):
    """Yield UniProtKB entries from the streamed JSON download.

    The ``/uniprotkb/stream`` endpoint with ``format=json`` returns a
    single JSON object ``{"results": [ ...entries... ]}``. The whole file
    is low tens of MB, so a one-shot ``json.load`` is fine here.
    """
    with open(src, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    results = data.get("results")
    if results is None:
        raise RuntimeError(
            "UniProt stream download has no 'results' array — the download may "
            "be truncated or the API response shape changed. Re-run with "
            "--force-download."
        )
    yield from results


def build_db(src: Path, db_path: Path) -> tuple[int, int]:
    """Parse the streamed proteome into SQLite.

    Returns ``(entries_inserted, index_rows)``. ``entries_inserted`` is the
    count of distinct UniProt entries stored; ``index_rows`` is the total
    gene-symbol index rows (primary + synonyms, after collision skips).
    """
    print(f"Building {db_path} from {src}", flush=True)
    tmp_path = _dbbuild.staging_db_path(db_path)
    conn = _dbbuild.connect(tmp_path)
    conn.execute("PRAGMA journal_mode = OFF")
    conn.execute("PRAGMA synchronous = OFF")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute(DDL)
    cur = conn.cursor()

    inserted = 0
    index_rows = 0
    skipped_no_symbol = 0
    collisions = 0
    claimed: set[str] = set()
    batch: list[tuple[str, str, str]] = []
    BATCH = 10_000
    t0 = time.perf_counter()

    for entry in _iter_entries(src):
        primary, synonyms = gene_symbols(entry)
        symbols: list[str] = []
        seen_upper: set[str] = set()
        for sym in ([primary] if primary else []) + synonyms:
            up = sym.upper()
            if up in seen_upper:
                continue
            seen_upper.add(up)
            symbols.append(up)
        if not symbols:
            skipped_no_symbol += 1
            continue

        payload = parse_entry(entry)
        payload.pop("gene", None)
        accession = payload.get("accession")
        payload_json = json.dumps(payload, separators=(",", ":"))

        any_new = False
        for up in symbols:
            if up in claimed:
                collisions += 1
                continue
            claimed.add(up)
            batch.append((up, accession, payload_json))
            index_rows += 1
            any_new = True
        if any_new:
            inserted += 1
        if len(batch) >= BATCH:
            cur.executemany(INSERT_SQL, batch)
            batch.clear()

    if batch:
        cur.executemany(INSERT_SQL, batch)

    print(
        f"  ✓ {inserted:,} entries / {index_rows:,} symbol rows "
        f"(skipped {skipped_no_symbol:,} symbol-less, "
        f"{collisions:,} symbol collisions) "
        f"in {time.perf_counter() - t0:.1f} s"
    )

    print("Creating index")
    t_idx = time.perf_counter()
    cur.execute("CREATE INDEX idx_uniprot_gene ON uniprot_entry (gene_symbol)")
    conn.commit()
    print(f"  ✓ index built in {time.perf_counter() - t_idx:.1f} s")

    conn.close()
    _dbbuild.publish(tmp_path, db_path)
    return inserted, index_rows


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="re-download the proteome even if the cached file is present",
    )
    args = parser.parse_args(argv)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    if not DOWNLOAD_PATH.exists() or args.force_download:
        download(STREAM_URL, DOWNLOAD_PATH)
    else:
        size_mb = DOWNLOAD_PATH.stat().st_size / (1024 * 1024)
        print(
            f"Using cached {DOWNLOAD_PATH} ({size_mb:.1f} MB) — "
            "pass --force-download to refresh."
        )
    inserted, index_rows = build_db(DOWNLOAD_PATH, DB_PATH)
    db_size_mb = DB_PATH.stat().st_size / (1024 * 1024)
    print(
        f"\nDone. {inserted:,} entries ({index_rows:,} symbol rows) written to "
        f"{DB_PATH} ({db_size_mb:.1f} MB)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
