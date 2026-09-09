#!/usr/bin/env python3
"""build_biogrid_db.py — build the local BioGRID interactions cache.

Downloads BioGRID's bulk tab3 archive for Homo sapiens, filters to
physical interactions only, and loads them into ``data/biogrid.db``.

BioGRID interaction data is freely available for academic use — see
https://thebiogrid.org/ — and the bulk download requires no licence
or API key.

BioGRID publishes a new release roughly monthly. Re-run this script
**quarterly** to keep the local cache reasonably fresh for a curation
aid. The interactions table is dropped and rebuilt every run, so the
script is safely re-runnable.

Usage::

    python3 build_biogrid_db.py                 # use cached download if present
    python3 build_biogrid_db.py --force-download
"""

from __future__ import annotations

import csv
import shutil
import ssl
import subprocess
import sys
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

import _dbbuild

URL = (
    "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/"
    "BIOGRID-ORGANISM-LATEST.tab3.zip"
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DB_PATH = DATA_DIR / "biogrid.db"
DOWNLOAD_PATH = DATA_DIR / "biogrid-organism-latest-tab3.zip"

WANTED_COLS = [
    "Official Symbol Interactor A",
    "Official Symbol Interactor B",
    "Experimental System",
    "Experimental System Type",
    "Author",
    "Publication Source",
    "Throughput",
]

DDL = """
CREATE TABLE interactions (
    symbol_a                  TEXT NOT NULL,
    symbol_b                  TEXT NOT NULL,
    experimental_system       TEXT,
    experimental_system_type  TEXT,
    author                    TEXT,
    pmid                      TEXT,
    throughput                TEXT
)
"""

INSERT_SQL = "INSERT INTO interactions VALUES (?,?,?,?,?,?,?)"


def _download_urllib(url: str, dest: Path) -> None:
    """Stream a URL into ``dest`` using the stdlib (certifi CA bundle).

    Sends a browser-like User-Agent: BioGRID's download host returns
    ``HTTP 403 Forbidden`` to the default ``Python-urllib/x.y`` agent, so
    the bare urlopen fails even though the file is public. The header makes
    the request look like an ordinary client and is accepted."""
    t0 = time.perf_counter()
    bytes_seen = 0
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": (
                "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"
            )
        },
    )
    with urllib.request.urlopen(req) as r, open(dest, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            bytes_seen += len(chunk)
            if bytes_seen % (1 << 25) < (1 << 20):
                print(f"  …{bytes_seen / (1024 * 1024):.1f} MB", flush=True)
    size_mb = dest.stat().st_size / (1024 * 1024)
    print(f"  ✓ {size_mb:.1f} MB in {time.perf_counter() - t0:.1f} s")


def _download_curl(url: str, dest: Path) -> None:
    """Fallback for networks with enterprise SSL inspection: shell out to
    system ``curl``, which uses Secure Transport / the macOS keychain
    and so trusts the corporate intercept root automatically."""
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
    """Download ``url`` to ``dest``, falling back to curl on TLS failure."""
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


def _pmid_from_pubsource(raw: str) -> str:
    """BioGRID stores ``Publication Source`` as ``PUBMED:12345``. Strip
    the prefix if present; pass through unchanged otherwise (a small
    number of records carry DOI: or other prefixes)."""
    if not raw:
        return ""
    s = raw.strip()
    if s.upper().startswith("PUBMED:"):
        return s.split(":", 1)[1].strip()
    return s


def _pick_tab3_member(zf: zipfile.ZipFile) -> str:
    """Find the actual tab3 .txt file inside the zip (the filename
    embeds the BioGRID release number so we can't hard-code it)."""
    candidates = [
        n for n in zf.namelist()
        if n.lower().endswith(".tab3.txt") and "homo_sapiens" in n.lower()
    ]
    if not candidates:
        raise RuntimeError(
            f"No *.tab3.txt member found in {zf.filename!r}. "
            f"Archive contents: {zf.namelist()}"
        )
    return candidates[0]


def build_db(src: Path, db_path: Path) -> int:
    """Parse the zipped tab3 file, keep physical interactions, load
    into SQLite. Returns rows inserted."""
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
    skipped_genetic = 0
    skipped_missing_symbol = 0
    t0 = time.perf_counter()

    with zipfile.ZipFile(src) as zf:
        member = _pick_tab3_member(zf)
        print(f"  reading {member}", flush=True)
        with zf.open(member) as raw:
            text = (line.decode("utf-8", errors="replace") for line in raw)
            reader = csv.reader(text, delimiter="\t")
            header = next(reader)
            if header and header[0].startswith("#"):
                header[0] = header[0].lstrip("#").strip()
            idx = {name: i for i, name in enumerate(header)}
            missing = [c for c in WANTED_COLS if c not in idx]
            if missing:
                raise RuntimeError(
                    f"tab3 header missing expected columns: {missing}. "
                    "BioGRID may have changed the schema."
                )

            type_i = idx["Experimental System Type"]
            a_i    = idx["Official Symbol Interactor A"]
            b_i    = idx["Official Symbol Interactor B"]
            sys_i  = idx["Experimental System"]
            auth_i = idx["Author"]
            pub_i  = idx["Publication Source"]
            thr_i  = idx["Throughput"]

            for row in reader:
                try:
                    sys_type = row[type_i]
                except IndexError:
                    continue
                if sys_type.lower() != "physical":
                    skipped_genetic += 1
                    continue
                try:
                    a = row[a_i].strip()
                    b = row[b_i].strip()
                except IndexError:
                    continue
                if not a or not b or a == "-" or b == "-":
                    skipped_missing_symbol += 1
                    continue
                rec = (
                    a.upper(),
                    b.upper(),
                    row[sys_i],
                    sys_type,
                    row[auth_i] if auth_i < len(row) else "",
                    _pmid_from_pubsource(row[pub_i] if pub_i < len(row) else ""),
                    row[thr_i] if thr_i < len(row) else "",
                )
                batch.append(rec)
                if len(batch) >= BATCH:
                    cur.executemany(INSERT_SQL, batch)
                    inserted += len(batch)
                    batch.clear()
                    if inserted % 200_000 == 0:
                        print(
                            f"  …{inserted:,} rows "
                            f"({time.perf_counter() - t0:.1f} s)",
                            flush=True,
                        )
            if batch:
                cur.executemany(INSERT_SQL, batch)
                inserted += len(batch)

    print(
        f"  ✓ {inserted:,} physical-interaction rows inserted "
        f"(skipped {skipped_genetic:,} genetic, "
        f"{skipped_missing_symbol:,} missing-symbol) "
        f"in {time.perf_counter() - t0:.1f} s"
    )

    print("Creating indexes")
    t_idx = time.perf_counter()
    cur.execute("CREATE INDEX idx_interactions_a ON interactions (symbol_a)")
    cur.execute("CREATE INDEX idx_interactions_b ON interactions (symbol_b)")
    conn.commit()
    print(f"  ✓ indexes built in {time.perf_counter() - t_idx:.1f} s")

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
    print(f"\nDone. {count:,} rows written to {DB_PATH} ({db_size_mb:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
