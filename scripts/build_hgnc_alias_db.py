#!/usr/bin/env python3
"""build_hgnc_alias_db.py — build the local HGNC gene-symbol alias map.

Downloads HGNC's "complete set" TSV (the canonical genenames.org gene
nomenclature dump) and distils it into a compact JSON the backend uses to
canonicalise an outdated / alias / previous gene symbol a curator types
into the current HGNC-approved symbol. Without this map, a stale symbol
(e.g. ``CSX`` for ``NKX2-5``) silently produces empty evidence because the
downstream clients query the wrong gene.

The backend client (``backend/clients/hgnc_alias.py``) reads this JSON in
place of a live HGNC lookup, returning a NEUTRAL fail-open result when the
file is absent so a missing optional asset never blocks a curation.

One JSON file is written to ``backend/data/hgnc_alias_map.json`` with this
exact shape::

    {"approved": ["A1BG", "A2M", ...],            # sorted approved symbols (UPPER)
     "aliases": {"ALIASUPPER": "ApprovedSymbol"}}  # alias/prev (UPPER) -> approved

Only rows with ``status == "Approved"`` are kept. The ``alias_symbol`` and
``prev_symbol`` columns are pipe-delimited; each entry is uppercased and
mapped to its approved symbol. When an alias resolves to MORE THAN ONE
approved symbol it is AMBIGUOUS and OMITTED entirely (we never guess) — the
count of omitted-ambiguous aliases is reported in the build log. An approved
symbol is never recorded as its own alias.

HGNC publishes a fresh complete set roughly monthly; re-running this script
keeps the map current. Re-running is safe — the JSON is rebuilt every run.

Data sources / licence:
  - HGNC complete set (EMBL-EBI / genenames.org) — CC0 1.0 (public domain).

Usage::

    python3 build_hgnc_alias_db.py                 # use a cached TSV if present
    python3 build_hgnc_alias_db.py --force-download

The cached download lives at ``backend/data/hgnc_complete_set.txt``; delete
it (or pass --force-download) to refresh it on the next run.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import ssl
import subprocess
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

USER_AGENT = "HeartVar/1.0"

HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/"
    "hgnc_complete_set.txt"
)
HGNC_FALLBACK_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/"
    "hgnc_complete_set.txt"
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "backend" / "data"
TSV_PATH = DATA_DIR / "hgnc_complete_set.txt"
OUT_PATH = DATA_DIR / "hgnc_alias_map.json"

COL_SYMBOL = "symbol"
COL_ALIAS = "alias_symbol"
COL_PREV = "prev_symbol"
COL_STATUS = "status"
STATUS_APPROVED = "Approved"
PIPE = "|"


def _download_urllib(url: str, dest: Path) -> None:
    """Stream a URL into ``dest`` using the stdlib (certifi CA bundle),
    following redirects and sending the project User-Agent so genenames.org /
    EBI can attribute the traffic."""
    t0 = time.perf_counter()
    bytes_seen = 0
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            bytes_seen += len(chunk)
            if bytes_seen % (1 << 23) < (1 << 20):
                print(f"  …{bytes_seen / (1024 * 1024):.1f} MB", flush=True)
    size_mb = dest.stat().st_size / (1024 * 1024)
    print(f"  ✓ {size_mb:.1f} MB in {time.perf_counter() - t0:.1f} s")


def _download_curl(url: str, dest: Path) -> None:
    """Fallback for networks with enterprise SSL inspection: shell out to
    system ``curl``, which uses Secure Transport / the macOS keychain and so
    trusts the corporate intercept root automatically."""
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


def download(url: str, dest: Path, fallback_url: str | None = None) -> None:
    """Download ``url`` to ``dest`` (following redirects), falling back to
    curl on TLS failure and to ``fallback_url`` if the primary host is
    unreachable."""
    print(f"Downloading {url}\n  → {dest}", flush=True)
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        _download_urllib(url, dest)
        return
    except (ssl.SSLError, urllib.error.URLError) as e:
        msg = str(e)
        if "CERTIFICATE_VERIFY_FAILED" in msg or "SSL" in msg:
            if dest.exists():
                dest.unlink()
            _download_curl(url, dest)
            return
        if fallback_url:
            print(f"  …primary host failed ({e}); trying {fallback_url}", flush=True)
            if dest.exists():
                dest.unlink()
            download(fallback_url, dest)
            return
        raise


def _ensure(path: Path, url: str, force: bool, fallback_url: str | None = None) -> None:
    if path.exists() and not force:
        size_mb = path.stat().st_size / (1024 * 1024)
        print(
            f"Using cached {path} ({size_mb:.1f} MB) — "
            "pass --force-download to refresh."
        )
        return
    download(url, path, fallback_url=fallback_url)


def _split_pipes(value: str) -> list[str]:
    """Split a pipe-delimited HGNC cell into trimmed, non-empty tokens.
    HGNC sometimes wraps multi-value cells in double quotes; csv handles
    the quoting, so we only strip whitespace here."""
    if not value:
        return []
    return [tok.strip() for tok in value.split(PIPE) if tok.strip()]


def build_map(src: Path) -> dict:
    """Parse the HGNC complete-set TSV into the compact alias map.

    Two passes: (1) collect every Approved symbol (uppercased) and the
    candidate alias→approved links; (2) drop ambiguous aliases (mapping to
    more than one distinct approved symbol) and aliases that collide with an
    approved symbol. Returns the JSON-ready dict + reports counts.
    """
    print(f"Parsing HGNC complete set from {src}", flush=True)
    approved_upper: set[str] = set()
    candidates: dict[str, set[str]] = defaultdict(set)

    with src.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        missing = [c for c in (COL_SYMBOL, COL_ALIAS, COL_PREV, COL_STATUS)
                   if c not in idx]
        if missing:
            raise RuntimeError(
                f"hgnc_complete_set header missing expected columns: {missing}. "
                f"HGNC may have changed the schema. Header was: {header}"
            )
        sym_i = idx[COL_SYMBOL]
        alias_i = idx[COL_ALIAS]
        prev_i = idx[COL_PREV]
        status_i = idx[COL_STATUS]
        max_i = max(sym_i, alias_i, prev_i, status_i)

        for row in reader:
            if len(row) <= max_i:
                continue
            if row[status_i].strip() != STATUS_APPROVED:
                continue
            symbol = row[sym_i].strip()
            if not symbol:
                continue
            approved_upper.add(symbol.upper())
            for tok in _split_pipes(row[alias_i]) + _split_pipes(row[prev_i]):
                candidates[tok.upper()].add(symbol)

    aliases: dict[str, str] = {}
    ambiguous_omitted = 0
    self_approved_skipped = 0
    for alias_upper, targets in candidates.items():
        if alias_upper in approved_upper:
            self_approved_skipped += 1
            continue
        if len(targets) > 1:
            ambiguous_omitted += 1
            continue
        aliases[alias_upper] = next(iter(targets))

    approved_sorted = sorted(approved_upper)
    print(f"  ✓ {len(approved_sorted):,} approved symbols")
    print(f"  ✓ {len(aliases):,} unambiguous aliases")
    print(f"  • {ambiguous_omitted:,} aliases omitted (ambiguous, >1 approved)")
    print(f"  • {self_approved_skipped:,} alias tokens skipped "
          "(already an approved symbol)")
    return {"approved": approved_sorted, "aliases": aliases}


def write_map(data: dict, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(out_path.name + ".tmp")
    with tmp_path.open("w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, sort_keys=True,
                  separators=(",", ":"))
    os.replace(tmp_path, out_path)
    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(f"  ✓ wrote {out_path} ({size_mb:.2f} MB)")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--force-download",
        action="store_true",
        help="redownload the HGNC complete set even if a cached copy exists",
    )
    args = ap.parse_args(argv)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    _ensure(TSV_PATH, HGNC_URL, args.force_download, fallback_url=HGNC_FALLBACK_URL)

    data = build_map(TSV_PATH)
    write_map(data, OUT_PATH)
    print(
        f"\nDone. {len(data['approved']):,} approved symbols + "
        f"{len(data['aliases']):,} aliases written to {OUT_PATH}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
