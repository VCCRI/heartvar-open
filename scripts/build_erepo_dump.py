"""Download the ClinGen Evidence Repository (eRepo) full dump for HeartVar.

Fetches the complete eRepo classifications TSV (~12.6k records across every
VCEP) to ``backend/data/erepo_all.tsv``. The production eRepo client
(``backend/clients/erepo_client.py``) reads this local dump when present and
surfaces the matching VCEP verdict to the curator as a read-only reference;
when the dump is absent it falls back to the per-gene live API (best-effort,
~25 records/gene). Building the dump at deploy gives the curator-reference its
full-fidelity, offline path and removes the per-curation live call.

Usage:
    python scripts/build_erepo_dump.py

The TSV is cached and re-downloaded only when older than 30 days. Refresh it
on the same monthly cadence as the other local caches (see the README).
Path overridable for the client via the ``HEARTVAR_EREPO_TSV`` env var.
"""

from __future__ import annotations

import time
from pathlib import Path

import httpx

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TSV_CACHE = PROJECT_ROOT / "backend" / "data" / "erepo_all.tsv"

EREPO_URL = "https://erepo.clinicalgenome.org/evrepo/api/classifications/all?format=tabbed"
CACHE_MAX_AGE_SECONDS = 30 * 24 * 60 * 60
DOWNLOAD_TIMEOUT_SECONDS = 600.0


def _cache_is_fresh(path: Path) -> bool:
    if not path.exists():
        return False
    age = time.time() - path.stat().st_mtime
    return age < CACHE_MAX_AGE_SECONDS


def _download_tsv(dest: Path) -> None:
    print(f"Downloading eRepo TSV from {EREPO_URL}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    bytes_seen = 0
    started = time.perf_counter()
    last_report = started
    with httpx.Client(timeout=DOWNLOAD_TIMEOUT_SECONDS, follow_redirects=True) as client:
        with client.stream("GET", EREPO_URL) as resp:
            resp.raise_for_status()
            total = int(resp.headers.get("content-length") or 0)
            with tmp.open("wb") as f:
                for chunk in resp.iter_bytes(chunk_size=64 * 1024):
                    f.write(chunk)
                    bytes_seen += len(chunk)
                    now = time.perf_counter()
                    if now - last_report >= 1.0:
                        mb = bytes_seen / (1024 * 1024)
                        if total:
                            pct = 100.0 * bytes_seen / total
                            total_mb = total / (1024 * 1024)
                            print(f"  …{mb:.1f}/{total_mb:.1f} MB ({pct:.0f}%)")
                        else:
                            print(f"  …{mb:.1f} MB")
                        last_report = now
    tmp.replace(dest)
    elapsed = time.perf_counter() - started
    print(f"  done — {bytes_seen / (1024 * 1024):.1f} MB in {elapsed:.0f}s "
          f"→ {dest.relative_to(PROJECT_ROOT)}")


def main() -> int:
    if _cache_is_fresh(TSV_CACHE):
        age_days = (time.time() - TSV_CACHE.stat().st_mtime) / 86400
        print(f"Using cached TSV ({TSV_CACHE.relative_to(PROJECT_ROOT)}, "
              f"{age_days:.1f} days old) — pass nothing to force, or delete it to re-download.")
        return 0
    _download_tsv(TSV_CACHE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
