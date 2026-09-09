#!/usr/bin/env python3
"""Shared cardiac-gene panel + GRCh38 interval resolver for the local-slice
build scripts (gnomAD frequency, SpliceAI).

The panel is the union of the CHDgene gene list and the VCEP-applicability
gene set already shipped in ``backend/data/`` (~200 genes). Variant-level
slices (gnomAD allele frequencies, SpliceAI delta scores) are restricted to
these genes' GRCh38 spans (padded for near-splice / deep-intronic variants);
anything off-panel falls back to the live API in the client — rare for a
cardiac-curation tool, and the throttle handles it.

GRCh38 intervals are resolved once via the Ensembl REST
``/lookup/symbol/homo_sapiens/{symbol}`` endpoint and cached to a BED at
``data/cardiac_panel.grch38.bed`` (gitignored, build-at-deploy) so the slice
builds don't re-resolve every run. Resolution is paced + sends a descriptive
User-Agent (good citizen); it is a one-off at deploy.

Importable by the sibling ``scripts/build_*.py`` slicers:

    import sys; from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from _cardiac_panel import panel_intervals
"""

from __future__ import annotations

import json
import os
import ssl
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BACKEND_DATA = PROJECT_ROOT / "backend" / "data"
DATA_DIR = PROJECT_ROOT / "data"
DEFAULT_BED = DATA_DIR / "cardiac_panel.grch38.bed"

ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/{sym}?content-type=application/json"
USER_AGENT = "HeartVar/1.0"
DEFAULT_PAD = 5000
REQUEST_PAUSE_SECONDS = 0.2

_MAIN_CONTIGS = {str(c) for c in range(1, 23)} | {"X", "Y", "MT"}


def cardiac_gene_symbols() -> list[str]:
    """Sorted union of CHDgene + VCEP-applicability gene symbols (~200)."""
    syms: set[str] = set()
    chd = json.loads((BACKEND_DATA / "chdgene_genes.json").read_text())
    for g in chd.get("genes", []):
        s = (g.get("symbol") or "").strip().upper()
        if s:
            syms.add(s)
    vcep = json.loads((BACKEND_DATA / "vcep_criteria_spec.json").read_text())
    for s in (vcep.get("genes") or {}):
        s = (s or "").strip().upper()
        if s:
            syms.add(s)
    return sorted(syms)


def _http_get_json(url: str, *, retries: int = 4) -> dict | None:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT,
                                               "Accept": "application/json"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except (ssl.SSLError, urllib.error.URLError) as e:
            msg = str(e)
            if "CERTIFICATE_VERIFY_FAILED" in msg or "SSL" in msg:
                try:
                    out = subprocess.run(
                        ["curl", "--location", "--fail", "--silent", "--show-error",
                         "-H", f"User-Agent: {USER_AGENT}", "-H", "Accept: application/json",
                         url],
                        capture_output=True, text=True, timeout=60,
                    )
                    if out.returncode == 0 and out.stdout.strip():
                        return json.loads(out.stdout)
                except (subprocess.SubprocessError, ValueError):
                    pass
            if attempt < retries - 1:
                wait = 2 ** attempt
                print(f"[panel] Ensembl request failed ({e}); retry {attempt + 1}/{retries - 1} "
                      f"in {wait} s…", file=sys.stderr)
                time.sleep(wait)
                continue
            return None
        except TimeoutError as e:
            if attempt < retries - 1:
                wait = 2 ** attempt
                print(f"[panel] Ensembl read timed out; retry {attempt + 1}/{retries - 1} "
                      f"in {wait} s…", file=sys.stderr)
                time.sleep(wait)
                continue
            return None
        except (urllib.error.HTTPError, ValueError):
            return None
    return None


def _resolve_one(symbol: str) -> tuple[str, int, int] | None:
    """Resolve a gene symbol to (chrom, start, end) on GRCh38, else None."""
    data = _http_get_json(ENSEMBL_LOOKUP.format(sym=symbol))
    if not isinstance(data, dict):
        return None
    chrom = str(data.get("seq_region_name") or "")
    start = data.get("start")
    end = data.get("end")
    if chrom not in _MAIN_CONTIGS or not isinstance(start, int) or not isinstance(end, int):
        return None
    return chrom, start, end


def _write_bed(path: Path, rows: list[tuple[str, int, int, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{c}\t{s}\t{e}\t{g}" for (c, s, e, g) in rows]
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text("\n".join(lines) + "\n")
    os.replace(tmp, path)


def _read_bed(path: Path) -> list[tuple[str, int, int, str]]:
    rows: list[tuple[str, int, int, str]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        rows.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    return rows


def panel_intervals(
    pad: int = DEFAULT_PAD,
    bed_path: Path | None = None,
    force: bool = False,
) -> list[tuple[str, int, int, str]]:
    """Return the cardiac-panel GRCh38 intervals as (chrom, start, end, gene),
    padded by ``pad`` bp each side and clamped at 1.

    Loads the cached BED at ``bed_path`` (default ``data/cardiac_panel.grch38.bed``)
    when present and ``force`` is False; otherwise resolves every panel symbol
    via Ensembl REST (paced, contact UA) and writes the cache. The cached BED
    already includes the padding.
    """
    bed_path = bed_path or DEFAULT_BED
    if bed_path.is_file() and not force:
        return _read_bed(bed_path)

    rows: list[tuple[str, int, int, str]] = []
    symbols = cardiac_gene_symbols()
    missed: list[str] = []
    for sym in symbols:
        iv = _resolve_one(sym)
        if iv is None:
            missed.append(sym)
        else:
            chrom, start, end = iv
            rows.append((chrom, max(1, start - pad), end + pad, sym))
        time.sleep(REQUEST_PAUSE_SECONDS)
    if missed:
        print(f"[panel] WARNING: {len(missed)} symbols unresolved: {missed[:15]}"
              f"{' …' if len(missed) > 15 else ''}", file=sys.stderr)
    if not rows:
        raise RuntimeError("cardiac panel resolved to ZERO intervals — aborting")
    _write_bed(bed_path, rows)
    print(f"[panel] resolved {len(rows)}/{len(symbols)} gene intervals → {bed_path}")
    return rows


if __name__ == "__main__":
    ivs = panel_intervals(force="--force" in sys.argv)
    print(f"{len(ivs)} cardiac-panel intervals; first 3: {ivs[:3]}")
