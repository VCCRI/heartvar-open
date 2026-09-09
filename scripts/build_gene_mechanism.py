#!/usr/bin/env python3
"""build_gene_mechanism.py — build the local ClinGen gene-dosage cache.

Downloads ClinGen's public **Gene-Dosage** curation list and writes one
record per cardiovascular-panel gene (filtered to ``cvd_gene_panel.json``)
into ``backend/data/gene_mechanism.json``. This file is the *only* mechanism
signal HeartVar can't already derive at runtime — the ClinGen
haploinsufficiency / triplosensitivity scores. The runtime classifier
(``backend/acmg/hard_coded.gene_mechanism``) merges this with the live
signals it already holds: VCEP PVS1/PP2 applicability, the CHDgene
dominant-negative inheritance flag, the hand ``_DOMINANT_NEGATIVE_GENES``
list, and gnomAD missense constraint.

We store ONLY the ClinGen dosage fields here — never the VCEP / CHDgene /
hand-list signals — to avoid duplicating data the backend already loads
(those would drift against their own source files).

Haploinsufficiency (and triplosensitivity) score semantics (ClinGen):

    0   No evidence available
    1   Little evidence for dosage pathogenicity
    2   Some evidence for dosage pathogenicity
    3   Sufficient evidence for dosage pathogenicity  -> LoF established
    30  Gene associated with autosomal recessive phenotype (biallelic LoF)
    40  Dosage sensitivity unlikely  -> evidence AGAINST haploinsufficiency

The ClinGen curation list is a small TSV (~1500 dosage-curated genes, a few
hundred KB) and needs no auth. Dosage curations change slowly; a one-off
build is enough — re-run periodically. Re-running is safe (the JSON is
rewritten every time).

Usage::

    python3 scripts/build_gene_mechanism.py                 # cached download
    python3 scripts/build_gene_mechanism.py --force-download

The cached TSV lives at ``backend/data/ClinGen_gene_curation_list_GRCh38.tsv``;
delete it (or pass --force-download) to refresh.

Data source: ClinGen Gene-Dosage (https://clinicalgenome.org/), released
under CC0. Please attribute "ClinGen Dosage Sensitivity Map".
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
from datetime import datetime, timezone
from pathlib import Path

URL = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"

USER_AGENT = "HeartVar/1.0"

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "backend" / "data"
PANEL_FILE = DATA_DIR / "cvd_gene_panel.json"
OUT_PATH = DATA_DIR / "gene_mechanism.json"
DOWNLOAD_PATH = DATA_DIR / "ClinGen_gene_curation_list_GRCh38.tsv"

GENE_COL = "Gene Symbol"
HI_SCORE_COL = "Haploinsufficiency Score"
HI_DESC_COL = "Haploinsufficiency Description"
TS_SCORE_COL = "Triplosensitivity Score"
TS_DESC_COL = "Triplosensitivity Description"
GENE_ID_COL = "Gene ID"

_SCORE_GLOSS = {
    "30": "Gene associated with autosomal recessive phenotype",
    "40": "Dosage sensitivity unlikely",
}


def _download_urllib(url: str, dest: Path) -> None:
    """Stream ``url`` into ``dest`` with the stdlib (certifi CA bundle)."""
    t0 = time.perf_counter()
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req) as r, open(dest, "wb") as f:
        shutil.copyfileobj(r, f)
    size_kb = dest.stat().st_size / 1024
    print(f"  ✓ {size_kb:.0f} KB in {time.perf_counter() - t0:.1f} s")


def _download_curl(url: str, dest: Path) -> None:
    """Fallback: system curl (trusts the macOS keychain / enterprise roots)."""
    curl = shutil.which("curl")
    if not curl:
        raise RuntimeError(
            "Python's TLS chain rejected the server certificate (likely an "
            "enterprise inspection proxy) and `curl` isn't on PATH."
        )
    print(f"  …falling back to {curl} (uses system keychain)", flush=True)
    cmd = [
        curl, "--location", "--fail", "--retry", "3",
        "--silent", "--show-error", "--user-agent", USER_AGENT,
        "--output", str(dest), url,
    ]
    if subprocess.run(cmd).returncode != 0:
        raise RuntimeError("curl failed to download the ClinGen curation list")
    size_kb = dest.stat().st_size / 1024
    print(f"  ✓ {size_kb:.0f} KB")


def download(url: str, dest: Path) -> None:
    print(f"Downloading {url}\n  → {dest}", flush=True)
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        _download_urllib(url, dest)
    except (ssl.SSLError, urllib.error.URLError) as e:
        if "CERTIFICATE_VERIFY_FAILED" in str(e) or "SSL" in str(e):
            if dest.exists():
                dest.unlink()
            _download_curl(url, dest)
        else:
            raise


def _load_panel_genes() -> set[str]:
    payload = json.loads(PANEL_FILE.read_text())
    genes = payload.get("genes") or []
    return {g.strip().upper() for g in genes if g and g.strip()}


def _find_header(src: Path) -> tuple[list[str], int]:
    """Return (header_columns, data_start_byte_offset_unused).

    The ClinGen file is a comment block followed by a header row that itself
    starts with '#'. We scan for the first comment line that contains
    'Gene Symbol' and 'Haploinsufficiency Score' — that is the header — strip
    its leading '#', and split on tabs.
    """
    with src.open("r", encoding="utf-8") as fh:
        for line in fh:
            if GENE_COL in line and HI_SCORE_COL in line:
                cols = line.lstrip("#").rstrip("\n").split("\t")
                return [c.strip() for c in cols], 0
    raise RuntimeError(
        f"Could not locate the ClinGen header row (looked for "
        f"'{GENE_COL}' + '{HI_SCORE_COL}'). The ClinGen schema may have "
        f"changed."
    )


def _clean(s: str | None) -> str | None:
    if s is None:
        return None
    s = s.strip()
    return s or None


def build(src: Path, out_path: Path, panel: set[str]) -> dict:
    print(f"Building {out_path} from {src}", flush=True)
    header, _ = _find_header(src)
    idx = {name: i for i, name in enumerate(header)}
    required = [GENE_COL, HI_SCORE_COL, HI_DESC_COL, TS_SCORE_COL, TS_DESC_COL]
    missing = [c for c in required if c not in idx]
    if missing:
        raise RuntimeError(
            f"ClinGen curation list missing expected columns: {missing}. "
            f"Header seen: {header}"
        )
    max_col = max(idx[c] for c in required + [GENE_ID_COL] if c in idx)

    genes: dict[str, dict] = {}
    rows_seen = 0
    with src.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) <= max_col:
                continue
            sym = (row[idx[GENE_COL]] or "").strip().upper()
            if not sym or sym not in panel:
                continue
            rows_seen += 1
            hi_raw = _clean(row[idx[HI_SCORE_COL]])
            ts_raw = _clean(row[idx[TS_SCORE_COL]])
            gene_id = _clean(row[idx[GENE_ID_COL]]) if GENE_ID_COL in idx else None
            rec: dict = {
                "hi_score": hi_raw,
                "hi_desc": _clean(row[idx[HI_DESC_COL]]) or _SCORE_GLOSS.get(hi_raw or ""),
                "ts_score": ts_raw,
                "ts_desc": _clean(row[idx[TS_DESC_COL]]) or _SCORE_GLOSS.get(ts_raw or ""),
            }
            if gene_id:
                rec["clingen_url"] = (
                    f"https://search.clinicalgenome.org/kb/gene-dosage/HGNC?search={sym}"
                )
            genes[sym] = rec

    payload = {
        "meta": {
            "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "source": "ClinGen Gene-Dosage (Dosage Sensitivity Map)",
            "source_url": URL,
            "license": "CC0",
            "panel_file": PANEL_FILE.name,
            "panel_gene_count": len(panel),
            "gene_count": len(genes),
            "note": (
                "ClinGen haploinsufficiency / triplosensitivity scores ONLY. "
                "Mechanism labels are computed at runtime by merging these with "
                "VCEP PVS1/PP2 applicability, the CHDgene dominant-negative flag, "
                "the hand DN list, and gnomAD constraint."
            ),
        },
        "genes": dict(sorted(genes.items())),
    }
    tmp_path = out_path.with_name(out_path.name + ".tmp")
    tmp_path.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(tmp_path, out_path)
    print(
        f"  ✓ {len(genes):,} panel genes with ClinGen dosage records "
        f"(scanned {rows_seen:,} matching rows)"
    )
    for sym in ("MYH7", "TTN", "GATA4", "NKX2-5", "TBX5"):
        r = genes.get(sym)
        if r:
            print(f"  spot-check {sym}: HI={r['hi_score']} ({r['hi_desc']})")
        else:
            print(f"  spot-check {sym}: (no ClinGen dosage record)")
    return payload


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="re-download the ClinGen curation list even if cached",
    )
    args = parser.parse_args(argv)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    if not DOWNLOAD_PATH.exists() or args.force_download:
        download(URL, DOWNLOAD_PATH)
    else:
        print(
            f"Using cached {DOWNLOAD_PATH.name} — pass --force-download to refresh."
        )
    panel = _load_panel_genes()
    print(f"Loaded {len(panel):,} cardiovascular-panel genes from {PANEL_FILE.name}")
    build(DOWNLOAD_PATH, OUT_PATH, panel)
    print(f"\nDone → {OUT_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
