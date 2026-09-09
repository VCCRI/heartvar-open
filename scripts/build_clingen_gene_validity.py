#!/usr/bin/env python3
"""build_clingen_gene_validity.py — build the local ClinGen Gene-Disease
Validity cache.

Downloads ClinGen's public **Gene-Disease Validity** curation CSV and writes
one record *per (gene, disease) assertion* for every cardiovascular-panel gene
(filtered to ``cvd_gene_panel.json``) into
``backend/data/clingen_gene_validity.json``.

This is the authoritative, **per-(gene,disease)-pair** mode-of-inheritance +
evidence-tier layer used by the deterministic carrier-status feature
(``backend/acmg/hard_coded.gene_inheritance_modes`` / ``carrier_status``).
ClinGen GV is the only source that records MOI *and* a strength tier at the
gene-disease granularity needed to resolve dual autosomal-dominant/recessive
cardiac genes (e.g. CASQ2, KCNQ1, TTN, JPH2) safely. The runtime rollup merges
this with the live signals HeartVar already holds (GenCC submissions, CHDgene
inheritance codes, PanelApp MOI, ClinGen Gene-Dosage HI==30).

We store ONLY the ClinGen GV fields here (gene / disease / MONDO / MOI /
classification / date / GCEP) — never the GenCC / CHDgene / PanelApp signals —
to avoid duplicating data the backend already loads (those would drift against
their own source files).

ClinGen GV vocabulary (verified against the live download):

    MOI:            AR, AD, XL, SD (semidominant), MT (mitochondrial),
                    UD (undetermined)
    CLASSIFICATION: Definitive, Strong, Moderate, Limited, Disputed,
                    Refuted, "No Known Disease Relationship"

The "adequate evidence" tier used by the carrier-status rollup is
{Definitive, Strong, Moderate}; Limited / Disputed / Refuted / No-Known are
dropped from the determinative path (mirrors how GenCC's Supportive tier is
excluded). This builder keeps ALL rows verbatim — the tier gate is applied at
runtime so the raw evidence stays auditable.

Usage::

    python3 scripts/build_clingen_gene_validity.py                 # cached
    python3 scripts/build_clingen_gene_validity.py --force-download

The cached CSV lives at ``backend/data/clingen_gene_validity_source.csv``;
delete it (or pass --force-download) to refresh. Re-running is safe (the JSON
is rewritten every time).

Data source: ClinGen Gene-Disease Validity
(https://search.clinicalgenome.org/kb/gene-validity), released under CC0.
Please attribute "ClinGen Gene-Disease Validity".
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import ssl
import subprocess
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

URL = "https://search.clinicalgenome.org/kb/gene-validity/download"

USER_AGENT = "HeartVar/1.0"

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "backend" / "data"
PANEL_FILE = DATA_DIR / "cvd_gene_panel.json"
OUT_PATH = DATA_DIR / "clingen_gene_validity.json"
DOWNLOAD_PATH = DATA_DIR / "clingen_gene_validity_source.csv"

GENE_COL = "GENE SYMBOL"
DISEASE_COL = "DISEASE LABEL"
MONDO_COL = "DISEASE ID (MONDO)"
MOI_COL = "MOI"
CLASS_COL = "CLASSIFICATION"
DATE_COL = "CLASSIFICATION DATE"
GCEP_COL = "GCEP"

REQUIRED = [GENE_COL, DISEASE_COL, MOI_COL, CLASS_COL]


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
        raise RuntimeError("curl failed to download the ClinGen GV CSV")
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


def _clean(s: str | None) -> str | None:
    if s is None:
        return None
    s = s.strip()
    return s or None


def _find_header(rows: list[list[str]]) -> tuple[int, dict[str, int], str | None]:
    """Return (header_row_index, {col_name: idx}, file_created_string).

    Scans for the row that contains GENE SYMBOL + CLASSIFICATION + MOI (the
    header), and harvests the "FILE CREATED:" metadata line for provenance.
    """
    file_created = None
    for r in rows[:8]:
        if r and r[0].upper().startswith("FILE CREATED"):
            file_created = r[0].split(":", 1)[1].strip() if ":" in r[0] else r[0]
    for i, r in enumerate(rows):
        upper = [c.strip().upper() for c in r]
        if GENE_COL in upper and CLASS_COL in upper and MOI_COL in upper:
            idx = {c.strip(): j for j, c in enumerate(r)}
            return i, idx, file_created
    raise RuntimeError(
        f"Could not locate the ClinGen GV header row (looked for "
        f"'{GENE_COL}' + '{CLASS_COL}' + '{MOI_COL}'). The schema may have "
        f"changed."
    )


def build(src: Path, out_path: Path, panel: set[str]) -> dict:
    print(f"Building {out_path} from {src}", flush=True)
    with src.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh))
    hdr_i, idx, file_created = _find_header(rows)
    missing = [c for c in REQUIRED if c not in idx]
    if missing:
        raise RuntimeError(
            f"ClinGen GV CSV missing expected columns: {missing}. "
            f"Header seen: {rows[hdr_i]}"
        )
    max_col = max(idx[c] for c in idx)

    genes: dict[str, list[dict]] = {}
    assertions = 0
    for row in rows[hdr_i + 1:]:
        if not row or not row[0].strip():
            continue
        if row[0].startswith("+"):
            continue
        if len(row) <= max_col:
            continue
        sym = (row[idx[GENE_COL]] or "").strip().upper()
        if not sym or sym not in panel:
            continue
        rec = {
            "disease": _clean(row[idx[DISEASE_COL]]),
            "mondo": _clean(row[idx[MONDO_COL]]) if MONDO_COL in idx else None,
            "moi": _clean(row[idx[MOI_COL]]),
            "classification": _clean(row[idx[CLASS_COL]]),
            "date": _clean(row[idx[DATE_COL]]) if DATE_COL in idx else None,
            "gcep": _clean(row[idx[GCEP_COL]]) if GCEP_COL in idx else None,
        }
        genes.setdefault(sym, []).append(rec)
        assertions += 1

    payload = {
        "meta": {
            "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "source": "ClinGen Gene-Disease Validity",
            "source_url": URL,
            "source_file_created": file_created,
            "license": "CC0",
            "panel_file": PANEL_FILE.name,
            "panel_gene_count": len(panel),
            "gene_count": len(genes),
            "assertion_count": assertions,
            "moi_vocab": ["AR", "AD", "XL", "SD", "MT", "UD"],
            "classification_vocab": [
                "Definitive", "Strong", "Moderate", "Limited",
                "Disputed", "Refuted", "No Known Disease Relationship",
            ],
            "note": (
                "ClinGen Gene-Disease Validity assertions ONLY (per gene-disease "
                "pair). Adequate-evidence tier {Definitive, Strong, Moderate} and "
                "the dual-AD/AR resolution are applied at runtime by "
                "gene_inheritance_modes / carrier_status; raw rows kept verbatim "
                "so the evidence stays auditable."
            ),
        },
        "genes": {k: genes[k] for k in sorted(genes)},
    }
    tmp_path = out_path.with_name(out_path.name + ".tmp")
    tmp_path.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(tmp_path, out_path)
    print(
        f"  ✓ {len(genes):,} panel genes, {assertions:,} gene-disease assertions"
    )
    for sym in ("CASQ2", "KCNQ1", "TTN", "TECRL", "GLA", "MYBPC3"):
        recs = genes.get(sym)
        if recs:
            modes = ", ".join(
                f"{r['moi']}:{r['classification']}" for r in recs
            )
            print(f"  spot-check {sym}: {modes}")
        else:
            print(f"  spot-check {sym}: (no ClinGen GV record)")
    return payload


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--force-download",
        action="store_true",
        help="re-download the ClinGen GV CSV even if a cached copy exists",
    )
    args = ap.parse_args(argv)

    if args.force_download and DOWNLOAD_PATH.exists():
        DOWNLOAD_PATH.unlink()
    if not DOWNLOAD_PATH.exists():
        download(URL, DOWNLOAD_PATH)
    else:
        print(f"Using cached {DOWNLOAD_PATH} (pass --force-download to refresh)")

    panel = _load_panel_genes()
    build(DOWNLOAD_PATH, OUT_PATH, panel)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
