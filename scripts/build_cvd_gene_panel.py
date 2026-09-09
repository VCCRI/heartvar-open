#!/usr/bin/env python3
"""build_cvd_gene_panel.py — build the master all-cardiovascular-disease gene list.

This is the MASTER INPUT for every panel-scoped mirror build (gnomAD-freq DB,
SpliceAI slice, Open Targets DB, offline literature). It is the union of:

  1. Every gene on a cardiovascular PanelApp Australia panel at GREEN or AMBER
     confidence (the panel set is exactly what ``panelapp._categorize_panel``
     keeps — cardiomyopathy, arrhythmia / channelopathy, congenital heart
     disease, aortic / vascular / connective-tissue, coronary, inherited lipid
     disorders and pulmonary hypertension).
  2. The ChDGene curated congenital-heart-disease gene list
     (``backend/data/chdgene_genes.json``).

The crawl is reused verbatim from ``build_panelapp_snapshot.py`` so the gene
set can never drift from the snapshot the runtime client reads. RED (confidence
level 1) and "no list" (0) genes are excluded — they are not yet panel-grade
evidence and would bloat the mirrors with low-value targets.

Output: ``backend/data/cvd_gene_panel.json`` — a sorted list of HGNC symbols
plus provenance metadata. Bare gene symbols are facts (not the curated panel
data PanelApp asks us not to redistribute), so this file IS committed; it must
be reproducible and versioned because the downstream builds depend on it.

Usage::

    python3 scripts/build_cvd_gene_panel.py
    python3 scripts/build_cvd_gene_panel.py --min-confidence green   # green only
    python3 scripts/build_cvd_gene_panel.py --max-panels 5           # smoke test
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import httpx

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CHDGENE_PATH = PROJECT_ROOT / "backend" / "data" / "chdgene_genes.json"
OUT_PATH = PROJECT_ROOT / "backend" / "data" / "cvd_gene_panel.json"

sys.path.insert(0, str(PROJECT_ROOT))
from scripts.build_panelapp_snapshot import crawl  # noqa: E402

_CONFIDENCE_FLOOR = {"green": {"3"}, "amber": {"3", "2"}}


def _chdgene_symbols() -> list[str]:
    """HGNC symbols from the committed ChDGene list."""
    try:
        data = json.loads(CHDGENE_PATH.read_text())
    except (OSError, json.JSONDecodeError) as e:
        print(f"[cvd-panel] WARNING: could not read ChDGene list: {e!r}", file=sys.stderr)
        return []
    return [
        (g.get("symbol") or "").strip().upper()
        for g in data.get("genes", [])
        if (g.get("symbol") or "").strip()
    ]


def _panelapp_symbols(snapshot: dict[str, list[dict]], allowed: set[str]) -> list[str]:
    """Symbols with at least one record meeting the confidence floor."""
    return [
        gene
        for gene, records in snapshot.items()
        if any(str(r.get("confidence_level")) in allowed for r in records)
    ]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-confidence", choices=("green", "amber"), default="amber",
                    help="lowest PanelApp confidence to include (default: amber)")
    ap.add_argument("--max-panels", type=int, default=None,
                    help="cap panels fetched (for a quick smoke test)")
    args = ap.parse_args()

    allowed = _CONFIDENCE_FLOOR[args.min_confidence]

    try:
        snapshot = crawl(max_panels=args.max_panels)
    except httpx.HTTPError as e:
        print(f"[cvd-panel] PanelApp crawl failed: {e!r}", file=sys.stderr)
        return 1

    panelapp = _panelapp_symbols(snapshot, allowed)
    chdgene = _chdgene_symbols()
    genes = sorted(set(panelapp) | set(chdgene))
    if not genes:
        print("[cvd-panel] no genes collected — aborting", file=sys.stderr)
        return 1

    out = {
        "meta": {
            "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "min_confidence": args.min_confidence,
            "gene_count": len(genes),
            "sources": {
                "panelapp_aus": {
                    "url": "https://panelapp-aus.org/api/v1/panels/",
                    "genes_contributed": len(set(panelapp)),
                },
                "chdgene": {
                    "url": "https://chdgene.victorchang.edu.au",
                    "genes_contributed": len(set(chdgene)),
                },
            },
        },
        "genes": genes,
    }
    tmp_path = OUT_PATH.with_name(OUT_PATH.name + ".tmp")
    tmp_path.write_text(json.dumps(out, ensure_ascii=False, indent=0) + "\n")
    os.replace(tmp_path, OUT_PATH)
    print(f"[cvd-panel] wrote {OUT_PATH} ({len(genes):,} genes; "
          f"PanelApp {len(set(panelapp)):,} / ChDGene {len(set(chdgene)):,}, "
          f"min_confidence={args.min_confidence})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
