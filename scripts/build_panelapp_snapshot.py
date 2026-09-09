#!/usr/bin/env python3
"""build_panelapp_snapshot.py — build the local PanelApp Australia snapshot.

Crawls the public PanelApp Australia REST API once, keeps the cardiovascular
panels (same name-keyword filter the backend client uses), and inverts them
into a ``{GENE_UPPER: [gene-on-panel record, …]}`` map written to
``data/panelapp_aus_snapshot.json``. The backend ``fetch_panelapp`` reads this
snapshot in place of a live ``panelapp-aus.org`` call per curation; the live
API remains the fallback when the snapshot is absent.

Each stored record matches the live ``/genes/`` item shape the client
post-processes::

    {"panel": {"id", "name", "version"},
     "confidence_level", "mode_of_inheritance", "phenotypes", "entity_status"}

LICENSE NOTE: PanelApp panels are "publicly available to browse, download and
query" but the project found NO explicit redistribution grant for the curated
panel data (the GEL codebase is Apache-2.0; that doesn't cover the data). So
this snapshot is BUILD-AT-DEPLOY and gitignored — do not commit it. Local use
is clearly within "download and query". Cite Martin et al., Nat Genet 2019.

Crawl etiquette: a descriptive User-Agent is sent and requests are
paced; PanelApp publishes no hard rate limit but this keeps us a good citizen.
Re-run weekly (UCSC mirrors PanelApp weekly).

Usage::

    python3 scripts/build_panelapp_snapshot.py
    python3 scripts/build_panelapp_snapshot.py --max-panels 5   # smoke test
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import httpx

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
OUT_PATH = DATA_DIR / "panelapp_aus_snapshot.json"

API = "https://panelapp-aus.org/api/v1"
USER_AGENT = "HeartVar/1.0"
REQUEST_PAUSE_SECONDS = 0.2

sys.path.insert(0, str(PROJECT_ROOT))
from backend.clients.panelapp import _categorize_panel  # noqa: E402


def _get_json(client: httpx.Client, url: str, params: dict | None = None) -> dict:
    r = client.get(url, params=params, timeout=60.0)
    r.raise_for_status()
    return r.json()


def crawl(max_panels: int | None = None) -> dict[str, list[dict]]:
    """Return {GENE_UPPER: [record, …]} for every cardiovascular panel."""
    snapshot: dict[str, list[dict]] = {}
    n_panels_seen = 0
    n_panels_kept = 0
    with httpx.Client(headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
                      follow_redirects=True) as client:
        url: str | None = f"{API}/panels/"
        cardiac_panel_ids: list[int] = []
        while url:
            payload = _get_json(client, url)
            for panel in payload.get("results", []) or []:
                n_panels_seen += 1
                if _categorize_panel(panel.get("name", "") or "") is not None:
                    cardiac_panel_ids.append(panel.get("id"))
            url = payload.get("next")
            time.sleep(REQUEST_PAUSE_SECONDS)
        print(f"[panelapp] {n_panels_seen} panels listed; "
              f"{len(cardiac_panel_ids)} cardiovascular")

        for pid in cardiac_panel_ids:
            if max_panels is not None and n_panels_kept >= max_panels:
                break
            detail = _get_json(client, f"{API}/panels/{pid}/")
            panel_ctx = {
                "id": detail.get("id"),
                "name": detail.get("name", "") or "",
                "version": detail.get("version"),
            }
            for g in detail.get("genes", []) or []:
                symbol = (g.get("entity_name") or "").strip().upper()
                if not symbol:
                    continue
                snapshot.setdefault(symbol, []).append({
                    "panel": panel_ctx,
                    "confidence_level": g.get("confidence_level"),
                    "mode_of_inheritance": g.get("mode_of_inheritance"),
                    "phenotypes": g.get("phenotypes") or [],
                    "entity_status": g.get("entity_status"),
                })
            n_panels_kept += 1
            time.sleep(REQUEST_PAUSE_SECONDS)
        print(f"[panelapp] fetched {n_panels_kept} panels; "
              f"{len(snapshot)} genes indexed")
    return snapshot


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-panels", type=int, default=None,
                    help="cap panels fetched (for a quick smoke test)")
    args = ap.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    try:
        genes = crawl(max_panels=args.max_panels)
    except httpx.HTTPError as e:
        print(f"[panelapp] crawl failed: {e!r}", file=sys.stderr)
        return 1
    if not genes:
        print("[panelapp] no cardiovascular genes collected — aborting",
              file=sys.stderr)
        return 1
    out = {
        "meta": {
            "source": f"{API}/panels/",
            "fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "gene_count": len(genes),
        },
        "genes": genes,
    }
    tmp_path = OUT_PATH.with_name(OUT_PATH.name + ".tmp")
    tmp_path.write_text(json.dumps(out, ensure_ascii=False) + "\n")
    os.replace(tmp_path, OUT_PATH)
    print(f"[panelapp] wrote {OUT_PATH} ({len(genes):,} genes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
