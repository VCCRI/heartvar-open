"""CHDgene lookup client.

CHDgene (chdgene.victorchang.edu.au) is the Victor Chang lab's manually
curated list of high-confidence CHD genes. Membership is binary: if a gene
is listed, it is established as CHD-associated. We don't tier within the
list — the CHD classification, inheritance modes, and reference count are
contextual info, not a confidence gate.

We cache the gene list as a static JSON file (`backend/data/chdgene_genes.json`).
Regenerate by running:

    python -m backend.clients.chdgene refresh
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from ._http_retry import _with_connect_cap

DATA_FILE = Path(__file__).resolve().parent.parent / "data" / "chdgene_genes.json"
SOURCE_URL = "https://chdgene.victorchang.edu.au/"


_INDEX: dict[str, dict] | None = None
_META: dict | None = None


def _load() -> tuple[dict[str, dict], dict]:
    global _INDEX, _META
    if _INDEX is None:
        if not DATA_FILE.exists():
            _INDEX, _META = {}, {"error": f"missing {DATA_FILE.name} — run 'python -m backend.clients.chdgene refresh'"}
            return _INDEX, _META
        payload = json.loads(DATA_FILE.read_text())
        _META = payload.get("meta", {})
        _INDEX = {g["symbol"].upper(): g for g in payload.get("genes", [])}
    return _INDEX, _META


def lookup(gene: str) -> dict:
    """Return CHDgene record for a gene symbol, or a `listed: False` record.

    Synchronous — the file is local. Wrap in `run_local` from the
    streaming endpoint if you want it to never block the loop (currently it
    won't, since the JSON is small and reads complete in microseconds).
    """
    index, meta = _load()
    if not index:
        return {"ok": False, "error": meta.get("error", "CHDgene data not available")}
    entry = index.get(gene.strip().upper())
    if entry is None:
        return {
            "ok": True,
            "listed": False,
            "gene": gene,
            "data_version": meta.get("fetched_at"),
        }
    return {
        "ok": True,
        "listed": True,
        "gene": entry["symbol"],
        "url": entry.get("url"),
        "chd_classification": entry.get("chd_classification") or [],
        "extra_cardiac_phenotype": entry.get("extra_cardiac_phenotype", False),
        "inheritance": entry.get("inheritance") or [],
        "supporting_references_count": int(entry.get("supporting_references_count") or 0),
        "data_version": meta.get("fetched_at"),
    }


def refresh() -> dict:
    """Fetch the live CHDgene table and rewrite the static JSON file.

    Requires `httpx` and `beautifulsoup4` (already in requirements.txt).
    """
    import httpx
    from bs4 import BeautifulSoup

    r = httpx.get(SOURCE_URL, timeout=_with_connect_cap(30.0), follow_redirects=True)
    r.raise_for_status()
    soup = BeautifulSoup(r.text, "lxml")
    table = soup.find("table")
    if table is None:
        raise RuntimeError("CHDgene homepage: no <table> found")

    genes: list[dict] = []
    for row in table.find_all("tr"):
        tds = row.find_all("td")
        if len(tds) != 6:
            continue
        link = tds[0].find("a")
        if not link:
            continue
        symbol = link.get_text(strip=True)
        href = link.get("href") or ""
        chd = [b.get_text(strip=True) for b in tds[1].find_all("span", class_="badge")]
        extra_cardiac = tds[2].get("data-text", "0") == "1"
        inh = [b.get_text(strip=True) for b in tds[3].find_all("span", class_="badge")]
        ranking_stars = int(tds[4].get("data-text", "0") or 0)
        refs_text = tds[5].get_text(strip=True)
        refs_count = int(refs_text) if refs_text.isdigit() else 0
        genes.append({
            "symbol": symbol,
            "url": f"https://chdgene.victorchang.edu.au{href}" if href.startswith("/") else href,
            "chd_classification": chd,
            "extra_cardiac_phenotype": extra_cardiac,
            "inheritance": inh,
            "ranking_stars": ranking_stars,
            "supporting_references_count": refs_count,
        })

    payload = {
        "meta": {
            "source_url": SOURCE_URL,
            "fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "gene_count": len(genes),
            "columns": {
                "ranking_stars": "uniform 5/5 — all genes in CHDgene meet the high-confidence tier",
                "supporting_references_count": "0–6 — count of curated supporting publications; higher = stronger evidence base",
            },
        },
        "genes": sorted(genes, key=lambda g: g["symbol"]),
    }

    DATA_FILE.parent.mkdir(parents=True, exist_ok=True)
    DATA_FILE.write_text(json.dumps(payload, indent=2) + "\n")

    global _INDEX, _META
    _INDEX = None
    _META = None

    return {"ok": True, "gene_count": len(genes), "path": str(DATA_FILE)}


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "refresh":
        res = refresh()
        print(f"Wrote {res['gene_count']} genes to {res['path']}")
    else:
        print(__doc__)
