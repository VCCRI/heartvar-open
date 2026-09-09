#!/usr/bin/env python3
"""build_gene_id_map.py — compact gene-ID cross-reference that SHIPS.

Writes ``backend/data/gene_id_map.json.gz``: EntrezGene id → Ensembl gene id and
gene symbol → Ensembl gene id, for every HGNC gene that has an Ensembl id.

Why this exists rather than reading hgnc_complete_set.txt at runtime
-------------------------------------------------------------------
``backend/clients/_gene_ids.py`` needs an Ensembl gene id at curation time,
because VEP's RefSeq transcript set (all NM_/NR_/XM_ input) reports ``gene_id``
as an EntrezGene id — see that module. The obvious source, the 17 MB
``hgnc_complete_set.txt``, is BOTH gitignored and listed in ``.dockerignore``:
it is a build-time cache that never reaches production. Reading it at runtime
therefore works on a dev machine and silently fails-open in the deployed
container, which would leave Open Targets reporting "unavailable" for every
RefSeq curation while looking fine locally.

The derived alias map (``hgnc_alias_map.json``, ~1.4 MB) solves the same problem
by being provisioned onto the /app/data Azure Files share and pointed at with
``HGNC_ALIAS_MAP_PATH``. This map deliberately does NOT follow that pattern: it
is tracked and gzipped (~0.5 MB) so it is COPYd into the image with the rest of
``backend/data/``, and a deploy therefore needs nothing placed on the share
first. Gene→Ensembl mappings are near-static, so a file baked into the image
ages far more gracefully than an alias map would.

Values are stored as the ENSG integer suffix (``ENSG00000049540`` → ``49540``)
which, with gzip, keeps the shipped artefact ~0.5 MB.

Input: ``backend/data/hgnc_complete_set.txt`` — run ``build_hgnc_alias_db.py``
first if absent, which downloads it. Cadence: VERSIONED, alongside the alias
map (wired into ``build_all.sh``). Idempotent.

Run:
    .venv/bin/python scripts/build_gene_id_map.py
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "backend" / "data"
TSV_PATH = DATA_DIR / "hgnc_complete_set.txt"
OUT_PATH = DATA_DIR / "gene_id_map.json.gz"

COL_SYMBOL = "symbol"
COL_ENSG = "ensembl_gene_id"
COL_ENTREZ = "entrez_id"


def build() -> dict:
    if not TSV_PATH.exists():
        raise SystemExit(
            f"{TSV_PATH} not found — run scripts/build_hgnc_alias_db.py first "
            "(it downloads the HGNC complete set)."
        )
    symbol_to: dict[str, int] = {}
    entrez_to: dict[str, int] = {}
    with TSV_PATH.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_sym = header.index(COL_SYMBOL)
            i_ensg = header.index(COL_ENSG)
            i_entrez = header.index(COL_ENTREZ)
        except ValueError as e:
            raise SystemExit(
                f"hgnc_complete_set.txt missing an expected column: {e}"
            ) from e
        width = max(i_sym, i_ensg, i_entrez) + 1
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < width:
                continue
            ensg = f[i_ensg].strip()
            if not ensg.startswith("ENSG"):
                continue
            try:
                num = int(ensg[4:])
            except ValueError:
                continue
            symbol = f[i_sym].strip().upper()
            if symbol:
                symbol_to.setdefault(symbol, num)
            entrez = f[i_entrez].strip()
            if entrez:
                entrez_to.setdefault(entrez, num)
    if not symbol_to or not entrez_to:
        raise SystemExit("parsed 0 mappings — refusing to write an empty map")
    return {
        "format": "ENSG integer suffix; ENSG00000049540 is stored as 49540",
        "source": "HGNC complete set (symbol, ensembl_gene_id, entrez_id)",
        "symbol_to_ensg": symbol_to,
        "entrez_to_ensg": entrez_to,
    }


def main() -> None:
    data = build()
    payload = json.dumps(data, separators=(",", ":"), sort_keys=True)
    tmp = OUT_PATH.with_suffix(OUT_PATH.suffix + ".tmp")
    with gzip.GzipFile(filename="", mode="wb", fileobj=tmp.open("wb"), mtime=0) as gz:
        gz.write(payload.encode("utf-8"))
    tmp.replace(OUT_PATH)
    print(
        f"{len(data['symbol_to_ensg']):,} symbols + "
        f"{len(data['entrez_to_ensg']):,} Entrez ids → {OUT_PATH} "
        f"({OUT_PATH.stat().st_size / 1024:.0f} KB gzipped)"
    )


if __name__ == "__main__":
    sys.exit(main())
