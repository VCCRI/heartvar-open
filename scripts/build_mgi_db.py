#!/usr/bin/env python3
"""build_mgi_db.py — build the local MGI mouse-ortholog + phenotype cache.

Downloads MGI's public flat reports from The Jackson Laboratory and the
Mammalian Phenotype (MP) ontology, then writes them into ``data/mgi.db``.
The backend MGI client reads from this database in place of the three
live HTTP lookups (genenames.org HGNC id → Alliance orthologs → Alliance
phenotypes), falling back to the live calls when the DB is absent or the
gene is not found locally.

Two tables are built:

  gene_ortholog(human_symbol PK, hgnc_id, mgi_id, mouse_symbol)
      Built from HOM_MouseHumanSequence.rpt. Human (NCBI Taxon 9606) and
      mouse (10090) rows that share a "DB Class Key" describe the same
      orthology group: the human row carries the HGNC id, the mouse row
      carries the MGI marker id + mouse symbol.

  gene_phenotype(mgi_id, mp_id, phenotype_statement, pubmed_id)
      Built from MGI_GenePheno.rpt (genotype→phenotype annotations with
      MP ids + PubMed ids keyed on the mouse marker MGI id), joined to MP
      term names parsed from mp.obo so phenotype_statement is the human-
      readable MP term label the client returns.

MGI publishes these reports weekly; re-running this script monthly keeps
the cache reasonably fresh for a curation aid. Re-running is safe — both
tables are dropped and rebuilt every run.

Data sources / licence:
  - MGI flat reports (The Jackson Laboratory) — CC-BY 4.0
  - Mammalian Phenotype Ontology (mp.obo) — CC-BY 4.0

Usage::

    python3 build_mgi_db.py                 # use cached downloads if present
    python3 build_mgi_db.py --force-download

The cached downloads live at ``data/HOM_MouseHumanSequence.rpt``,
``data/MGI_GenePheno.rpt`` and ``data/mp.obo``; delete one (or pass
--force-download) to refresh it on the next run.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sqlite3
import ssl
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import _dbbuild

USER_AGENT = "HeartVar/1.0"

HOM_URL = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
GENEPHENO_URL = "https://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt"
MP_OBO_URL = "http://purl.obolibrary.org/obo/mp.obo"
MP_OBO_FALLBACK_URL = (
    "https://raw.githubusercontent.com/obophenotype/"
    "mammalian-phenotype-ontology/master/mp.obo"
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DB_PATH = DATA_DIR / "mgi.db"
HOM_PATH = DATA_DIR / "HOM_MouseHumanSequence.rpt"
GENEPHENO_PATH = DATA_DIR / "MGI_GenePheno.rpt"
MP_OBO_PATH = DATA_DIR / "mp.obo"

HOM_HUMAN_TAXON = "9606"
HOM_MOUSE_TAXON = "10090"
HOM_WANTED_COLS = [
    "DB Class Key",
    "Common Organism Name",
    "NCBI Taxon ID",
    "Symbol",
    "Mouse MGI ID",
    "HGNC ID",
]

DDL_ORTHOLOG = """
CREATE TABLE gene_ortholog (
    human_symbol  TEXT PRIMARY KEY,
    hgnc_id       TEXT,
    mgi_id        TEXT,
    mouse_symbol  TEXT
)
"""

DDL_PHENOTYPE = """
CREATE TABLE gene_phenotype (
    mgi_id               TEXT NOT NULL,
    mp_id                TEXT,
    phenotype_statement  TEXT,
    pubmed_id            TEXT
)
"""


def _download_urllib(url: str, dest: Path) -> None:
    """Stream a URL into ``dest`` using the stdlib (certifi CA bundle),
    sending the project User-Agent so MGI/JAX can attribute the traffic."""
    t0 = time.perf_counter()
    bytes_seen = 0
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req) as r, open(dest, "wb") as f:
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
    system ``curl``, which uses Secure Transport / the macOS keychain and
    so trusts the corporate intercept root automatically."""
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
    """Download ``url`` to ``dest``, falling back to curl on TLS failure
    and to ``fallback_url`` if the primary host is unreachable."""
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


def parse_mp_obo(path: Path) -> dict[str, str]:
    """Walk mp.obo and emit ``MP:NNNNNNN`` → primary ``name`` mapping.

    Same OBO [Term]/id/name parsing as ``build_hpo_labels.py``: skips
    [Typedef] stanzas and entries marked ``is_obsolete: true``; only the
    first ``id:`` line in each [Term] is honoured.
    """
    out: dict[str, str] = {}
    in_term = False
    cur_id: str | None = None
    cur_name: str | None = None
    cur_obsolete = False

    def _flush() -> None:
        if in_term and cur_id and cur_name and not cur_obsolete:
            out[cur_id] = cur_name

    with path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            if line.startswith("[Term]"):
                _flush()
                in_term = True
                cur_id = None
                cur_name = None
                cur_obsolete = False
                continue
            if line.startswith("[") and line.endswith("]"):
                _flush()
                in_term = False
                cur_id = None
                cur_name = None
                cur_obsolete = False
                continue
            if not in_term:
                continue
            if line.startswith("id: ") and cur_id is None:
                cur_id = line[4:].strip()
            elif line.startswith("name: ") and cur_name is None:
                cur_name = line[6:].strip()
            elif line.startswith("is_obsolete: ") and line.endswith("true"):
                cur_obsolete = True
        _flush()
    return {k: v for k, v in out.items() if k.startswith("MP:")}


def build_orthologs(src: Path, conn: sqlite3.Connection) -> int:
    """Group HOM rows by 'DB Class Key', pair the human (9606) row with
    the mouse (10090) row in the same group, and emit one ortholog row per
    human symbol. Returns rows inserted."""
    print(f"Parsing orthologs from {src}", flush=True)
    groups: dict[str, dict[str, dict[str, str]]] = {}
    with src.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if header and header[0].startswith("#"):
            header[0] = header[0].lstrip("#").strip()
        idx = {name: i for i, name in enumerate(header)}
        missing = [c for c in HOM_WANTED_COLS if c not in idx]
        if missing:
            raise RuntimeError(
                f"HOM_MouseHumanSequence header missing expected columns: "
                f"{missing}. MGI may have changed the schema. "
                f"Header was: {header}"
            )
        key_i = idx["DB Class Key"]
        taxon_i = idx["NCBI Taxon ID"]
        symbol_i = idx["Symbol"]
        mgi_i = idx["Mouse MGI ID"]
        hgnc_i = idx["HGNC ID"]

        for row in reader:
            if len(row) <= max(key_i, taxon_i, symbol_i, mgi_i, hgnc_i):
                continue
            key = row[key_i].strip()
            taxon = row[taxon_i].strip()
            symbol = row[symbol_i].strip()
            if not key or not symbol:
                continue
            grp = groups.setdefault(key, {})
            if taxon == HOM_HUMAN_TAXON:
                grp["human"] = {
                    "symbol": symbol,
                    "hgnc_id": row[hgnc_i].strip() or None,
                }
            elif taxon == HOM_MOUSE_TAXON:
                grp["mouse"] = {
                    "symbol": symbol,
                    "mgi_id": row[mgi_i].strip() or None,
                }

    rows: list[tuple] = []
    seen: set[str] = set()
    for grp in groups.values():
        human = grp.get("human")
        mouse = grp.get("mouse")
        if not human or not mouse:
            continue
        human_symbol = human["symbol"]
        if human_symbol in seen:
            continue
        seen.add(human_symbol)
        rows.append((
            human_symbol,
            human.get("hgnc_id"),
            mouse.get("mgi_id"),
            mouse.get("symbol"),
        ))

    conn.executemany(
        "INSERT OR REPLACE INTO gene_ortholog VALUES (?,?,?,?)", rows
    )
    print(f"  ✓ {len(rows):,} human→mouse ortholog rows")
    return len(rows)


GENEPHENO_MP_COL = 4
GENEPHENO_PUBMED_COL = 5
GENEPHENO_MARKER_COL = 6
GENEPHENO_MIN_COLS = 7


def build_phenotypes(
    src: Path, mp_names: dict[str, str], conn: sqlite3.Connection
) -> int:
    """Stream MGI_GenePheno.rpt, resolve each MP id to its term name via
    ``mp_names``, and emit (mgi_id, mp_id, phenotype_statement, pubmed_id)
    rows keyed on the mouse marker MGI id. Returns rows inserted."""
    print(f"Parsing phenotypes from {src}", flush=True)
    BATCH = 50_000
    batch: list[tuple] = []
    inserted = 0
    skipped_no_name = 0
    cur = conn.cursor()
    with src.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < GENEPHENO_MIN_COLS:
                continue
            marker = row[GENEPHENO_MARKER_COL].strip()
            mp_id = row[GENEPHENO_MP_COL].strip()
            pubmed = row[GENEPHENO_PUBMED_COL].strip() or None
            if not marker or not mp_id:
                continue
            statement = mp_names.get(mp_id)
            if not statement:
                skipped_no_name += 1
                continue
            batch.append((marker, mp_id, statement, pubmed))
            if len(batch) >= BATCH:
                cur.executemany(
                    "INSERT INTO gene_phenotype VALUES (?,?,?,?)", batch
                )
                inserted += len(batch)
                batch.clear()
        if batch:
            cur.executemany(
                "INSERT INTO gene_phenotype VALUES (?,?,?,?)", batch
            )
            inserted += len(batch)
    print(
        f"  ✓ {inserted:,} phenotype rows inserted "
        f"(skipped {skipped_no_name:,} with no MP term name)"
    )
    return inserted


def build_db(db_path: Path) -> tuple[int, int]:
    """Rebuild ``mgi.db`` from the cached downloads. Returns
    (ortholog_rows, phenotype_rows)."""
    print(f"Building {db_path}", flush=True)
    tmp_path = _dbbuild.staging_db_path(db_path)
    conn = _dbbuild.connect(tmp_path)
    conn.execute("PRAGMA journal_mode = OFF")
    conn.execute("PRAGMA synchronous = OFF")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute(DDL_ORTHOLOG)
    conn.execute(DDL_PHENOTYPE)

    ortho_count = build_orthologs(HOM_PATH, conn)

    print("Parsing MP ontology term names", flush=True)
    mp_names = parse_mp_obo(MP_OBO_PATH)
    print(f"  ✓ {len(mp_names):,} MP term names")
    if not mp_names:
        raise RuntimeError(
            "No MP: terms parsed from mp.obo — the phenotype join would be "
            "empty. Aborting."
        )

    pheno_count = build_phenotypes(GENEPHENO_PATH, mp_names, conn)

    print("Creating indexes")
    t_idx = time.perf_counter()
    conn.execute(
        "CREATE INDEX idx_gene_phenotype_mgi ON gene_phenotype (mgi_id)"
    )
    conn.commit()
    print(f"  ✓ indexes built in {time.perf_counter() - t_idx:.1f} s")
    conn.close()
    _dbbuild.publish(tmp_path, db_path)
    return ortho_count, pheno_count


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--force-download",
        action="store_true",
        help="redownload the MGI reports + mp.obo even if cached copies exist",
    )
    args = ap.parse_args(argv)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    force = args.force_download
    _ensure(HOM_PATH, HOM_URL, force)
    _ensure(GENEPHENO_PATH, GENEPHENO_URL, force)
    _ensure(MP_OBO_PATH, MP_OBO_URL, force, fallback_url=MP_OBO_FALLBACK_URL)

    ortho_count, pheno_count = build_db(DB_PATH)
    db_size_mb = DB_PATH.stat().st_size / (1024 * 1024)
    print(
        f"\nDone. {ortho_count:,} ortholog rows + {pheno_count:,} phenotype "
        f"rows written to {DB_PATH} ({db_size_mb:.1f} MB)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
