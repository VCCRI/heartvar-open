#!/usr/bin/env python3
"""build_opentargets_db.py — offline Open Targets gene-disease snapshot.

Open Targets is otherwise the one about-sources entry that hits the network on
EVERY curation (live GraphQL). This builder snapshots, per cardiac-panel gene,
the top-50 associated diseases + datatype scores into ``data/opentargets.db`` so
backend/clients/opentargets_client.py can read them offline. All the per-disease
HPO->EFO matching and top-N fallback run client-side over these raw rows, so the
snapshot only needs the raw ``target`` object (approvedSymbol + the
associatedDiseases.rows) per Ensembl gene id — exactly what the live ``data``
field returns.

Panel genes: backend/data/cvd_gene_panel.json. Symbol -> Ensembl gene id:
backend/data/hgnc_complete_set.txt (``ensembl_gene_id`` column). Live Open
Targets GraphQL is queried at BUILD time only (paced), never at curation.

Cadence: upstream is VERSIONED (Open Targets platform releases ~quarterly), but
build_all.sh classifies this `monthly` so the scheduled job picks up a new release
without anyone having to notice it shipped (decision 2026-08-13). A rebuild costs
one paced GraphQL sweep of the ~677-gene panel.
Idempotent (drop + rebuild). Writes ``data/opentargets.db`` (honours the
``OPENTARGETS_DB_PATH`` override the reader also uses). Built on local disk and
atomically published to the mount via _dbbuild (SMB-safe).
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _dbbuild  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
BACKEND_DATA = REPO_ROOT / "backend" / "data"
DATA_DIR = REPO_ROOT / "data"
DB_PATH = Path(os.environ.get("OPENTARGETS_DB_PATH") or (DATA_DIR / "opentargets.db"))
PANEL_JSON = BACKEND_DATA / "cvd_gene_panel.json"
HGNC_TSV = BACKEND_DATA / "hgnc_complete_set.txt"

GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
QUERY = """
query TopDiseases($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    associatedDiseases(page: {size: 50, index: 0}) {
      rows { disease { id name } score datatypeScores { id score } }
    }
  }
}
""".strip()

_PACE = 0.25
_UA = "heartvar-databuild/1.0 (variant-curation; heartvar@victorchang.edu.au)"
_MAX_RETRY = 3


def _load_panel_symbols() -> list[str]:
    data = json.loads(PANEL_JSON.read_text())
    genes = data.get("genes") or []
    return sorted({str(g).strip().upper() for g in genes if str(g).strip()})


def _ensure_hgnc_tsv() -> None:
    """Make sure the HGNC complete set is on disk, fetching it if it is not.

    ``hgnc_complete_set.txt`` is both gitignored AND dockerignored, so it does not
    exist in a fresh container — it appears only because
    ``build_hgnc_alias_db.py`` downloads it. This builder used to assume that had
    already happened and died on a bare ``FileNotFoundError`` when it had not,
    which is exactly how every Open Targets build failed up to 2026-08-13
    (build_all.sh ran this step *before* hgnc_alias).

    The ordering in build_all.sh is fixed, but an implicit cross-script
    dependency on a gitignored intermediate is fragile — ``--skip hgnc_alias``,
    or a future reordering, would silently reintroduce the same failure. So
    fetch it here rather than depend on run order. Reuses the downloader from
    build_hgnc_alias_db (retries, curl fallback on TLS failure, mirror URL)
    instead of duplicating it.
    """
    if HGNC_TSV.exists() and HGNC_TSV.stat().st_size > 0:
        return
    print(
        f"[opentargets] {HGNC_TSV.name} absent — fetching it "
        "(normally build_hgnc_alias_db.py provides this)",
        flush=True,
    )
    try:
        from build_hgnc_alias_db import HGNC_FALLBACK_URL, HGNC_URL, download
    except ImportError as e:  # pragma: no cover — same directory, on sys.path
        raise SystemExit(
            f"cannot import the HGNC downloader ({e}). Run "
            "scripts/build_hgnc_alias_db.py first — it downloads "
            f"{HGNC_TSV.name}."
        ) from e
    try:
        download(HGNC_URL, HGNC_TSV, fallback_url=HGNC_FALLBACK_URL)
    except Exception as e:  # noqa: BLE001 — turn any failure into an actionable message
        raise SystemExit(
            f"could not download {HGNC_TSV.name} ({e}). Run "
            "scripts/build_hgnc_alias_db.py first, then re-run this builder."
        ) from e


def _symbol_to_ensembl() -> dict[str, str]:
    """Map uppercased gene symbol -> Ensembl gene id from the HGNC complete set."""
    _ensure_hgnc_tsv()
    m: dict[str, str] = {}
    with HGNC_TSV.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_sym = header.index("symbol")
            i_ensg = header.index("ensembl_gene_id")
        except ValueError:
            raise SystemExit(
                "hgnc_complete_set.txt missing symbol / ensembl_gene_id columns"
            )
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= max(i_sym, i_ensg):
                continue
            sym = cols[i_sym].strip().upper()
            ensg = cols[i_ensg].strip()
            if sym and ensg.startswith("ENSG"):
                m.setdefault(sym, ensg)
    return m


def _query_ot(ensembl_id: str) -> dict | None:
    """POST the GraphQL query and return the ``target`` object (or None)."""
    body = json.dumps({"query": QUERY, "variables": {"ensemblId": ensembl_id}}).encode()
    last: Exception | None = None
    for attempt in range(_MAX_RETRY):
        try:
            req = urllib.request.Request(
                GRAPHQL,
                data=body,
                headers={"Content-Type": "application/json", "User-Agent": _UA},
            )
            with urllib.request.urlopen(req, timeout=30) as r:
                payload = json.loads(r.read().decode())
            if payload.get("errors"):
                print(
                    f"    [ot] GraphQL errors for {ensembl_id}: {payload['errors']}",
                    file=sys.stderr,
                )
            return (payload.get("data") or {}).get("target")
        except (urllib.error.URLError, TimeoutError, ValueError) as e:
            last = e
            time.sleep(_PACE * (attempt + 2))
    raise RuntimeError(f"query failed after {_MAX_RETRY} attempts: {last!r}")


def main() -> int:
    symbols = _load_panel_symbols()
    sym2ensg = _symbol_to_ensembl()
    ensg_by_sym = {s: sym2ensg[s] for s in symbols if s in sym2ensg}
    missing = [s for s in symbols if s not in sym2ensg]
    print(
        f"[opentargets] panel={len(symbols)} genes; resolved ENSG for "
        f"{len(ensg_by_sym)}; unresolved={len(missing)}"
    )
    if missing:
        shown = ", ".join(missing[:15])
        print(f"[opentargets] no ENSG (skipped): {shown}{' …' if len(missing) > 15 else ''}")

    tmp = _dbbuild.staging_db_path(DB_PATH)
    con = _dbbuild.connect(tmp)
    con.execute("PRAGMA journal_mode = OFF")
    con.execute("PRAGMA synchronous = OFF")
    con.execute("DROP TABLE IF EXISTS opentargets")
    con.execute(
        "CREATE TABLE opentargets ("
        "ensembl_id TEXT PRIMARY KEY, symbol TEXT, payload TEXT)"
    )

    ok = fail = 0
    total = len(ensg_by_sym)
    for i, (sym, ensg) in enumerate(sorted(ensg_by_sym.items()), 1):
        try:
            target = _query_ot(ensg)
        except Exception as e:  # noqa: BLE001 — log + continue; one gene shouldn't abort
            print(f"    [ot] {sym} ({ensg}) FAILED: {e!r}", file=sys.stderr)
            fail += 1
            time.sleep(_PACE)
            continue
        if target:
            con.execute(
                "INSERT OR REPLACE INTO opentargets (ensembl_id, symbol, payload)"
                " VALUES (?, ?, ?)",
                (ensg, sym, json.dumps(target, separators=(",", ":"))),
            )
            ok += 1
        if i % 50 == 0:
            con.commit()
            print(f"    [ot] {i}/{total} …")
        time.sleep(_PACE)

    con.commit()
    con.execute("CREATE INDEX IF NOT EXISTS ix_ot_symbol ON opentargets(symbol)")
    con.close()
    _dbbuild.publish(tmp, DB_PATH)
    print(f"[opentargets] wrote {DB_PATH} — stored={ok}, failed={fail}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
