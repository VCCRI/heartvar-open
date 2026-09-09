#!/usr/bin/env python3
"""build_cdot_db.py — convert the cdot transcript JSONs into a SQLite database.

WHY THIS EXISTS, and it is a hard deployability blocker rather than a tidy-up.
``cdot.hgvs.dataproviders.JSONDataProvider`` parses its transcript files fully
into memory. Measured on 2026-08-28 with the two GRCh38 files (61 MB RefSeq +
43 MB Ensembl 113, 906,754 transcripts):

    baseline                  13 MB
    after imports             41 MB
    after JSONDataProvider  4480 MB      <-- container limit is 4 GiB

So the resolver could not be deployed at all. Resolution itself is cheap (0.2 ms
per variant once loaded); the cost is entirely the in-memory transcript set.

WHAT WAS REJECTED, and why. Restricting the transcript set to the cardiac panel
brought it to 498 MB — but it also drops any variant whose gene is outside the
198-gene panel to the REST fallback. That
trades correctness for memory when the problem is the STORAGE, not the SCOPE. A
gene being off-panel is not a reason to annotate it worse.

WHAT THIS DOES INSTEAD. cdot routes every lookup through one retrieval hook —
``_get_transcript(tx_ac)`` — plus ``_get_transcript_ids_for_gene(gene)``. So the
transcript set can live on disk, indexed, and be fetched a row at a time. Full
coverage, and memory that does not scale with the transcript count.

The row payload is the cdot record VERBATIM. Nothing is reshaped or dropped,
because every field is something cdot's own code may read — the exon list, the
cigar strings, cds bounds, the MANE tag, the contig. Storing a subset would mean
re-deriving cdot's semantics here, which is the class of mistake that has cost
this project a week.

Usage:
    python3 scripts/build_cdot_db.py [--out DATA_DIR] [--keep-json]
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import sqlite3
import sys
import tempfile
import time
from pathlib import Path

import ijson

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = Path(os.environ.get("DATA_DIR") or (PROJECT_ROOT / "data"))

CDOT_DATA_VERSION = os.environ.get("CDOT_DATA_VERSION", "0.2.34")
_FALLBACK_VEP_RELEASE = "113"


def _vep_release() -> str:
    """The VEP release to match, derived from the installed binary.

    ⚠ THE FALLBACK IS NOT BENIGN, which is why it is announced rather than
    swallowed. This value picks WHICH Ensembl transcript file the database is
    built from, and a transcript present in cdot but absent from the VEP cache
    resolves to coordinates that cannot then be annotated. Silently defaulting to
    113 after a failed read is exactly the release-skew this function exists to
    prevent — so every path that does not derive the real value says so on
    stderr, loudly enough to appear in the build log.
    """
    env = os.environ.get("VEP_RELEASE_NUM", "").strip()
    if env:
        return env

    import shutil

    binary = shutil.which("vep")
    if not binary:
        print(f"!!! NOTE  vep is not on PATH, so the VEP release could not be "
              f"derived. Assuming {_FALLBACK_VEP_RELEASE}. Expected on a dev "
              f"machine; in the builder image this means something is wrong.",
              file=sys.stderr)
        return _FALLBACK_VEP_RELEASE

    constants = Path(binary).parent / "modules/Bio/EnsEMBL/VEP/Constants.pm"
    try:
        text = constants.read_text(encoding="utf-8")
    except OSError as exc:
        print(f"!!! NOTE  could not read {constants} ({exc.__class__.__name__}: "
              f"{exc}). Assuming VEP release {_FALLBACK_VEP_RELEASE} — if the "
              f"image is on a different release, the Ensembl transcript file "
              f"will not match the cache.", file=sys.stderr)
        return _FALLBACK_VEP_RELEASE

    for line in text.splitlines():
        if "VEP_VERSION" in line and "=" in line:
            tail = line.split("=", 1)[1].strip().rstrip(";")
            if tail.isdigit():
                return tail

    print(f"!!! NOTE  {constants} contains no parseable VEP_VERSION. Assuming "
          f"{_FALLBACK_VEP_RELEASE}.", file=sys.stderr)
    return _FALLBACK_VEP_RELEASE


SCHEMA = """
CREATE TABLE IF NOT EXISTS transcripts (
    accession TEXT PRIMARY KEY,
    gene      TEXT,
    -- The cdot record verbatim, gzipped. Verbatim because cdot's own code reads
    -- fields this script has no business having an opinion about.
    data      BLOB NOT NULL
);
-- get_tx_for_gene walks every transcript of a gene, so this index is what keeps
-- that from becoming a table scan.
CREATE INDEX IF NOT EXISTS idx_transcripts_gene ON transcripts(gene);

CREATE TABLE IF NOT EXISTS genes (
    -- Keyed EXACTLY as cdot keys them (Entrez id for RefSeq, ENSG for Ensembl),
    -- not by symbol. JSONDataProvider does a bare `self.genes.get(gene)`, so
    -- matching its key space is what keeps _get_gene behaving identically —
    -- including the fact that a SYMBOL lookup misses in both. Re-keying by
    -- symbol here would be a behaviour change dressed up as a fix.
    key   TEXT PRIMARY KEY,
    data  BLOB NOT NULL
);

CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT
);
"""


def _iter_transcripts(path: Path):
    """Stream (kind, key, gene, record) out of one cdot file.

    ⚠⚠ THIS USED json.load, AND THE ASSUMPTION IN ITS DOCSTRING WAS WRONG.
    It said the whole-document parse was "a BUILD-time cost paid once in a
    container with the whole box to itself". The container does NOT have the box
    to itself. Measured on the real Azure Container App Job 2026-08-29:

        staging (local disk): /tmp/cdot_build_1wwaqkoe/cdot_transcripts.db
        build_cdot_transcripts.sh: line 114: 48 Killed  ... build_cdot_db.py

    SIGKILL, six seconds in — the OOM killer. Measured cost of json.load on the
    RefSeq file alone, in this image: 18 MB baseline -> **3601 MB**. So the
    BUILDER needed essentially the same memory as the 4.5 GB in-memory provider
    that this whole script exists to replace. The problem was moved, not solved.

    ijson streams it instead. Measured on the same file, same image:

        json.load    3601 MB   22 s
        ijson         22 MB     5 s      (yajl2_c backend, 518,810 transcripts)

    164x less memory and faster, because nothing is retained between records.

    ⚠ FIDELITY CHECKED, NOT ASSUMED. Each record is stored VERBATIM (json.dumps
    then gzip) and cdot's own c.->g. arithmetic reads it back, so an int quietly
    becoming a float — or a Decimal appearing, which ijson yields by default —
    would move coordinates. 20,000 records compared against the json.load
    result: 0 mismatches, and json.dumps never raised, so no Decimals escaped.
    ``use_float=True`` is set anyway, so a float in some future cdot release
    arrives as a float rather than as a Decimal that would fail to serialise.

    Two passes over the file, one per top-level key. At 5 s each that is cheaper
    than holding either collection in memory.
    """
    opener = gzip.open if path.name.endswith(".gz") else open
    for kind, top_key in (("tx", "transcripts"), ("gene", "genes")):
        with opener(path, "rb") as fh:
            for key, record in ijson.kvitems(fh, top_key, use_float=True):
                gene = (record.get("gene_name") or None) if kind == "tx" else None
                yield kind, key, gene, record


def build(out_dir: Path, files: list[Path], *, keep_json: bool) -> int:
    out_dir.mkdir(parents=True, exist_ok=True)
    final = out_dir / "cdot_transcripts.db"

    #     !!! build failed: OperationalError('database is locked')
    staging = Path(tempfile.mkdtemp(prefix="cdot_build_",
                                    dir=os.environ.get("CDOT_BUILD_TMPDIR") or None))
    tmp = staging / "cdot_transcripts.db"
    print(f"    staging (local disk): {tmp}")

    total = 0
    try:
        conn = sqlite3.connect(str(tmp))
        conn.executescript(SCHEMA)
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")

        for path in files:
            if not path.is_file():
                print(f"!!! missing {path}", file=sys.stderr)
                return 1
            started = time.time()
            tx_rows: list = []
            gene_rows: list = []
            count = genes = 0

            def _flush():
                if tx_rows:
                    conn.executemany(
                        "INSERT OR REPLACE INTO transcripts(accession, gene, data) "
                        "VALUES (?, ?, ?)", tx_rows)
                    tx_rows.clear()
                if gene_rows:
                    conn.executemany(
                        "INSERT OR REPLACE INTO genes(key, data) VALUES (?, ?)",
                        gene_rows)
                    gene_rows.clear()

            for kind, key, gene, record in _iter_transcripts(path):
                blob = gzip.compress(
                    json.dumps(record, separators=(",", ":")).encode("utf-8"), 1)
                if kind == "tx":
                    tx_rows.append((key, gene, blob)); count += 1
                else:
                    gene_rows.append((key, blob)); genes += 1
                if len(tx_rows) + len(gene_rows) >= 5000:
                    _flush()
            _flush()
            conn.commit()
            total += count
            print(f"    {path.name}: {count:,} transcripts, {genes:,} genes "
                  f"({time.time() - started:.0f}s)")

        conn.execute("INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)",
                     ("cdot_data_version", CDOT_DATA_VERSION))
        conn.execute("INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)",
                     ("vep_release", _vep_release()))
        conn.execute("INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)",
                     ("transcript_count", str(total)))
        conn.execute("INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)",
                     ("built_utc", time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                 time.gmtime())))
        conn.commit()
        conn.execute("ANALYZE")
        conn.commit()
        conn.close()
    except Exception as exc:  # pragma: no cover - build-time failure
        shutil.rmtree(staging, ignore_errors=True)
        print(f"!!! build failed: {exc!r}", file=sys.stderr)
        return 1

    if total == 0:
        shutil.rmtree(staging, ignore_errors=True)
        print("!!! no transcripts written", file=sys.stderr)
        return 1

    part = final.with_name(final.name + ".part")
    try:
        shutil.copyfile(tmp, part)
        os.replace(part, final)
    except OSError as exc:
        part.unlink(missing_ok=True)
        shutil.rmtree(staging, ignore_errors=True)
        print(f"!!! could not publish the database to {final}: {exc!r}", file=sys.stderr)
        return 1
    shutil.rmtree(staging, ignore_errors=True)
    size_mb = final.stat().st_size / 1e6
    print(f"  wrote {final} — {total:,} transcripts, {size_mb:.0f} MB")

    if not keep_json:
        for path in files:
            try:
                path.unlink()
                print(f"    removed {path.name} (superseded by the .db)")
            except OSError as exc:
                print(f"    NOTE  could not remove {path.name} "
                      f"({exc.__class__.__name__}: {exc}); {path.stat().st_size / 1e6:.0f} MB "
                      f"left on the mount. Harmless, but it will not self-clean.",
                      file=sys.stderr)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(DATA_DIR / "cdot"),
                    help="directory to write cdot_transcripts.db into")
    ap.add_argument("--keep-json", action="store_true",
                    help="keep the source JSONs (default: remove them, the .db "
                         "supersedes them and they are 104 MB)")
    ap.add_argument("--json-dir", default=None,
                    help="where the downloaded cdot JSONs are (default: --out)")
    args = ap.parse_args()

    out_dir = Path(args.out)
    json_dir = Path(args.json_dir) if args.json_dir else out_dir
    rel = _vep_release()
    files = [
        json_dir / f"cdot-{CDOT_DATA_VERSION}.refseq.GRCh38.json.gz",
        json_dir / f"cdot-{CDOT_DATA_VERSION}.Homo_sapiens_GRCh38_Ensembl_{rel}.gtf.json.gz",
    ]
    print("==> cdot transcript database")
    print(f"    cdot data version : {CDOT_DATA_VERSION}")
    print(f"    Ensembl release   : {rel} (matched to the VEP cache)")
    print(f"    source JSONs      : {json_dir}")
    print(f"    output            : {out_dir / 'cdot_transcripts.db'}")
    return build(out_dir, files, keep_json=args.keep_json)


if __name__ == "__main__":
    raise SystemExit(main())
