"""Guards on scripts/build_cdot_db.py — chiefly WHERE it writes.

THE BUG THIS EXISTS FOR. The builder ran `tempfile.mkstemp(dir=out_dir)` and
opened SQLite on it, i.e. directly on /app/data — Azure Files over SMB.
Measured on the real mount 2026-08-29, 9 seconds into a build that takes ~50 s
locally:

    !!! build failed: OperationalError('database is locked')

Azure Files does not support the byte-range locks SQLite takes when writing.
This is the one class of defect a local end-to-end run cannot find, because
local filesystems lock fine — so the property has to be asserted structurally:
the database is BUILT on local disk and PUBLISHED to the destination.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import sqlite3
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "build_cdot_db.py"


def _load():
    spec = importlib.util.spec_from_file_location("build_cdot_db", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _tiny_cdot(path: Path) -> None:
    blob = {
        "transcripts": {
            "NM_000257.4": {"gene_name": "MYH7", "chrom": "14", "genome_builds": {}},
        },
        "genes": {"4625": {"gene_symbol": "MYH7"}},
    }
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        json.dump(blob, fh)


def test_the_database_is_not_built_on_the_destination_filesystem(tmp_path, monkeypatch):
    """⚠ THE ONE THAT MATTERS. Every sqlite3.connect during the build must
    target a path OUTSIDE out_dir, because out_dir is an SMB mount in
    production and SQLite cannot take write locks there."""
    mod = _load()
    out_dir = tmp_path / "mount" / "cdot"
    out_dir.mkdir(parents=True)
    src = out_dir / "cdot-tiny.refseq.GRCh38.json.gz"
    _tiny_cdot(src)

    connected: list[str] = []
    real_connect = sqlite3.connect

    def _spy(target, *a, **k):
        connected.append(str(target))
        return real_connect(target, *a, **k)

    monkeypatch.setattr(mod.sqlite3, "connect", _spy)

    rc = mod.build(out_dir, [src], keep_json=True)
    assert rc == 0, "the tiny build should succeed"
    assert connected, "no SQLite connection was made — the test proves nothing"
    for target in connected:
        assert not str(target).startswith(str(out_dir)), (
            f"SQLite was opened for WRITING at {target}, which is under the "
            "destination directory. In production that is Azure Files over SMB "
            "and the build dies with OperationalError('database is locked')."
        )


def test_the_finished_database_lands_in_out_dir_and_opens(tmp_path):
    """Publishing must still put a working database where the app expects it —
    a copy across a filesystem boundary plus a rename, since os.replace cannot
    cross devices."""
    mod = _load()
    out_dir = tmp_path / "cdot"
    out_dir.mkdir()
    src = out_dir / "cdot-tiny.refseq.GRCh38.json.gz"
    _tiny_cdot(src)

    assert mod.build(out_dir, [src], keep_json=True) == 0
    final = out_dir / "cdot_transcripts.db"
    assert final.is_file(), "the published database is missing"
    conn = sqlite3.connect(f"file:{final}?mode=ro", uri=True)
    try:
        rows = conn.execute("SELECT COUNT(*) FROM transcripts").fetchone()[0]
    finally:
        conn.close()
    assert rows >= 1


def test_no_staging_leftovers_in_the_destination(tmp_path):
    """A `.part` or `.cdot_*` file left on the share would be mistaken for the
    real thing by anything that checks for presence rather than usability."""
    mod = _load()
    out_dir = tmp_path / "cdot"
    out_dir.mkdir()
    src = out_dir / "cdot-tiny.refseq.GRCh38.json.gz"
    _tiny_cdot(src)
    assert mod.build(out_dir, [src], keep_json=True) == 0
    leftovers = [p.name for p in out_dir.iterdir()
                 if p.name.endswith(".part") or p.name.startswith(".cdot_")]
    assert not leftovers, f"staging files left on the destination: {leftovers}"


def test_the_json_is_streamed_not_loaded_whole():
    """⚠ THE BUILDER MUST NOT NEED THE MEMORY IT EXISTS TO SAVE.

    json.load on the cdot RefSeq file peaks at 3601 MB (measured in the builder
    image: 18 MB baseline -> 3601 MB, 518,810 transcripts). The Azure Container
    App Job OOM-killed the build six seconds in:

        build_cdot_transcripts.sh: line 114: 48 Killed ... build_cdot_db.py

    So the BUILDER needed essentially the same memory as the 4.5 GB in-memory
    JSONDataProvider that clients/cdot_sqlite exists to replace — the problem
    had been moved, not solved. ijson streams the same file in 22 MB and 5 s.

    Asserted structurally because the failure needs a memory-capped container to
    reproduce, and every dev machine has enough RAM to hide it."""
    text = SCRIPT.read_text()
    assert "ijson" in text, "the builder must stream the cdot JSON"
    assert "json.load(" not in text, (
        "regression: json.load holds the whole document (3601 MB measured) and "
        "the Azure job has less than that"
    )
    assert "use_float=True" in text, (
        "ijson yields Decimal by default; records are re-serialised verbatim "
        "into the .db, and a Decimal would fail to serialise"
    )


def test_ijson_is_pinned_in_requirements():
    """It is a real runtime dependency of the data build, and the builder image
    installs requirements.txt — an unpinned or missing entry is an ImportError
    in a monthly job nobody is watching."""
    req = (ROOT / "requirements.txt").read_text()
    assert "ijson==" in req


def test_streamed_records_are_identical_to_a_whole_document_parse(tmp_path):
    """FIDELITY, on the round trip the database actually performs.

    Records are stored verbatim (json.dumps then gzip) and cdot's own c.->g.
    arithmetic reads them back, so an int becoming a float would move
    coordinates. Against the real files 20,000 records compared clean; this
    pins the property on a fixture so it cannot regress silently."""
    mod = _load()
    src = tmp_path / "cdot-tiny.refseq.GRCh38.json.gz"
    blob = {
        "transcripts": {
            "NM_1.1": {"gene_name": "MYH7", "start": 23429278, "end": 23429279,
                       "exons": [[1, 2, 3], [4, 5, 6]], "flag": True,
                       "score": 0.886, "nothing": None},
        },
        "genes": {"4625": {"gene_symbol": "MYH7", "map_location": "14q11.2"}},
    }
    with gzip.open(src, "wt", encoding="utf-8") as fh:
        json.dump(blob, fh)

    streamed = {(k, key): rec for k, key, _g, rec in mod._iter_transcripts(src)}
    assert json.dumps(streamed[("tx", "NM_1.1")], sort_keys=True) == \
        json.dumps(blob["transcripts"]["NM_1.1"], sort_keys=True)
    assert json.dumps(streamed[("gene", "4625")], sort_keys=True) == \
        json.dumps(blob["genes"]["4625"], sort_keys=True)
    rec = streamed[("tx", "NM_1.1")]
    assert isinstance(rec["start"], int) and not isinstance(rec["start"], bool)
    assert isinstance(rec["exons"][0][0], int)
    assert isinstance(rec["score"], float)


def test_the_staging_location_is_overridable():
    """A container whose local disk cannot hold ~600 MB needs a way out, and it
    must be documented in the script rather than discovered."""
    text = SCRIPT.read_text()
    assert "CDOT_BUILD_TMPDIR" in text
    assert "mkdtemp" in text, "the build must stage in its own directory"
    assert "mkstemp(dir=str(out_dir))" not in text, (
        "regression: the build is creating its SQLite file on the destination "
        "filesystem again"
    )


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
