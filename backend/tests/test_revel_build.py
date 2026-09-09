"""Guards on scripts/_build_revel_data.sh — peak disk and output correctness.

Why this file exists: the REVEL step, not the 25 GB cache tarball, is what
actually set the size of the Azure Files share. The original script expanded the
667 MB zip and then wrote THREE full uncompressed copies (~8 GB each) without
deleting the earlier ones, then sorted with scratch on the same volume — ~25 GB
peak to produce a 1.5 GB file. Collapsing it to a single pipeline drops that to
~10 GB, which is the difference between needing ~77 GB of share and ~60 GB.

The end-to-end tests build a tiny synthetic REVEL zip, so they exercise the real
pipeline — column handling, header, '.'-position dropping, sort order, tabix — in
well under a second.
"""
from __future__ import annotations

import gzip
import os
import re
import subprocess
import zipfile
from pathlib import Path
from shutil import which

import pytest

BUILDER = Path(__file__).resolve().parents[2] / "scripts" / "_build_revel_data.sh"
SCRIPT = BUILDER.read_text()
CODE = "\n".join(l for l in SCRIPT.splitlines() if not l.lstrip().startswith("#"))

HEADER = "chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid"
ROWS = [
    "14,23900000,23424081,G,A,R,Q,0.913,ENST00000355349",
    "1,1000,2000,C,T,A,V,0.100,ENST00000000001",
    "14,23800000,23424000,T,C,L,P,0.500,ENST00000355349",
    "2,5000,.,A,G,M,T,0.200,ENST00000000002",
]

_TOOLS = ("bgzip", "tabix", "unzip", "curl")
_HAVE_TOOLS = all(which(t) for t in _TOOLS)


def _fake_zip(path: Path) -> None:
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("revel_with_transcript_ids", "\n".join([HEADER] + ROWS) + "\n")


def _build(tmp_path: Path, assembly: str) -> tuple[subprocess.CompletedProcess, Path]:
    zip_path = tmp_path / "revel.zip"
    _fake_zip(zip_path)
    out_dir = tmp_path / f"revel_{assembly}"
    res = subprocess.run(
        ["bash", str(BUILDER), str(out_dir), assembly],
        capture_output=True, text=True, timeout=180,
        env={"PATH": os.environ["PATH"], "REVEL_ZIP_URL": zip_path.as_uri()},
    )
    return res, out_dir


def _rows(out: Path) -> list[list[str]]:
    lines = gzip.decompress(out.read_bytes()).decode().strip().split("\n")
    return [ln.split("\t") for ln in lines[1:]], lines[0]


def test_no_intermediate_uncompressed_copies():
    """The three temp copies are the bug. None of them may survive."""
    for dead in ('"$tmp/tabbed_revel.tsv"', '"$tmp/new_tabbed_revel.tsv"'):
        assert dead not in SCRIPT, (
            f"{dead} is an intermediate full copy of an ~8 GB file — the pipeline "
            "must stream instead of materialising it."
        )


def test_single_pipeline_feeds_bgzip():
    """The sorted stream must reach bgzip directly, with no temp file between.

    Asserted on the brace group rather than on `sort | bgzip` adjacency: the
    header is printed before the sorted body, so the pipe into bgzip comes off
    the closing brace, not off sort itself.
    """
    assert re.search(r'^\}\s*\|\s*bgzip -c > "\$tmp/out\.tsv\.gz"', SCRIPT,
                     re.MULTILINE), (
        "the brace group holding header + sorted body must pipe straight to bgzip"
    )
    assert not re.search(r'bgzip -c > "\$OUT_TSV"', CODE), (
        "bgzip must not write the FINAL path: a killed build then leaves a "
        "truncated file at the finished name, which every later run skips"
    )
    assert not re.search(r'sort[^\n]*>\s*"?\$tmp', SCRIPT), (
        "sort must not spill to a temp file — that is the 8 GB copy being removed"
    )


def test_sort_scratch_is_redirectable():
    """On Azure Files the sort scratch is the slow, expensive part; a container
    with roomy local disk should be able to use it instead."""
    assert "REVEL_SORT_TMP" in SCRIPT
    assert re.search(r'sort\s+[^\n]*-T\s*"?\$', SCRIPT)


def test_sort_memory_is_bounded():
    """An unbounded sort in a container is an OOM kill, and this deployment has
    already OOM-killed the builder (exit 137, in a recycle loop)."""
    assert re.search(r"sort[^\n]*-S\s", SCRIPT)


def test_scratch_dir_is_created_not_assumed():
    """The trap rm -rf's the scratch dir, so it must be one we made ourselves —
    with a fixed path, REVEL_SORT_TMP=/tmp would have the script delete /tmp."""
    assert re.search(r"mktemp -d[^\n]*SORT|SORT_TMP=\"\$\(mktemp -d", SCRIPT), (
        "the sort scratch must come from mktemp -d, not a fixed path"
    )


def test_idempotent_skip_is_preserved():
    assert re.search(r'-s "\$OUT_TSV"', SCRIPT)


@pytest.mark.skipif(not _HAVE_TOOLS, reason=f"needs {', '.join(_TOOLS)}")
def test_end_to_end_builds_a_sorted_tabixed_file(tmp_path):
    res, out_dir = _build(tmp_path, "GRCh38")
    assert res.returncode == 0, res.stderr

    out = out_dir / "new_tabbed_revel_grch38.tsv.gz"
    assert out.exists(), res.stderr
    assert (out_dir / "new_tabbed_revel_grch38.tsv.gz.tbi").exists()

    body, header = _rows(out)
    assert header.startswith("#chr\t"), f"header must survive, commented, first: {header!r}"
    assert "," not in header, "commas must become tabs"

    assert all(r[2] != "." for r in body), "rows with '.' GRCh38 pos must be dropped"
    assert len(body) == 3, f"expected 3 rows, got {body}"
    chr14 = [int(r[2]) for r in body if r[0] == "14"]
    assert chr14 == sorted(chr14), f"chr14 positions not sorted: {chr14}"
    assert any(r[7] == "0.913" for r in body), body


@pytest.mark.skipif(not _HAVE_TOOLS, reason=f"needs {', '.join(_TOOLS)}")
def test_grch37_indexes_the_hg19_column(tmp_path):
    """Position column differs by assembly — 2 for GRCh37, 3 for GRCh38."""
    res, out_dir = _build(tmp_path, "GRCh37")
    assert res.returncode == 0, res.stderr
    body, _ = _rows(out_dir / "new_tabbed_revel_grch37.tsv.gz")
    assert len(body) == 4, f"expected 4 rows for GRCh37, got {body}"


@pytest.mark.skipif(not _HAVE_TOOLS, reason=f"needs {', '.join(_TOOLS)}")
def test_rerun_is_a_no_op(tmp_path):
    """build_all.sh may re-enter this on any deploy; a rebuild would be ~10 GB of
    pointless work."""
    res, out_dir = _build(tmp_path, "GRCh38")
    assert res.returncode == 0, res.stderr
    out = out_dir / "new_tabbed_revel_grch38.tsv.gz"
    before = out.stat().st_mtime_ns

    zip_path = tmp_path / "revel.zip"
    again = subprocess.run(
        ["bash", str(BUILDER), str(out_dir), "GRCh38"],
        capture_output=True, text=True, timeout=180,
        env={"PATH": os.environ["PATH"], "REVEL_ZIP_URL": zip_path.as_uri()},
    )
    assert again.returncode == 0, again.stderr
    assert "skipping" in again.stdout.lower(), again.stdout
    assert out.stat().st_mtime_ns == before, "the output was rewritten"


@pytest.mark.skipif(not _HAVE_TOOLS, reason=f"needs {', '.join(_TOOLS)}")
def test_scratch_is_cleaned_up(tmp_path):
    """~10 GB of sort scratch left behind on the share would be worse than the
    problem this change fixes."""
    res, out_dir = _build(tmp_path, "GRCh38")
    assert res.returncode == 0, res.stderr
    leftovers = [p.name for p in out_dir.iterdir() if p.name.startswith(".")]
    assert leftovers == [], f"scratch left behind: {leftovers}"


@pytest.mark.skipif(not _HAVE_TOOLS, reason=f"needs {', '.join(_TOOLS)}")
def test_a_truncated_output_is_rebuilt_not_skipped(tmp_path):
    res, out_dir = _build(tmp_path, "GRCh38")
    assert res.returncode == 0, res.stderr
    out = out_dir / "new_tabbed_revel_grch38.tsv.gz"

    good = out.read_bytes()
    out.write_bytes(good[: len(good) // 2])

    zip_path = tmp_path / "revel.zip"
    again = subprocess.run(
        ["bash", str(BUILDER), str(out_dir), "GRCh38"],
        capture_output=True, text=True, timeout=180,
        env={"PATH": os.environ["PATH"], "REVEL_ZIP_URL": zip_path.as_uri()},
    )
    assert again.returncode == 0, again.stdout + again.stderr
    assert "UNUSABLE" in again.stderr, again.stdout + again.stderr
    assert out.read_bytes() == good, "the truncated file was not replaced"


@pytest.mark.skipif(not _HAVE_TOOLS, reason=f"needs {', '.join(_TOOLS)}")
def test_a_file_indexed_on_the_other_assembly_is_rebuilt(tmp_path):
    """THE FILE ACTUALLY ON THE MOUNT on 2026-08-26: `bgzip -t` OK, index
    readable, correct header, and no row at the GRCh38 coordinate of MYH7 R403Q.
    A structurally perfect file indexed on hg19_pos under a grch38 filename. No
    integrity check can see that, so the skip must ask whether the index answers
    in THIS assembly's coordinates."""
    rows = [
        "14,23500000,90000000,G,A,R,Q,0.913,ENST00000355349",
        "14,23600000,90100000,T,C,L,P,0.500,ENST00000355349",
    ]
    zip_path = tmp_path / "revel.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("revel_with_transcript_ids", "\n".join([HEADER] + rows) + "\n")

    out_dir = tmp_path / "revel"
    env = {"PATH": os.environ["PATH"], "REVEL_ZIP_URL": zip_path.as_uri()}
    res = subprocess.run(["bash", str(BUILDER), str(out_dir), "GRCh37"],
                         capture_output=True, text=True, timeout=180, env=env)
    assert res.returncode == 0, res.stderr

    for suffix in ("", ".tbi"):
        (out_dir / f"new_tabbed_revel_grch37.tsv.gz{suffix}").rename(
            out_dir / f"new_tabbed_revel_grch38.tsv.gz{suffix}")

    again = subprocess.run(["bash", str(BUILDER), str(out_dir), "GRCh38"],
                           capture_output=True, text=True, timeout=180, env=env)
    assert again.returncode == 0, again.stdout + again.stderr
    assert "UNUSABLE" in again.stderr, again.stdout + again.stderr
    got = subprocess.run(
        ["tabix", str(out_dir / "new_tabbed_revel_grch38.tsv.gz"), "14:90000000-90000000"],
        capture_output=True, text=True)
    assert got.stdout.strip(), "rebuilt file does not answer on grch38_pos"


@pytest.mark.skipif(not _HAVE_TOOLS, reason=f"needs {', '.join(_TOOLS)}")
def test_force_revel_rebuilds_a_usable_file(tmp_path):
    res, out_dir = _build(tmp_path, "GRCh38")
    assert res.returncode == 0, res.stderr
    out = out_dir / "new_tabbed_revel_grch38.tsv.gz"
    before = out.stat().st_mtime_ns

    zip_path = tmp_path / "revel.zip"
    again = subprocess.run(
        ["bash", str(BUILDER), str(out_dir), "GRCh38"],
        capture_output=True, text=True, timeout=180,
        env={"PATH": os.environ["PATH"], "REVEL_ZIP_URL": zip_path.as_uri(),
             "FORCE_REVEL": "1"},
    )
    assert again.returncode == 0, again.stderr
    assert "skipping" not in again.stdout.lower(), again.stdout
    assert out.stat().st_mtime_ns != before, "FORCE_REVEL did not rebuild"


@pytest.mark.skipif(not _HAVE_TOOLS, reason=f"needs {', '.join(_TOOLS)}")
def test_a_failed_build_leaves_nothing_at_the_final_path(tmp_path):
    """Atomic publish: an interrupted run must leave nothing for the next run to
    mistake for finished output."""
    out_dir = tmp_path / "revel_GRCh38"
    res = subprocess.run(
        ["bash", str(BUILDER), str(out_dir), "GRCh38"],
        capture_output=True, text=True, timeout=180,
        env={"PATH": os.environ["PATH"],
             "REVEL_ZIP_URL": (tmp_path / "does-not-exist.zip").as_uri()},
    )
    assert res.returncode != 0, res.stdout
    assert not (out_dir / "new_tabbed_revel_grch38.tsv.gz").exists()


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
