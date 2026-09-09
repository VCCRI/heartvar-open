"""Unit tests for backend.clients.spliceai — local-first SpliceAI lookup.

Covers, with NO network:

  (a) BUILD ROUNDTRIP: write a tiny synthetic *source* SpliceAI VCF, run the
      real scripts/build_spliceai_db.py slicer against it (with a fixture
      panel BED), then read the produced slice back through the client and
      assert the success shape + max_delta. Proves the parse/bgzip/tabix/
      read roundtrip without the 27 GB download.

  (b) LOCAL hit: a hand-built bgzipped+tabix VCF with a high DS_DG returns
      ok=True, correct max_delta, and the documented scores_per_transcript
      shape — local path wins (live fallback is sabotaged to assert it).

  (c) IN-PANEL MISS: an in-panel position with no SpliceAI record returns
      {ok:False, lookup_failed:False, reason:"No SpliceAI score at this
      position"} (matches the live empty-scores branch).

  (d) LIVE FALLBACK: slice absent → live Broad API (mocked with
      httpx.MockTransport) returns the original success shape; off-panel and
      GRCh37 also route to live; indel routes to the not_applicable skip.

Runnable with pytest (``python -m pytest backend/tests/test_spliceai.py``) or
directly (``python -m backend.tests.test_spliceai``).
"""
from __future__ import annotations

import asyncio
import importlib.util
import os
import tempfile
from pathlib import Path

import httpx
import pysam

import backend.clients.spliceai as spliceai

_REAL_ASYNC_CLIENT = httpx.AsyncClient

_INFO_HIGH = "SpliceAI=T|MYH7|0.02|0.01|0.97|0.03|-7|12|2|-25"

_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _write_bgzip_tabix_vcf(tmpdir: Path, records: list[str], contig_prefix: str = "") -> Path:
    """Write a plain VCF with ``records`` (each a tab-joined VCF data line),
    bgzip + tabix-index it, return the .vcf.gz path."""
    plain = tmpdir / "slice.vcf"
    plain.write_text(_VCF_HEADER + "".join(
        r if r.endswith("\n") else r + "\n" for r in records))
    gz = tmpdir / "slice.vcf.gz"
    pysam.tabix_compress(str(plain), str(gz), force=True)
    pysam.tabix_index(str(gz), preset="vcf", force=True)
    return gz


def _write_panel_bed(tmpdir: Path, rows: list[tuple[str, int, int, str]]) -> Path:
    """Write a fixture cardiac-panel BED (chrom, start, end, gene)."""
    bed = tmpdir / "panel.bed"
    bed.write_text("\n".join(f"{c}\t{s}\t{e}\t{g}" for (c, s, e, g) in rows) + "\n")
    return bed


def _set_env(slice_path: str | None, bed_path: str | None) -> dict[str, str | None]:
    """Point the client's env overrides at the fixtures; return the prior
    values so the caller can restore them. Also reset the panel cache."""
    prev = {
        "SPLICEAI_DB_PATH": os.environ.get("SPLICEAI_DB_PATH"),
        "SPLICEAI_PANEL_BED_PATH": os.environ.get("SPLICEAI_PANEL_BED_PATH"),
    }
    if slice_path is None:
        os.environ.pop("SPLICEAI_DB_PATH", None)
    else:
        os.environ["SPLICEAI_DB_PATH"] = slice_path
    if bed_path is None:
        os.environ.pop("SPLICEAI_PANEL_BED_PATH", None)
    else:
        os.environ["SPLICEAI_PANEL_BED_PATH"] = bed_path
    spliceai._reset_panel_cache()
    return prev


def _restore_env(prev: dict[str, str | None]) -> None:
    for k, v in prev.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v
    spliceai._reset_panel_cache()


def _no_live_client(*_a, **_k):
    raise AssertionError(
        "the Broad SpliceAI API must NOT be called — this deployment is "
        "local-only unless SPLICEAI_ALLOW_LIVE is set"
    )


def _no_live(*_a, **_k):
    raise AssertionError("live fallback must NOT run on a local-answerable lookup")


def _load_build_module():
    """Import scripts/build_spliceai_db.py by path (scripts/ isn't a package)."""
    build_path = (
        Path(__file__).resolve().parent.parent.parent / "scripts" / "build_spliceai_db.py"
    )
    spec = importlib.util.spec_from_file_location("build_spliceai_db", build_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_build_roundtrip_then_client_read():
    """Run the real build_slice() over a synthetic source VCF, then read the
    produced slice through fetch_spliceai. Proves parse + bgzip + tabix +
    client-read end to end with no download."""
    tmp = Path(tempfile.mkdtemp(prefix="spliceai_build_"))
    src_records = [
        f"14\t23424000\t.\tC\tT\t.\t.\t{_INFO_HIGH}\n",
        f"14\t99999999\t.\tA\tG\t.\t.\t{_INFO_HIGH}\n",
    ]
    source_gz = _write_bgzip_tabix_vcf(tmp, src_records)

    bed = _write_panel_bed(tmp, [("14", 23423000, 23425000, "MYH7")])

    build = _load_build_module()
    out_vcf = tmp / "spliceai_cardiac.masked.grch38.vcf.gz"
    orig_panel = build.panel_intervals
    build.panel_intervals = lambda pad=0, **_k: [("14", 23423000, 23425000, "MYH7")]
    try:
        n = build.build_slice(str(source_gz), out_vcf, pad=0)
    finally:
        build.panel_intervals = orig_panel

    assert n == 1, f"expected exactly 1 in-panel record sliced, got {n}"
    assert out_vcf.is_file() and Path(str(out_vcf) + ".tbi").is_file()

    prev = _set_env(str(out_vcf), str(bed))
    orig_live = spliceai._fetch_spliceai_live
    spliceai._fetch_spliceai_live = _no_live
    try:
        res = asyncio.run(spliceai.fetch_spliceai("14-23424000-C-T", "GRCh38"))
    finally:
        spliceai._fetch_spliceai_live = orig_live
        _restore_env(prev)

    assert res["ok"] is True and res["lookup_failed"] is False, res
    assert res["variant_id"] == "14-23424000-C-T", res
    assert res["assembly"] == "GRCh38", res
    assert abs(res["max_delta"] - 0.97) < 1e-9, res
    assert res["model_message"] == "High confidence splice impact", res
    assert len(res["scores_per_transcript"]) == 1, res
    entry = res["scores_per_transcript"][0]
    assert entry["gene"] == "MYH7", entry
    assert entry["transcript_id"] is None, entry
    assert abs(entry["DS_DG"] - 0.97) < 1e-9, entry
    assert abs(entry["max_delta"] - 0.97) < 1e-9, entry
    assert entry["DP_DG"] == 2 and entry["DP_AG"] == -7, entry
    assert set(entry) == {
        "gene", "transcript_id",
        "DS_AG", "DS_AL", "DS_DG", "DS_DL",
        "DP_AG", "DP_AL", "DP_DG", "DP_DL", "max_delta",
    }, entry
    print("test_build_roundtrip_then_client_read: PASS")


def test_local_hit_returns_correct_shape():
    tmp = Path(tempfile.mkdtemp(prefix="spliceai_local_"))
    gz = _write_bgzip_tabix_vcf(tmp, [f"14\t23424000\t.\tC\tT\t.\t.\t{_INFO_HIGH}\n"])
    bed = _write_panel_bed(tmp, [("14", 23423000, 23425000, "MYH7")])
    prev = _set_env(str(gz), str(bed))
    orig_live = spliceai._fetch_spliceai_live
    spliceai._fetch_spliceai_live = _no_live
    try:
        res = asyncio.run(spliceai.fetch_spliceai("14-23424000-C-T", "GRCh38"))
    finally:
        spliceai._fetch_spliceai_live = orig_live
        _restore_env(prev)
    assert res["ok"] is True, res
    assert abs(res["max_delta"] - 0.97) < 1e-9, res
    assert set(res) == {
        "ok", "lookup_failed", "variant_id", "assembly",
        "max_delta", "model_message", "scores_per_transcript",
    }, res
    print("test_local_hit_returns_correct_shape: PASS")


def test_local_hit_lowercase_assembly_and_chr_prefix_bed():
    """Case-insensitive assembly + a chr-prefixed BED still resolves a bare-
    contig variant id."""
    tmp = Path(tempfile.mkdtemp(prefix="spliceai_chr_"))
    gz = _write_bgzip_tabix_vcf(tmp, [f"14\t23424000\t.\tC\tT\t.\t.\t{_INFO_HIGH}\n"])
    bed = _write_panel_bed(tmp, [("chr14", 23423000, 23425000, "MYH7")])
    prev = _set_env(str(gz), str(bed))
    try:
        res = asyncio.run(spliceai.fetch_spliceai("14-23424000-C-T", "grch38"))
    finally:
        _restore_env(prev)
    assert res["ok"] is True and abs(res["max_delta"] - 0.97) < 1e-9, res
    print("test_local_hit_lowercase_assembly_and_chr_prefix_bed: PASS")


def test_in_panel_miss_returns_no_score():
    """An in-panel position with no SpliceAI record at the queried locus →
    {ok:False, lookup_failed:False, reason:"No SpliceAI score at this
    position"} (NOT a live fallback)."""
    tmp = Path(tempfile.mkdtemp(prefix="spliceai_miss_"))
    gz = _write_bgzip_tabix_vcf(tmp, [f"14\t23424000\t.\tC\tT\t.\t.\t{_INFO_HIGH}\n"])
    bed = _write_panel_bed(tmp, [("14", 23423000, 23425000, "MYH7")])
    prev = _set_env(str(gz), str(bed))
    orig_live = spliceai._fetch_spliceai_live
    spliceai._fetch_spliceai_live = _no_live
    try:
        res = asyncio.run(spliceai.fetch_spliceai("14-23424500-A-G", "GRCh38"))
    finally:
        spliceai._fetch_spliceai_live = orig_live
        _restore_env(prev)
    assert res["ok"] is False and res["lookup_failed"] is False, res
    assert res["reason"] == "No SpliceAI score at this position", res
    print("test_in_panel_miss_returns_no_score: PASS")


def test_empty_slice_reports_lookup_failed_and_never_calls_live():
    """REGRESSION, 2026-09-02. An EMPTY slice must FAIL LOUDLY, offline.

    Production returned {ok:False, lookup_failed:False, reason:"No SpliceAI
    score at this position"} for EVERY position probed — MYH7 and MYBPC3,
    including a canonical splice site scoring 0.890 in the repo's own slice —
    while the identical code answered ok=True locally. Cause: a slice with no
    usable contigs produced an empty ``query_contigs``, and that branch returned
    ``_no_score``. With lookup_failed=False the UI painted the row green and
    BP7 / PP3 / BP4 / PVS1 read a total data outage as a real negative.

    The fix must NOT be a live fallback: this deployment is offline by design,
    so a broken slice is reported as lookup_failed=True and the Broad API is
    never called."""
    tmp = Path(tempfile.mkdtemp(prefix="spliceai_empty_"))
    gz = _write_bgzip_tabix_vcf(tmp, [])
    bed = _write_panel_bed(tmp, [("14", 23423000, 23425000, "MYH7")])
    prev = _set_env(str(gz), str(bed))
    orig_live = spliceai._fetch_spliceai_live
    spliceai._fetch_spliceai_live = _no_live
    try:
        res = asyncio.run(spliceai.fetch_spliceai("14-23424000-C-T", "GRCh38"))
    finally:
        spliceai._fetch_spliceai_live = orig_live
        _restore_env(prev)
    assert res["ok"] is False, res
    assert res["lookup_failed"] is True, "a broken slice must be a LOUD failure"
    assert res.get("reason") != "No SpliceAI score at this position", res
    assert "NO CONTIGS" in res.get("error", ""), res
    print("test_empty_slice_reports_lookup_failed_and_never_calls_live: PASS")


def test_slice_missing_the_queried_contig_reports_lookup_failed():
    """Same class: the slice is readable and non-empty but carries no record for
    the queried contig, while the panel BED says the position IS in-panel. That
    is a broken/partial slice, not an answer — and still no live call."""
    tmp = Path(tempfile.mkdtemp(prefix="spliceai_contig_"))
    gz = _write_bgzip_tabix_vcf(tmp, [f"11\t47331882\t.\tC\tT\t.\t.\t{_INFO_HIGH}\n"])
    bed = _write_panel_bed(tmp, [("14", 23423000, 23425000, "MYH7"),
                                 ("11", 47331000, 47332000, "MYBPC3")])
    prev = _set_env(str(gz), str(bed))
    orig_live = spliceai._fetch_spliceai_live
    spliceai._fetch_spliceai_live = _no_live
    try:
        res = asyncio.run(spliceai.fetch_spliceai("14-23424000-C-T", "GRCh38"))
    finally:
        spliceai._fetch_spliceai_live = orig_live
        _restore_env(prev)
    assert res["ok"] is False and res["lookup_failed"] is True, res
    assert "does not carry contig" in res.get("error", ""), res
    print("test_slice_missing_the_queried_contig_reports_lookup_failed: PASS")


def test_in_panel_position_match_allele_mismatch_is_no_score():
    """Record exists at the position but its ALT (and SpliceAI ALLELE) differs
    from the variant → in-panel miss, not a wrong-score success."""
    tmp = Path(tempfile.mkdtemp(prefix="spliceai_allele_"))
    gz = _write_bgzip_tabix_vcf(tmp, [f"14\t23424000\t.\tC\tT\t.\t.\t{_INFO_HIGH}\n"])
    bed = _write_panel_bed(tmp, [("14", 23423000, 23425000, "MYH7")])
    prev = _set_env(str(gz), str(bed))
    orig_live = spliceai._fetch_spliceai_live
    spliceai._fetch_spliceai_live = _no_live
    try:
        res = asyncio.run(spliceai.fetch_spliceai("14-23424000-C-A", "GRCh38"))
    finally:
        spliceai._fetch_spliceai_live = orig_live
        _restore_env(prev)
    assert res["ok"] is False and res["lookup_failed"] is False, res
    assert res["reason"] == "No SpliceAI score at this position", res
    print("test_in_panel_position_match_allele_mismatch_is_no_score: PASS")


def _make_live_client(scores: list[dict] | None = None, error: str | None = None):
    """AsyncClient wired to a MockTransport emulating the Broad SpliceAI
    Lookup endpoint."""
    def handler(request: httpx.Request) -> httpx.Response:
        if error is not None:
            return httpx.Response(200, json={"variant": "x", "error": error})
        return httpx.Response(200, json={"variant": "x", "scores": scores or []})
    return _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler), follow_redirects=True)


def test_slice_absent_is_a_loud_local_failure_not_a_live_call():
    """Slice file absent → lookup_failed=True, and the Broad API is NOT called.

    DECISION: no live SpliceAI lookups; the slice must work
    offline. A missing slice is the DEPLOYMENT GAP that caused the 2026-09-02
    production outage, so it has to be loud rather than quietly degraded."""
    prev = _set_env("/nonexistent/heartvar-test/spliceai.vcf.gz", None)
    orig_client = spliceai.httpx.AsyncClient
    spliceai.httpx.AsyncClient = _no_live_client
    try:
        res = asyncio.run(spliceai.fetch_spliceai("14-23424000-C-T", "GRCh38"))
    finally:
        spliceai.httpx.AsyncClient = orig_client
        _restore_env(prev)
    assert res["ok"] is False and res["lookup_failed"] is True, res
    assert "slice not found" in res["reason"], res
    print("test_slice_absent_is_a_loud_local_failure_not_a_live_call: PASS")


def test_slice_absent_still_uses_live_when_explicitly_allowed():
    """The escape hatch stays reachable and tested: SPLICEAI_ALLOW_LIVE=1
    restores the Broad fallback, so the decision is an env var, not a patch."""
    prev = _set_env("/nonexistent/heartvar-test/spliceai.vcf.gz", None)
    prev_live = os.environ.get("SPLICEAI_ALLOW_LIVE")
    os.environ["SPLICEAI_ALLOW_LIVE"] = "1"
    orig_client = spliceai.httpx.AsyncClient
    spliceai.httpx.AsyncClient = lambda *a, **k: _make_live_client(scores=[{
        "gene_name": "MYH7", "transcript_id": "ENST00000355349",
        "DS_AG": 0.0, "DS_AL": 0.0, "DS_DG": 0.91, "DS_DL": 0.0,
        "DP_AG": -7, "DP_AL": 12, "DP_DG": 2, "DP_DL": -25,
    }])
    try:
        res = asyncio.run(spliceai.fetch_spliceai("14-23424000-C-T", "GRCh38"))
    finally:
        spliceai.httpx.AsyncClient = orig_client
        if prev_live is None:
            os.environ.pop("SPLICEAI_ALLOW_LIVE", None)
        else:
            os.environ["SPLICEAI_ALLOW_LIVE"] = prev_live
        _restore_env(prev)
    assert res["ok"] is True and abs(res["max_delta"] - 0.91) < 1e-9, res
    print("test_slice_absent_still_uses_live_when_explicitly_allowed: PASS")


def test_off_panel_is_a_scope_answer_not_a_failure():
    """Off-panel is a TRUE STATEMENT ABOUT SCOPE, not a fault: the slice is
    panel-scoped by design. lookup_failed stays False so the row stays green
    and the criteria read "no splice score", not "broken input"."""
    tmp = Path(tempfile.mkdtemp(prefix="spliceai_offpanel_"))
    gz = _write_bgzip_tabix_vcf(tmp, [f"14\t23424000\t.\tC\tT\t.\t.\t{_INFO_HIGH}\n"])
    bed = _write_panel_bed(tmp, [("14", 23423000, 23425000, "MYH7")])
    prev = _set_env(str(gz), str(bed))
    orig_client = spliceai.httpx.AsyncClient
    spliceai.httpx.AsyncClient = _no_live_client
    try:
        res = asyncio.run(spliceai.fetch_spliceai("14-99999999-C-T", "GRCh38"))
    finally:
        spliceai.httpx.AsyncClient = orig_client
        _restore_env(prev)
    assert res["ok"] is False, res
    assert res["lookup_failed"] is False, "off-panel is scope, not a fault"
    assert "outside the cardiac panel" in res["reason"], res
    print("test_off_panel_is_a_scope_answer_not_a_failure: PASS")


def test_panel_bed_absent_is_a_loud_local_failure():
    """BED absent → membership undeterminable. We must not mislabel an in-panel
    miss as off-panel, so this is a fault, reported locally, with no live call."""
    tmp = Path(tempfile.mkdtemp(prefix="spliceai_nobed_"))
    gz = _write_bgzip_tabix_vcf(tmp, [f"14\t23424000\t.\tC\tT\t.\t.\t{_INFO_HIGH}\n"])
    prev = _set_env(str(gz), "/nonexistent/heartvar-test/panel.bed")
    orig_client = spliceai.httpx.AsyncClient
    spliceai.httpx.AsyncClient = _no_live_client
    try:
        res = asyncio.run(spliceai.fetch_spliceai("14-23424000-C-T", "GRCh38"))
    finally:
        spliceai.httpx.AsyncClient = orig_client
        _restore_env(prev)
    assert res["ok"] is False and res["lookup_failed"] is True, res
    assert "BED" in res["reason"] and "unreadable" in res["reason"], res
    print("test_panel_bed_absent_is_a_loud_local_failure: PASS")


def test_grch37_is_a_loud_local_failure():
    """GRCh37 → the slice is GRCh38-only. Report it; do not call the Broad."""
    tmp = Path(tempfile.mkdtemp(prefix="spliceai_37_"))
    gz = _write_bgzip_tabix_vcf(tmp, [f"14\t23424000\t.\tC\tT\t.\t.\t{_INFO_HIGH}\n"])
    bed = _write_panel_bed(tmp, [("14", 23423000, 23425000, "MYH7")])
    prev = _set_env(str(gz), str(bed))
    orig_client = spliceai.httpx.AsyncClient
    spliceai.httpx.AsyncClient = _no_live_client
    try:
        res = asyncio.run(spliceai.fetch_spliceai("14-23424000-C-T", "GRCh37"))
    finally:
        spliceai.httpx.AsyncClient = orig_client
        _restore_env(prev)
    assert res["ok"] is False and res["lookup_failed"] is True, res
    assert "GRCh38-only" in res["reason"], res
    print("test_grch37_is_a_loud_local_failure: PASS")


def test_indel_skip_preserved():
    """Indel id → skipped/not_applicable shape (no slice read, no network)."""
    tmp = Path(tempfile.mkdtemp(prefix="spliceai_indel_"))
    gz = _write_bgzip_tabix_vcf(tmp, [f"14\t23424000\t.\tC\tT\t.\t.\t{_INFO_HIGH}\n"])
    bed = _write_panel_bed(tmp, [("14", 23423000, 23425000, "MYH7")])
    prev = _set_env(str(gz), str(bed))
    try:
        res = asyncio.run(spliceai.fetch_spliceai("8-60845257-A--", "GRCh38"))
    finally:
        _restore_env(prev)
    assert res.get("skipped") is True and res.get("not_applicable") is True, res
    assert "variant_id" in res and "reason" in res, res
    print("test_indel_skip_preserved: PASS")


if __name__ == "__main__":
    test_build_roundtrip_then_client_read()
    test_local_hit_returns_correct_shape()
    test_local_hit_lowercase_assembly_and_chr_prefix_bed()
    test_in_panel_miss_returns_no_score()
    test_in_panel_position_match_allele_mismatch_is_no_score()
    test_live_fallback_when_slice_absent()
    test_off_panel_falls_through_to_live()
    test_panel_bed_absent_fails_safe_to_live()
    test_grch37_routes_to_live()
    test_indel_skip_preserved()
    print("ALL SPLICEAI TESTS PASS")
