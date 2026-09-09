"""Unit tests for backend.clients.gnomad (local-first frequency + constraint).

Covers every path of the refactored gnomAD client without touching the
network or needing the multi-TB gnomAD VCFs:

  (A) BUILD ROUNDTRIP: write a tiny synthetic bgzipped + tabix-indexed VCF with
      a couple of variants, run the freq build script's slice/parse over it
      (vcf_record_to_payload_fields + merge_payload), and assert the produced
      payload matches the GraphQL variant shape — proving the parsing+read
      roundtrip without the real download.

  (B) LOCAL FREQ HIT: build a tiny gnomad_freq.db (the production schema, JSON
      payload column), point GNOMAD_FREQ_DB_PATH at it, and assert fetch_gnomad
      returns variant_found=True + the EXACT stored blob — with the live
      fallback hard-blocked so a local hit can't silently fall through.

  (C) IN-PANEL MISS → ABSENT: a variant inside a (monkeypatched) panel interval
      that's NOT in an existing freq DB returns variant_found=False / ok=True
      with NO network call (drives PM2 "absent").

  (D) OFF-PANEL / DB-ABSENT → LIVE: a variant outside the panel (DB present) and
      a variant when the DB is absent both fall back to the mocked live GraphQL
      query, preserving the original variant_found/variant shape and the indel
      rsID resolved_by path.

  (E) LOCAL CONSTRAINT: the gene block is served from a temp constraint DB
      (constraint-only path AND alongside a freq hit); live fallback when the
      constraint DB is absent.

Runnable with pytest (``python -m pytest backend/tests/test_gnomad.py``) or
directly (``python backend/tests/test_gnomad.py``).
"""
from __future__ import annotations

import asyncio
import json
import os
import sqlite3
import sys
import tempfile
from pathlib import Path

import httpx
import pysam

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
_SCRIPTS = _REPO_ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

import backend.clients.gnomad as gnomad  # noqa: E402
import build_gnomad_freq_db as freq_build  # noqa: E402

_REAL_ASYNC_CLIENT = httpx.AsyncClient


_PAYLOAD = {
    "variantId": "14-23412711-G-T",
    "rsid": "rs12147570",
    "exome": None,
    "genome": {
        "ac": 19395,
        "an": 152196,
        "af": 0.127434,
        "ac_hom": 1362,
        "ac_hemi": 0,
        "filters": [],
        "populations": [
            {"id": "nfe", "ac": 10000, "an": 68042, "ac_hom": 800, "ac_hemi": 0},
            {"id": "afr", "ac": 2000, "an": 40000, "ac_hom": 100, "ac_hemi": 0},
        ],
        "faf95": {"popmax": 0.15555222, "popmax_population": "nfe"},
    },
}

_FREQ_DDL = """
CREATE TABLE gnomad_freq (
    variant_id  TEXT PRIMARY KEY,
    payload     TEXT
)
"""

_CONSTRAINT_DDL = """
CREATE TABLE gnomad_constraint (
    gene_symbol     TEXT,
    gene_id         TEXT,
    pLI             REAL,
    oe_lof          REAL,
    oe_lof_upper    REAL,
    oe_mis          REAL,
    oe_mis_upper    REAL,
    mis_z           REAL,
    syn_z           REAL
)
"""


def _make_freq_db(rows: list[tuple[str, dict]]) -> Path:
    """Write a tiny gnomad_freq.db with the production schema. ``rows`` is
    [(variant_id, payload_dict), …]."""
    tmpdir = tempfile.mkdtemp(prefix="gnomad_freq_test_")
    db_path = Path(tmpdir) / "gnomad_freq.db"
    conn = sqlite3.connect(db_path)
    conn.execute(_FREQ_DDL)
    conn.executemany(
        "INSERT INTO gnomad_freq VALUES (?, ?)",
        [(vid, json.dumps(p, separators=(",", ":"))) for vid, p in rows],
    )
    conn.commit()
    conn.close()
    return db_path


def _make_constraint_db(rows: list[tuple]) -> Path:
    """Write a tiny gnomad_constraint.db with the production schema. Each row is
    (gene_symbol, gene_id, pLI, oe_lof, oe_lof_upper, oe_mis, oe_mis_upper,
    mis_z, syn_z)."""
    tmpdir = tempfile.mkdtemp(prefix="gnomad_constraint_test_")
    db_path = Path(tmpdir) / "gnomad_constraint.db"
    conn = sqlite3.connect(db_path)
    conn.execute(_CONSTRAINT_DDL)
    conn.executemany("INSERT INTO gnomad_constraint VALUES (?,?,?,?,?,?,?,?,?)", rows)
    conn.execute(
        "CREATE INDEX idx_gnomad_constraint_gene_upper "
        "ON gnomad_constraint (UPPER(gene_symbol))"
    )
    conn.commit()
    conn.close()
    return db_path


class _EnvDBs:
    """Context manager that points the env-resolved DB paths at temp files and
    forces every panel-membership / live call to be explicit, then restores
    the originals + the panel cache on exit."""

    def __init__(self, freq_db: Path | None, constraint_db: Path | None):
        self.freq_db = freq_db
        self.constraint_db = constraint_db
        self._saved_env: dict[str, str | None] = {}
        self._saved_panel = None

    def __enter__(self):
        for key, val in (
            ("GNOMAD_FREQ_DB_PATH", self.freq_db),
            ("GNOMAD_CONSTRAINT_DB_PATH", self.constraint_db),
        ):
            self._saved_env[key] = os.environ.get(key)
            if val is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = str(val)
        self._saved_panel = gnomad._PANEL_BY_CHROM
        gnomad._PANEL_BY_CHROM = None
        return self

    def __exit__(self, *exc):
        for key, prev in self._saved_env.items():
            if prev is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = prev
        gnomad._PANEL_BY_CHROM = self._saved_panel
        return False


def _set_panel(intervals_by_chrom: dict[str, list[tuple[int, int]]] | None):
    """Inject a panel-membership map directly (bypasses Ensembl resolution).
    None means 'panel unavailable' (forces live fallback)."""
    gnomad._PANEL_BY_CHROM = False if intervals_by_chrom is None else intervals_by_chrom


def _block_live():
    """Replace every live entry point with a tripwire so a test that should
    stay local can't silently hit the network. Returns a restore callable."""
    saved = (
        gnomad._fetch_gnomad_live,
        gnomad._fetch_constraint_only_live,
        gnomad._fetch_constraint_block_live,
    )

    async def _boom(*a, **k):  # noqa: ANN001, ANN002
        raise AssertionError("live fallback must NOT run on this path")

    gnomad._fetch_gnomad_live = _boom
    gnomad._fetch_constraint_only_live = _boom
    gnomad._fetch_constraint_block_live = _boom

    def restore():
        (gnomad._fetch_gnomad_live,
         gnomad._fetch_constraint_only_live,
         gnomad._fetch_constraint_block_live) = saved

    return restore


def _write_synthetic_vcf(path: Path, records: list[str]) -> str:
    """Write a minimal VCF (header + ``records`` data lines), bgzip + tabix it
    via pysam. Returns the path to the .gz file (the .tbi sits beside it)."""
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr14>\n"
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"\">\n"
        "##INFO=<ID=nhomalt,Number=A,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=AC_XY,Number=A,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=fafmax_faf95_max,Number=1,Type=Float,Description=\"\">\n"
        "##INFO=<ID=fafmax_faf95_max_gen_anc,Number=1,Type=String,Description=\"\">\n"
        "##INFO=<ID=AC_nfe,Number=A,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=AN_nfe,Number=1,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=nhomalt_nfe,Number=A,Type=Integer,Description=\"\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    raw = path.with_suffix("")
    raw.write_text(header + "".join(r if r.endswith("\n") else r + "\n"
                                    for r in records))
    gz = pysam.tabix_index(str(raw), preset="vcf", force=True)
    return gz


def test_build_roundtrip_synthetic_vcf():
    """Build a tiny bgzipped+tabix VCF, slice it via the freq builder's parse,
    and assert the produced payload is the correct GraphQL shape — the
    parse+read roundtrip the real (huge) VCF would exercise."""
    tmpdir = Path(tempfile.mkdtemp(prefix="gnomad_vcf_test_"))
    gz = _write_synthetic_vcf(
        tmpdir / "tiny.vcf.gz",
        [
            "chr14\t23412711\trs12147570\tG\tT\t.\tPASS\t"
            "AC=19395;AN=152196;AF=0.127434;nhomalt=1362;AC_XY=0;"
            "fafmax_faf95_max=0.155552;fafmax_faf95_max_gen_anc=nfe;"
            "AC_nfe=10000;AN_nfe=68042;nhomalt_nfe=800",
            "chr14\t23412712\trs1\tC\tT\t.\tPASS\t"
            "AC=1;AN=152180;AF=7e-06;nhomalt=0;AC_XY=0",
        ],
    )
    tbx = pysam.TabixFile(gz)
    halves: dict[str, dict] = {}
    for line in tbx.fetch("chr14", 23412710, 23412713):
        freq_build._accumulate(line, halves, "genome")
    tbx.close()

    assert set(halves) == {"14-23412711-G-T", "14-23412712-C-T"}, halves

    common = halves["14-23412711-G-T"]
    payload = freq_build.merge_payload(
        "14-23412711-G-T", common.get("exome"), common.get("genome")
    )
    assert set(payload) == {"variantId", "rsid", "exome", "genome"}, payload
    assert payload["variantId"] == "14-23412711-G-T"
    assert payload["rsid"] == "rs12147570"
    assert payload["exome"] is None
    g = payload["genome"]
    assert set(g) == {
        "ac", "an", "af", "ac_hom", "ac_hemi", "filters", "populations", "faf95"
    }, g
    assert g["ac"] == 19395 and g["an"] == 152196
    assert abs(g["af"] - 0.127434) < 1e-9
    assert g["ac_hom"] == 1362
    assert g["ac_hemi"] == 0
    assert g["filters"] == []
    assert g["faf95"] == {"popmax": 0.155552, "popmax_population": "nfe"}
    assert g["populations"] == [
        {"id": "nfe", "ac": 10000, "an": 68042, "ac_hom": 800, "ac_hemi": 0},
    ], g["populations"]

    rare_half = halves["14-23412712-C-T"]
    rare = freq_build.merge_payload(
        "14-23412712-C-T", rare_half.get("exome"), rare_half.get("genome")
    )
    assert rare["genome"]["faf95"] == {"popmax": None, "popmax_population": None}
    assert rare["genome"]["populations"] == []
    print("test_build_roundtrip_synthetic_vcf: PASS")


def test_chrx_nonpar_hemi_maps_ac_xy():
    """On chrX non-PAR, AC_XY IS the hemizygote count and maps to ac_hemi; on
    autosomes the same AC_XY value must NOT (it's a male allele count). This
    locks the BS2-critical sex-chrom gating."""
    vid_x, half_x = freq_build.vcf_record_to_payload_fields(
        "chrX", 100000137, "rsX", "A", "C", "PASS",
        "AC=14;AN=183000;AF=7.6e-05;nhomalt=0;AC_XY=3;AC_nfe=5;AN_nfe=80000;"
        "nhomalt_nfe=0;AC_XY_nfe=2",
    )
    block_x = half_x["block"]
    assert vid_x == "X-100000137-A-C"
    assert block_x["ac_hemi"] == 3, block_x
    assert block_x["populations"][0]["ac_hemi"] == 2, block_x

    _vid_a, half_a = freq_build.vcf_record_to_payload_fields(
        "chr1", 5000000, "rsA", "A", "C", "PASS",
        "AC=14;AN=183000;AF=7.6e-05;nhomalt=0;AC_XY=3;AC_nfe=5;AN_nfe=80000;"
        "nhomalt_nfe=0;AC_XY_nfe=2",
    )
    block_a = half_a["block"]
    assert block_a["ac_hemi"] == 0, block_a
    assert block_a["populations"][0]["ac_hemi"] == 0, block_a

    _vid_par, half_par = freq_build.vcf_record_to_payload_fields(
        "chrX", 1000000, "rsPar", "A", "C", "PASS",
        "AC=14;AN=183000;AF=7.6e-05;nhomalt=0;AC_XY=3",
    )
    assert half_par["block"]["ac_hemi"] == 0, half_par["block"]
    print("test_chrx_nonpar_hemi_maps_ac_xy: PASS")


def test_build_db_end_to_end_via_db_then_client():
    """Drive merge_payload → gnomad_freq.db → fetch_gnomad to prove the FULL
    build→read chain produces a row the client serves byte-identically."""
    payload = freq_build.merge_payload(
        "14-23412711-G-T",
        None,
        {"rsid": "rs12147570", "block": {
            "ac": 19395, "an": 152196, "af": 0.127434, "ac_hom": 1362,
            "ac_hemi": 0, "filters": [],
            "populations": [
                {"id": "nfe", "ac": 10000, "an": 68042, "ac_hom": 800, "ac_hemi": None},
            ],
            "faf95": {"popmax": 0.155552, "popmax_population": "nfe"},
        }},
    )
    freq_db = _make_freq_db([("14-23412711-G-T", payload)])
    with _EnvDBs(freq_db=freq_db, constraint_db=None):
        _set_panel({"14": [(23000000, 24000000)]})
        restore = _block_live()
        async def _gene_block(gene):  # noqa: ANN001
            return {"gene_id": "ENSG00000092054", "symbol": gene,
                    "gnomad_constraint": {"pLI": 0.0, "oe_lof": 0.5,
                                          "oe_lof_upper": 0.6, "oe_mis": 0.7,
                                          "oe_mis_upper": 0.8, "mis_z": 3.0,
                                          "syn_z": 0.1}}
        gnomad._fetch_constraint_block_live = _gene_block
        try:
            res = asyncio.run(gnomad.fetch_gnomad("14-23412711-G-T", "MYH7"))
        finally:
            restore()
    assert res["ok"] is True and res["variant_found"] is True, res
    assert res["variant"] == payload, res
    print("test_build_db_end_to_end_via_db_then_client: PASS")


def test_local_freq_hit_returns_exact_blob():
    """A variant present in the local freq DB returns variant_found=True + the
    EXACT stored payload, with constraint served locally and NO live call."""
    freq_db = _make_freq_db([("14-23412711-G-T", _PAYLOAD)])
    constraint_db = _make_constraint_db([
        ("MYH7", "ENSG00000092054", 1.68e-19, 0.644, 0.662, 0.655, 0.682, 7.37, 0.68),
    ])
    with _EnvDBs(freq_db=freq_db, constraint_db=constraint_db):
        _set_panel({"14": [(23000000, 24000000)]})
        restore = _block_live()
        try:
            res = asyncio.run(gnomad.fetch_gnomad("14-23412711-G-T", "MYH7"))
        finally:
            restore()
    assert res["ok"] is True, res
    assert res["variant_found"] is True, res
    assert res["variant_id"] == "14-23412711-G-T", res
    assert res["variant"] == _PAYLOAD, res
    assert res["errors"] is None, res
    assert res["gene"]["gene_id"] == "ENSG00000092054", res
    assert res["gene"]["symbol"] == "MYH7", res
    c = res["gene"]["gnomad_constraint"]
    assert set(c) == {"pLI", "oe_lof", "oe_lof_upper", "oe_mis", "oe_mis_upper",
                      "mis_z", "syn_z"}, c
    assert c["oe_lof"] == 0.644 and c["oe_lof_upper"] == 0.662 and c["mis_z"] == 7.37, c
    assert set(res) == {"ok", "variant_id", "variant_found", "variant",
                        "gene", "errors"}, res
    print("test_local_freq_hit_returns_exact_blob: PASS")


def test_in_panel_miss_is_absent_no_network():
    """A variant inside a panel interval but absent from an existing freq DB
    returns variant_found=False / ok=True with NO live call (drives PM2)."""
    freq_db = _make_freq_db([("14-23412711-G-T", _PAYLOAD)])
    constraint_db = _make_constraint_db([
        ("MYH7", "ENSG00000092054", 1.68e-19, 0.644, 0.662, 0.655, 0.682, 7.37, 0.68),
    ])
    with _EnvDBs(freq_db=freq_db, constraint_db=constraint_db):
        _set_panel({"14": [(23000000, 24000000)]})
        restore = _block_live()
        try:
            res = asyncio.run(gnomad.fetch_gnomad("14-23412999-A-G", "MYH7"))
        finally:
            restore()
    assert res["ok"] is True, res
    assert res["variant_found"] is False, res
    assert res["variant"] is None, res
    assert res["variant_id"] == "14-23412999-A-G", res
    assert res["gene"]["gene_id"] == "ENSG00000092054", res
    print("test_in_panel_miss_is_absent_no_network: PASS")


def _make_live_client(variant_found: bool = True, by_rsid: bool = False):
    """AsyncClient wired to a MockTransport emulating the gnomAD GraphQL
    endpoint the live fallback POSTs to."""
    variant_blob = {
        "variantId": "1-99999-A-G",
        "rsid": "rs999" if not by_rsid else "rs777",
        "exome": None,
        "genome": {
            "ac": 5, "an": 100000, "af": 5e-05, "ac_hom": 0, "ac_hemi": 0,
            "filters": [], "populations": [],
            "faf95": {"popmax": 1e-05, "popmax_population": "nfe"},
        },
    }
    gene_blob = {
        "gene_id": "ENSG00000099999", "symbol": "OFFGENE",
        "gnomad_constraint": {"pLI": 0.9, "oe_lof": 0.3, "oe_lof_upper": 0.4,
                              "oe_mis": 0.8, "oe_mis_upper": 0.9, "mis_z": 2.0,
                              "syn_z": 0.0},
    }

    def handler(request: httpx.Request) -> httpx.Response:
        body = json.loads(request.content.decode())
        q = body.get("query", "")
        variables = body.get("variables", {})
        if "VariantByRsid" in q or "rsid:" in q:
            return httpx.Response(200, json={"data": {"variant": variant_blob}})
        data = {
            "variant": variant_blob if variant_found else None,
            "gene": gene_blob,
        }
        if "variantId" not in variables:
            data = {"gene": gene_blob}
        return httpx.Response(200, json={"data": data})

    return _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler))


def test_off_panel_falls_back_to_live():
    """A variant OUTSIDE the panel (DB present) routes to the live query."""
    freq_db = _make_freq_db([("14-23412711-G-T", _PAYLOAD)])
    with _EnvDBs(freq_db=freq_db, constraint_db=None):
        _set_panel({"14": [(23000000, 24000000)]})
        saved_client = gnomad.httpx.AsyncClient
        gnomad.httpx.AsyncClient = lambda *a, **k: _make_live_client(variant_found=True)
        try:
            res = asyncio.run(gnomad.fetch_gnomad("1-99999-A-G", "OFFGENE"))
        finally:
            gnomad.httpx.AsyncClient = saved_client
    assert res["ok"] is True and res["variant_found"] is True, res
    assert res["variant"]["variantId"] == "1-99999-A-G", res
    assert res["gene"]["gene_id"] == "ENSG00000099999", res
    print("test_off_panel_falls_back_to_live: PASS")


def test_db_absent_falls_back_to_live():
    """With the freq DB absent, any variant routes to the live query."""
    with _EnvDBs(freq_db=Path("/nonexistent/heartvar/gnomad_freq.db"),
                 constraint_db=None):
        _set_panel({"14": [(23000000, 24000000)]})
        saved_client = gnomad.httpx.AsyncClient
        gnomad.httpx.AsyncClient = lambda *a, **k: _make_live_client(variant_found=True)
        try:
            res = asyncio.run(gnomad.fetch_gnomad("14-23412711-G-T", "MYH7"))
        finally:
            gnomad.httpx.AsyncClient = saved_client
    assert res["ok"] is True and res["variant_found"] is True, res
    assert res["variant"]["variantId"] == "1-99999-A-G", res
    print("test_db_absent_falls_back_to_live: PASS")


def test_panel_unavailable_falls_back_to_live():
    """When the panel can't be resolved, an in-DB-miss must NOT be reported as
    absent — membership is undeterminable, so we fail safe to live."""
    freq_db = _make_freq_db([("14-23412711-G-T", _PAYLOAD)])
    with _EnvDBs(freq_db=freq_db, constraint_db=None):
        _set_panel(None)
        saved_client = gnomad.httpx.AsyncClient
        gnomad.httpx.AsyncClient = lambda *a, **k: _make_live_client(variant_found=True)
        try:
            res = asyncio.run(gnomad.fetch_gnomad("14-99999999-A-G", "MYH7"))
        finally:
            gnomad.httpx.AsyncClient = saved_client
    assert res["ok"] is True and res["variant_found"] is True, res
    assert res["variant"]["variantId"] == "1-99999-A-G", res
    print("test_panel_unavailable_falls_back_to_live: PASS")


def test_indel_rsid_fallback_stays_live():
    """The indel/dup rsID fallback path is preserved in the live fallback:
    direct lookup misses, the rsID re-query resolves it, resolved_by='rsid'."""
    with _EnvDBs(freq_db=Path("/nonexistent/heartvar/gnomad_freq.db"),
                 constraint_db=None):
        _set_panel({"14": [(23000000, 24000000)]})
        saved_client = gnomad.httpx.AsyncClient
        gnomad.httpx.AsyncClient = lambda *a, **k: _make_live_client(
            variant_found=False, by_rsid=True)

        import backend.clients.ensembl_vep as vep
        saved_recoder = vep.fetch_variant_recoder_rsid

        async def _fake_recoder(hgvs):  # noqa: ANN001
            return {"rsid": "rs777", "spdi": "NC_000014.9:1000:A:AG"}

        vep.fetch_variant_recoder_rsid = _fake_recoder
        try:
            res = asyncio.run(gnomad.fetch_gnomad(
                "14-99999-A-AG", "MYH7", indel_hgvs="NM_000257.4:c.100dupA"))
        finally:
            gnomad.httpx.AsyncClient = saved_client
            vep.fetch_variant_recoder_rsid = saved_recoder
    assert res["ok"] is True and res["variant_found"] is True, res
    assert res["resolved_by"] == "rsid", res
    assert res["indel_rsid"] == "rs777", res
    print("test_indel_rsid_fallback_stays_live: PASS")


def test_in_panel_indel_miss_recovers_via_rsid():
    """NEW: an in-panel indel missing from the local freq slice is no longer
    declared absent — the rsID fallback now runs on the LOCAL path too (was
    live-only) and recovers the frequency (fixes the BA1/PM2 indel coupling)."""
    freq_db = _make_freq_db([("14-23412711-G-T", _PAYLOAD)])
    with _EnvDBs(freq_db=freq_db, constraint_db=None):
        _set_panel({"14": [(23000000, 24000000)]})
        saved_client = gnomad.httpx.AsyncClient
        gnomad.httpx.AsyncClient = lambda *a, **k: _make_live_client(by_rsid=True)
        import backend.clients.ensembl_vep as vep
        saved_recoder = vep.fetch_variant_recoder_rsid

        async def _fake_recoder(hgvs):  # noqa: ANN001
            return {"rsid": "rs888", "spdi": "NC_000014.9:23412998:AAA:A"}

        vep.fetch_variant_recoder_rsid = _fake_recoder
        try:
            res = asyncio.run(gnomad.fetch_gnomad(
                "14-23412998-AAA-A", "MYH7", indel_hgvs="NM_000257.4:c.100_102del"))
        finally:
            gnomad.httpx.AsyncClient = saved_client
            vep.fetch_variant_recoder_rsid = saved_recoder
    assert res["ok"] is True and res["variant_found"] is True, res
    assert res["resolved_by"] == "rsid" and res["indel_rsid"] == "rs888", res
    print("test_in_panel_indel_miss_recovers_via_rsid: PASS")


def test_in_panel_indel_unresolved_marks_flag():
    """An in-panel indel with no resolvable rsID is flagged indel_unresolved
    (frequency UNKNOWN) instead of silently 'absent' — so the PM2 gate won't
    over-call. No rsID → no FREQUENCY network call.

    ⚠ A CONSTRAINT DB IS SUPPLIED, and it has to be. This fixture passed
    constraint_db=None while _block_live() tripwires the constraint fallback as
    well as the frequency one — so the assertion the test exists to make could
    never be reached: an in-panel frequency miss legitimately fetches gene
    constraint (gnomad.py:447), the tripwire fired first, and the failure read
    as "a live frequency call leaked" when nothing about the frequency path was
    wrong. The constraint live call is CORRECT behaviour with no local DB; what
    this test is about is the frequency path, so the fixture now serves
    constraint locally and leaves the tripwire meaning what it says."""
    freq_db = _make_freq_db([("14-23412711-G-T", _PAYLOAD)])
    constraint_db = _make_constraint_db([
        ("MYH7", "ENSG00000092054", 1.68e-19, 0.644, 0.662, 0.655, 0.682, 7.37, 0.68),
    ])
    with _EnvDBs(freq_db=freq_db, constraint_db=constraint_db):
        _set_panel({"14": [(23000000, 24000000)]})
        restore = _block_live()
        import backend.clients.ensembl_vep as vep
        saved_recoder = vep.fetch_variant_recoder_rsid

        async def _fake_recoder(hgvs):  # noqa: ANN001
            return {"rsid": None, "spdi": "NC_000014.9:23412998:AAA:A"}

        vep.fetch_variant_recoder_rsid = _fake_recoder
        try:
            res = asyncio.run(gnomad.fetch_gnomad(
                "14-23412998-AAA-A", "MYH7", indel_hgvs="NM_000257.4:c.100_102del"))
        finally:
            restore()
            vep.fetch_variant_recoder_rsid = saved_recoder
    assert res["variant_found"] is False and res.get("indel_unresolved") is True, res
    print("test_in_panel_indel_unresolved_marks_flag: PASS")


_FBN1_CONSTRAINT = [
    ("FBN1", "ENSG00000166147", 1.0, 0.1, 0.2, 0.7, 0.8, 3.5, 0.1),
]


def test_canonical_key_hit_returns_variant():
    """A canonical (vcf_string-derived) indel key that HITS the local freq DB
    returns variant_found=True with NO network — the canonical flag never
    perturbs a hit."""
    freq_db = _make_freq_db([("15-48425438-TG-T", _PAYLOAD)])
    constraint_db = _make_constraint_db(_FBN1_CONSTRAINT)
    with _EnvDBs(freq_db=freq_db, constraint_db=constraint_db):
        _set_panel({"15": [(48000000, 49000000)]})
        restore = _block_live()
        try:
            res = asyncio.run(gnomad.fetch_gnomad(
                "15-48425438-TG-T", "FBN1",
                indel_hgvs="NM_000138.5:c.7383del",
                variant_id_canonical=True))
        finally:
            restore()
    assert res["ok"] is True and res["variant_found"] is True, res
    assert res["variant"] == _PAYLOAD, res
    print("test_canonical_key_hit_returns_variant: PASS")


def test_canonical_indel_miss_is_confirmed_absent():
    """Canonical-key in-panel indel MISS + recoder SUCCEEDS reporting NO rsID →
    the indel is declared GENUINELY ABSENT (indel_unresolved NOT set), so
    PM2_Supporting can fire. This is the core WI2 recovery."""
    freq_db = _make_freq_db([("15-99999999-A-G", _PAYLOAD)])
    constraint_db = _make_constraint_db(_FBN1_CONSTRAINT)
    with _EnvDBs(freq_db=freq_db, constraint_db=constraint_db):
        _set_panel({"15": [(48000000, 49000000)]})
        restore = _block_live()
        import backend.clients.ensembl_vep as vep
        saved_recoder = vep.fetch_variant_recoder_rsid

        async def _fake_recoder(hgvs):  # noqa: ANN001
            return {"ok": True, "rsid": None, "spdi": "NC_000015.10:48425437:TG:T"}

        vep.fetch_variant_recoder_rsid = _fake_recoder
        try:
            res = asyncio.run(gnomad.fetch_gnomad(
                "15-48425438-TG-T", "FBN1",
                indel_hgvs="NM_000138.5:c.7383del",
                variant_id_canonical=True))
        finally:
            restore()
            vep.fetch_variant_recoder_rsid = saved_recoder
    assert res["variant_found"] is False, res
    assert res.get("indel_unresolved") is not True, res
    print("test_canonical_indel_miss_is_confirmed_absent: PASS")


def test_canonical_indel_miss_failed_recoder_stays_unresolved():
    """SAFETY: a FAILED recoder call (ok=False) is NOT evidence of absence —
    even with a canonical key the indel stays UNRESOLVED. Without this a common
    indel hitting a failed/offline recoder call would be wrongly declared
    absent → wrong PM2 + wrong BA1/BS1 suppression = catastrophic FP."""
    freq_db = _make_freq_db([("15-99999999-A-G", _PAYLOAD)])
    constraint_db = _make_constraint_db(_FBN1_CONSTRAINT)
    with _EnvDBs(freq_db=freq_db, constraint_db=constraint_db):
        _set_panel({"15": [(48000000, 49000000)]})
        restore = _block_live()
        import backend.clients.ensembl_vep as vep
        saved_recoder = vep.fetch_variant_recoder_rsid

        async def _fake_recoder(hgvs):  # noqa: ANN001
            return {"ok": False, "rsid": None, "spdi": None}

        vep.fetch_variant_recoder_rsid = _fake_recoder
        try:
            res = asyncio.run(gnomad.fetch_gnomad(
                "15-48425438-TG-T", "FBN1",
                indel_hgvs="NM_000138.5:c.7383del",
                variant_id_canonical=True))
        finally:
            restore()
            vep.fetch_variant_recoder_rsid = saved_recoder
    assert res["variant_found"] is False and res.get("indel_unresolved") is True, res
    print("test_canonical_indel_miss_failed_recoder_stays_unresolved: PASS")


def test_noncanonical_indel_miss_stays_unresolved():
    """A NON-canonical hand-built indel key can never match gnomAD's anchored
    key, so its miss proves nothing — the indel stays UNRESOLVED even when the
    recoder succeeds with no rsID (variant_id_canonical defaults False)."""
    freq_db = _make_freq_db([("15-99999999-A-G", _PAYLOAD)])
    constraint_db = _make_constraint_db(_FBN1_CONSTRAINT)
    with _EnvDBs(freq_db=freq_db, constraint_db=constraint_db):
        _set_panel({"15": [(48000000, 49000000)]})
        restore = _block_live()
        import backend.clients.ensembl_vep as vep
        saved_recoder = vep.fetch_variant_recoder_rsid

        async def _fake_recoder(hgvs):  # noqa: ANN001
            return {"ok": True, "rsid": None, "spdi": "NC_000015.10:48425437:TG:T"}

        vep.fetch_variant_recoder_rsid = _fake_recoder
        try:
            res = asyncio.run(gnomad.fetch_gnomad(
                "15-48425439-G-GA", "FBN1",
                indel_hgvs="NM_000138.5:c.7383del"))
        finally:
            restore()
            vep.fetch_variant_recoder_rsid = saved_recoder
    assert res["variant_found"] is False and res.get("indel_unresolved") is True, res
    print("test_noncanonical_indel_miss_stays_unresolved: PASS")


def test_constraint_only_local():
    """The constraint-only path (variant_id=None) serves the gene block from
    the local constraint DB with NO live call."""
    constraint_db = _make_constraint_db([
        ("TTN", "ENSG00000155657", 1.39e-89, 0.479, 0.490, 0.99, 1.01, 7.05, 0.5),
    ])
    with _EnvDBs(freq_db=None, constraint_db=constraint_db):
        restore = _block_live()
        try:
            res = asyncio.run(gnomad.fetch_gnomad(None, "ttn"))
        finally:
            restore()
    assert res["ok"] is True, res
    assert res["variant_id"] is None and res["variant_found"] is False, res
    assert res["variant"] is None, res
    assert res["gene"]["gene_id"] == "ENSG00000155657", res
    assert res["gene"]["symbol"] == "TTN", res
    assert res["gene"]["gnomad_constraint"]["oe_lof"] == 0.479, res
    assert res["gene"]["gnomad_constraint"]["oe_lof_upper"] == 0.490, res
    assert set(res) == {"ok", "variant_id", "variant_found", "variant", "gene"}, res
    print("test_constraint_only_local: PASS")


def test_constraint_only_live_fallback_when_db_absent():
    """Constraint-only with the constraint DB absent falls back to the live
    constraint query, preserving the original shape."""
    with _EnvDBs(freq_db=None, constraint_db=Path("/nonexistent/c.db")):
        saved_client = gnomad.httpx.AsyncClient
        gnomad.httpx.AsyncClient = lambda *a, **k: _make_live_client()
        try:
            res = asyncio.run(gnomad.fetch_gnomad(None, "OFFGENE"))
        finally:
            gnomad.httpx.AsyncClient = saved_client
    assert res["ok"] is True, res
    assert res["variant_id"] is None and res["variant_found"] is False, res
    assert res["gene"]["gene_id"] == "ENSG00000099999", res
    print("test_constraint_only_live_fallback_when_db_absent: PASS")


def test_constraint_gene_not_in_db_falls_back_live():
    """A gene absent from the local constraint DB (constraint-only path) falls
    through to the live constraint query rather than returning a null gene."""
    constraint_db = _make_constraint_db([
        ("MYH7", "ENSG00000092054", 1.68e-19, 0.644, 0.662, 0.655, 0.682, 7.37, 0.68),
    ])
    with _EnvDBs(freq_db=None, constraint_db=constraint_db):
        saved_client = gnomad.httpx.AsyncClient
        gnomad.httpx.AsyncClient = lambda *a, **k: _make_live_client()
        try:
            res = asyncio.run(gnomad.fetch_gnomad(None, "OFFGENE"))
        finally:
            gnomad.httpx.AsyncClient = saved_client
    assert res["ok"] is True and res["gene"]["gene_id"] == "ENSG00000099999", res
    print("test_constraint_gene_not_in_db_falls_back_live: PASS")


def _reconcile(result, indel_hgvs, canonical=False):
    return asyncio.run(
        gnomad._resolve_indel_freq_by_rsid(result, indel_hgvs, canonical))


def test_the_recoder_is_skipped_when_vep_is_offline():
    import backend.clients.ensembl_vep as vep
    saved = vep.fetch_variant_recoder_rsid
    prior = os.environ.get("HEARTVAR_VEP_OFFLINE")
    os.environ["HEARTVAR_VEP_OFFLINE"] = "1"

    async def _must_not_call(hgvs):  # noqa: ANN001
        raise AssertionError(
            "reached Ensembl variant_recoder with HEARTVAR_VEP_OFFLINE set")

    vep.fetch_variant_recoder_rsid = _must_not_call
    try:
        out = _reconcile({"variant_found": False, "ok": True},
                         "ENST00000423902.7:c.5058del")
    finally:
        vep.fetch_variant_recoder_rsid = saved
        if prior is None:
            os.environ.pop("HEARTVAR_VEP_OFFLINE", None)
        else:
            os.environ["HEARTVAR_VEP_OFFLINE"] = prior
    assert out["indel_unresolved"] is True, out
    assert out.get("variant_found") is not True, out
    print("test_the_recoder_is_skipped_when_vep_is_offline: PASS")


def test_the_recoder_is_still_used_when_the_annotation_path_is_live():
    """The gate is narrow on purpose — it is NOT HEARTVAR_OFFLINE_STRICT, which
    stays off because twelve clients want their live fallback on a local miss."""
    import backend.clients.ensembl_vep as vep
    saved = vep.fetch_variant_recoder_rsid
    prior = os.environ.get("HEARTVAR_VEP_OFFLINE")
    os.environ.pop("HEARTVAR_VEP_OFFLINE", None)
    called = {}

    async def _recoder(hgvs):  # noqa: ANN001
        called["hgvs"] = hgvs
        return {"ok": True, "rsid": None, "spdi": None}

    vep.fetch_variant_recoder_rsid = _recoder
    try:
        out = _reconcile({"variant_found": False, "ok": True},
                         "ENST00000423902.7:c.5058del", canonical=True)
    finally:
        vep.fetch_variant_recoder_rsid = saved
        if prior is not None:
            os.environ["HEARTVAR_VEP_OFFLINE"] = prior
    assert called["hgvs"] == "ENST00000423902.7:c.5058del", called
    assert out.get("indel_unresolved") is not True, out
    print("test_the_recoder_is_still_used_when_the_annotation_path_is_live: PASS")


if __name__ == "__main__":
    test_build_roundtrip_synthetic_vcf()
    test_chrx_nonpar_hemi_maps_ac_xy()
    test_build_db_end_to_end_via_db_then_client()
    test_local_freq_hit_returns_exact_blob()
    test_in_panel_miss_is_absent_no_network()
    test_off_panel_falls_back_to_live()
    test_db_absent_falls_back_to_live()
    test_panel_unavailable_falls_back_to_live()
    test_indel_rsid_fallback_stays_live()
    test_in_panel_indel_miss_recovers_via_rsid()
    test_in_panel_indel_unresolved_marks_flag()
    test_canonical_key_hit_returns_variant()
    test_canonical_indel_miss_is_confirmed_absent()
    test_canonical_indel_miss_failed_recoder_stays_unresolved()
    test_noncanonical_indel_miss_stays_unresolved()
    test_constraint_only_local()
    test_constraint_only_live_fallback_when_db_absent()
    test_constraint_gene_not_in_db_falls_back_live()
    test_the_recoder_is_skipped_when_vep_is_offline()
    test_the_recoder_is_still_used_when_the_annotation_path_is_live()
    print("ALL GNOMAD TESTS PASS")
