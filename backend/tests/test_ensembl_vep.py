"""Unit tests for backend.clients.ensembl_vep — the FIRST external call on
every curation, previously untested.

Two layers, both fully offline (no network, no DB):

  * Pure-function tests imported directly: _max_dbnsfp, _is_refseq,
    _strip_version, _pick_transcript_consequence, _complement,
    gnomad_variant_id, _vcf_input_for_region, strip_transcript_prefix.
    (parse_variant_input is intentionally NOT re-tested — it is covered by
    test_input_and_phenotype_parsing.py.)

  * HTTP-shape tests for fetch_vep / fetch_vep_by_coordinates driven through
    an httpx.MockTransport (mirrors test_pmcoa.py): the real
    httpx.AsyncClient is captured up front, ensembl_vep.httpx.AsyncClient is
    monkeypatched to a MockTransport-backed client keyed on str(request.url)
    substrings, and restored in a finally block. Each async entry point is
    driven with asyncio.run.

No pytest-only fixtures — runnable directly
(``python -m backend.tests.test_ensembl_vep``).
"""
from __future__ import annotations

import asyncio

import httpx

import backend.clients.ensembl_vep as vep
from backend.clients.ensembl_vep import (
    _max_dbnsfp,
    _is_refseq,
    _strip_version,
    _pick_transcript_consequence,
    _build_transcript_table,
    _complement,
    gnomad_variant_id,
    gnomad_variant_id_with_provenance,
    _id_from_vcf_string,
    _vcf_input_for_region,
    strip_transcript_prefix,
)

_REAL_ASYNC_CLIENT = httpx.AsyncClient


def test_max_dbnsfp_comma_string_takes_max():
    assert _max_dbnsfp("0.6,0.9") == 0.9
    assert _max_dbnsfp(".,0.651,0.651,.,.") == 0.651
    assert _max_dbnsfp("0.10, 0.40 ,0.25") == 0.40


def test_max_dbnsfp_none_and_dot_and_empty_return_none():
    assert _max_dbnsfp(None) is None
    assert _max_dbnsfp(".") is None
    assert _max_dbnsfp("") is None
    assert _max_dbnsfp(".,.,.") is None


def test_max_dbnsfp_plain_number_passthrough():
    assert _max_dbnsfp(0.7) == 0.7
    assert _max_dbnsfp(1) == 1.0
    assert isinstance(_max_dbnsfp(1), float)


def test_max_dbnsfp_non_numeric_tokens_skipped():
    assert _max_dbnsfp("foo,0.3,bar") == 0.3
    assert _max_dbnsfp("foo,bar") is None


def test_is_refseq_true_for_refseq_accessions():
    assert _is_refseq("NM_000257.4") is True
    assert _is_refseq("NR_001234") is True
    assert _is_refseq("XM_005343") is True
    assert _is_refseq("XR_999.1") is True
    assert _is_refseq("nm_000257.4") is True


def test_is_refseq_false_for_enst_and_bare_and_none():
    assert _is_refseq("ENST00000355349") is False
    assert _is_refseq("ENST00000355349.4") is False
    assert _is_refseq(None) is False
    assert _is_refseq("") is False
    assert _is_refseq("NM000257") is False


def test_strip_version():
    assert _strip_version("NM_000257.4") == "NM_000257"
    assert _strip_version("ENST00000355349.4") == "ENST00000355349"
    assert _strip_version("NM_000257") == "NM_000257"
    assert _strip_version(None) is None


def test_pick_transcript_empty_returns_empty_dict():
    assert _pick_transcript_consequence([], "NM_000257.4") == {}
    assert _pick_transcript_consequence([], None) == {}


def test_pick_transcript_supplied_base_version_match_wins():
    consequences = [
        {"transcript_id": "ENST00000000001", "mane_select": "NM_111111.1"},
        {"transcript_id": "NM_000257.3"},
    ]
    picked = _pick_transcript_consequence(consequences, "NM_000257.4")
    assert picked["transcript_id"] == "NM_000257.3"


def test_pick_transcript_mane_select_fallback():
    consequences = [
        {"transcript_id": "ENST00000000009", "canonical": 1},
        {"transcript_id": "ENST00000355349", "mane_select": "NM_000257.4"},
    ]
    picked = _pick_transcript_consequence(consequences, None)
    assert picked["transcript_id"] == "ENST00000355349"


def test_pick_transcript_canonical_then_first():
    consequences = [
        {"transcript_id": "ENST_A"},
        {"transcript_id": "ENST_B", "canonical": 1},
    ]
    assert _pick_transcript_consequence(consequences, None)["transcript_id"] == "ENST_B"
    plain = [{"transcript_id": "ENST_X"}, {"transcript_id": "ENST_Y"}]
    assert _pick_transcript_consequence(plain, None)["transcript_id"] == "ENST_X"


def test_pick_transcript_prefers_requested_gene_over_a_neighbours_mane():
    consequences = [
        {"transcript_id": "ENST00000358528", "gene_symbol": "CSDE1",
         "mane_select": "NM_001007553.3", "canonical": 1,
         "consequence_terms": ["downstream_gene_variant"]},
        {"transcript_id": "ENST00000369535", "gene_symbol": "NRAS",
         "mane_select": "NM_002524.5", "canonical": 1,
         "hgvsp": "ENSP00000358548.4:p.Gly12Ser",
         "consequence_terms": ["missense_variant"]},
    ]
    picked = _pick_transcript_consequence(consequences, None, gene="NRAS")
    assert picked["gene_symbol"] == "NRAS"
    assert picked["hgvsp"] == "ENSP00000358548.4:p.Gly12Ser"


def test_pick_transcript_gene_filter_is_case_insensitive():
    consequences = [
        {"transcript_id": "ENST_OTHER", "gene_symbol": "CSDE1", "mane_select": "NM_1.1"},
        {"transcript_id": "ENST_WANT", "gene_symbol": "NRAS", "mane_select": "NM_2.1"},
    ]
    picked = _pick_transcript_consequence(consequences, None, gene="nras")
    assert picked["transcript_id"] == "ENST_WANT"


def test_pick_transcript_gene_filter_fails_open():
    consequences = [{"transcript_id": "ENST_A", "gene_symbol": "CSDE1", "canonical": 1}]
    picked = _pick_transcript_consequence(consequences, None, gene="NRAS")
    assert picked["transcript_id"] == "ENST_A"


def test_pick_transcript_supplied_transcript_beats_gene_filter():
    consequences = [
        {"transcript_id": "ENST_NRAS", "gene_symbol": "NRAS", "mane_select": "NM_2.1"},
        {"transcript_id": "NM_001007553.3", "gene_symbol": "CSDE1"},
    ]
    picked = _pick_transcript_consequence(
        consequences, "NM_001007553.3", gene="NRAS")
    assert picked["transcript_id"] == "NM_001007553.3"


def test_pick_transcript_gene_none_preserves_legacy_order():
    consequences = [
        {"transcript_id": "ENST_CANON", "gene_symbol": "CSDE1", "canonical": 1},
        {"transcript_id": "ENST_MANE", "gene_symbol": "NRAS", "mane_select": "NM_2.1"},
    ]
    picked = _pick_transcript_consequence(consequences, None, gene=None)
    assert picked["transcript_id"] == "ENST_MANE"


def test_complement_single_base():
    assert _complement("A") == "T"
    assert _complement("C") == "G"
    assert _complement("G") == "C"
    assert _complement("T") == "A"


def test_complement_multibase_is_reverse_complement():
    assert _complement("AG") == "CT"
    assert _complement("AC") == "GT"
    assert _complement("ACGT") == "ACGT"


def test_gnomad_variant_id_not_ok_returns_none():
    assert gnomad_variant_id({"ok": False}) is None
    assert gnomad_variant_id({}) is None


def test_gnomad_variant_id_forward_passthrough():
    vep_res = {"ok": True, "forward_variant_id": "7-117548628-C-T",
               "seq_region_name": "7", "start": 117548628,
               "allele_string": "C/T", "strand": 1}
    assert gnomad_variant_id(vep_res) == "7-117548628-C-T"


def test_gnomad_variant_id_minus_strand_complemented():
    vep_res = {"ok": True, "seq_region_name": "14", "start": 23429278,
               "allele_string": "C/T", "strand": -1}
    assert gnomad_variant_id(vep_res) == "14-23429278-G-A"


def test_gnomad_variant_id_plus_strand_not_complemented():
    vep_res = {"ok": True, "seq_region_name": "7", "start": 117548628,
               "allele_string": "C/T", "strand": 1}
    assert gnomad_variant_id(vep_res) == "7-117548628-C-T"


def test_gnomad_variant_id_missing_slash_returns_none():
    vep_res = {"ok": True, "seq_region_name": "7", "start": 117548628,
               "allele_string": "CT", "strand": 1}
    assert gnomad_variant_id(vep_res) is None


def test_id_from_vcf_string_list_and_string():
    assert _id_from_vcf_string(["15-48425438-TG-T"]) == "15-48425438-TG-T"
    assert _id_from_vcf_string("15-48425438-TG-T") == "15-48425438-TG-T"
    assert _id_from_vcf_string(["chr15-48425438-TG-T"]) == "15-48425438-TG-T"
    assert _id_from_vcf_string([]) is None
    assert _id_from_vcf_string(None) is None
    assert _id_from_vcf_string("15-48425439-G--") is None
    assert _id_from_vcf_string("15-pos-A-G") is None
    assert _id_from_vcf_string("15-48425438-N-T") is None


def test_gnomad_variant_id_prefers_canonical_vcf_string():
    vep_res = {"ok": True, "vcf_string": ["15-48425438-TG-T"],
               "seq_region_name": "15", "start": 48425439,
               "allele_string": "G/-", "strand": 1}
    vid, canonical = gnomad_variant_id_with_provenance(vep_res)
    assert vid == "15-48425438-TG-T"
    assert canonical is True
    assert gnomad_variant_id(vep_res) == "15-48425438-TG-T"


def test_gnomad_variant_id_snv_vcf_string_matches_complement_no_regression():
    with_vcf = {"ok": True, "vcf_string": ["14-23429278-G-A"],
                "seq_region_name": "14", "start": 23429278,
                "allele_string": "C/T", "strand": -1}
    without_vcf = {"ok": True, "seq_region_name": "14", "start": 23429278,
                   "allele_string": "C/T", "strand": -1}
    vid_new, canonical = gnomad_variant_id_with_provenance(with_vcf)
    assert vid_new == gnomad_variant_id(without_vcf) == "14-23429278-G-A"
    assert canonical is True


def test_gnomad_variant_id_falls_back_when_vcf_string_absent():
    vep_res = {"ok": True, "seq_region_name": "14", "start": 23429278,
               "allele_string": "C/T", "strand": -1}
    vid, canonical = gnomad_variant_id_with_provenance(vep_res)
    assert vid == "14-23429278-G-A"
    assert canonical is False


def test_gnomad_variant_id_garbage_vcf_string_falls_through_noncanonical():
    vep_res = {"ok": True, "vcf_string": ["15-48425439-G--"],
               "forward_variant_id": "15-48425439-G-GA",
               "seq_region_name": "15", "start": 48425439,
               "allele_string": "G/GA", "strand": 1}
    vid, canonical = gnomad_variant_id_with_provenance(vep_res)
    assert vid == "15-48425439-G-GA"
    assert canonical is False


def test_vcf_input_for_region():
    assert _vcf_input_for_region("7", 117548628, "C", "T") == "7 117548628 . C T . . ."
    assert _vcf_input_for_region("X", 100, "A", "G") == "X 100 . A G . . ."


def test_strip_transcript_prefix_refseq_with_version():
    assert strip_transcript_prefix("NM_000257.4:c.1208G>A") == ("c.1208G>A", "NM_000257.4")


def test_strip_transcript_prefix_enst_without_version():
    assert strip_transcript_prefix("ENST00000355349:c.1208G>A") == ("c.1208G>A", "ENST00000355349")


def test_strip_transcript_prefix_colonless_accession():
    assert strip_transcript_prefix("NM_001136239.3c.1301A>T") == ("c.1301A>T", "NM_001136239.3")


def test_strip_transcript_prefix_bare_body_unchanged():
    assert strip_transcript_prefix("c.1208G>A") == ("c.1208G>A", None)
    assert strip_transcript_prefix("") == ("", None)


def _make_client(handler):
    """Build an AsyncClient wired to a MockTransport driving ``handler``.
    Uses the saved real AsyncClient so it survives the monkeypatch."""
    return _REAL_ASYNC_CLIENT(
        transport=httpx.MockTransport(handler), follow_redirects=True
    )


def _run_with_handler(handler, coro_factory):
    """Patch vep.httpx.AsyncClient to a MockTransport-backed client, run the
    coroutine produced by ``coro_factory()`` to completion, restore in finally.
    Each AsyncClient() call inside the production code yields a fresh
    MockTransport client bound to the same handler."""
    orig = vep.httpx.AsyncClient
    vep.httpx.AsyncClient = lambda *a, **k: _make_client(handler)
    try:
        return asyncio.run(coro_factory())
    finally:
        vep.httpx.AsyncClient = orig


def test_fetch_vep_success():
    def handler(request: httpx.Request) -> httpx.Response:
        u = str(request.url)
        if "/vep/human/hgvs" in u:
            return httpx.Response(200, json=[{
                "assembly_name": "GRCh38",
                "seq_region_name": "14",
                "start": 23424081,
                "end": 23424081,
                "allele_string": "G/A",
                "strand": -1,
                "most_severe_consequence": "missense_variant",
                "transcript_consequences": [{
                    "transcript_id": "ENST00000355349",
                    "gene_id": "ENSG00000092054",
                    "gene_symbol": "MYH7",
                    "biotype": "protein_coding",
                    "mane_select": "NM_000257.4",
                    "hgvsc": "ENST00000355349.4:c.1208G>A",
                    "hgvsp": "ENSP00000347507.3:p.Arg403Gln",
                    "consequence_terms": ["missense_variant"],
                    "exon": "13/40",
                    "impact": "MODERATE",
                    "revel_score": "0.6,0.9",
                    "cadd_phred": 28.5,
                }],
            }])
        if "/lookup/symbol" in u:
            return httpx.Response(200, json={"Transcript": []})
        if "/lookup/id" in u:
            return httpx.Response(404, json={})
        return httpx.Response(404, json={})

    res = _run_with_handler(
        handler, lambda: vep.fetch_vep("MYH7", "c.1208G>A")
    )
    assert res["ok"] is True, res
    assert res["most_severe_consequence"] == "missense_variant"
    assert res["hgvsp"] == "ENSP00000347507.3:p.Arg403Gln"
    assert res["revel_score"] == 0.9
    assert res["cadd_phred"] == 28.5
    assert res["gene_symbol"] == "MYH7"
    assert isinstance(res["nmd_escape"], bool)
    assert res["nmd_escape"] is False
    assert "exon_lengths" in res
    assert res["selected_transcript_id"] == "ENST00000355349"
    assert res["is_mane_select"] is True
    assert isinstance(res["transcript_consequences_all"], list)
    assert res["transcript_consequences_all"][0]["transcript_id"] == "ENST00000355349"
    assert res["transcript_consequences_all"][0]["is_picked"] is True
    assert res["consequence_differs_significantly"] is False
    assert res["transcript_set"] == "GENCODE"


def test_fetch_vep_failure_returns_error():
    def handler(request: httpx.Request) -> httpx.Response:
        u = str(request.url)
        if "/vep/human/hgvs" in u:
            return httpx.Response(400, text="Could not get a Transcript object")
        if "/lookup/symbol" in u:
            return httpx.Response(200, json={"Transcript": []})
        return httpx.Response(404, json={})

    res = _run_with_handler(
        handler, lambda: vep.fetch_vep("NOTAGENE", "c.999G>A")
    )
    assert res["ok"] is False, res
    assert res.get("error"), res
    assert "400" in res["error"]


def test_fetch_vep_refseq_versioned_fallback_when_versionless_400s():
    """Some RefSeq accessions (e.g. TPM1 NM_001365781.2) are NOT resolvable
    version-less — VEP 400s "Could not find transcript" — but resolve with the
    exact supplied version. The cascade must try the versioned form before
    falling through to the (wrong-frame) Ensembl canonical. Regression for the
    TPM1 c.215C>T hard-rejection."""
    calls: list[str] = []

    def handler(request: httpx.Request) -> httpx.Response:
        u = str(request.url)
        if "/vep/human/hgvs" in u:
            calls.append(u)
            if "NM_001365781.2" not in u and "NM_001365781" in u:
                return httpx.Response(400, text=(
                    "Unable to parse HGVS notation 'NM_001365781:c.215C>T': "
                    "Could not find transcript"))
            if "NM_001365781.2" in u:
                return httpx.Response(200, json=[{
                    "assembly_name": "GRCh38", "seq_region_name": "15",
                    "start": 63057067, "end": 63057067, "allele_string": "C/T",
                    "strand": 1, "most_severe_consequence": "missense_variant",
                    "transcript_consequences": [{
                        "transcript_id": "NM_001365781.2", "gene_symbol": "TPM1",
                        "biotype": "protein_coding", "source": "RefSeq",
                        "hgvsc": "NM_001365781.2:c.215C>T",
                        "consequence_terms": ["missense_variant"], "impact": "MODERATE",
                    }],
                }])
            return httpx.Response(400, text="reference allele mismatch")
        if "/lookup/symbol" in u:
            return httpx.Response(200, json={"Transcript": []})
        return httpx.Response(404, json={})

    res = _run_with_handler(
        handler,
        lambda: vep.fetch_vep("TPM1", "c.215C>T", supplied_transcript="NM_001365781.2"),
    )
    assert res["ok"] is True, res
    assert res["most_severe_consequence"] == "missense_variant"
    assert any("NM_001365781.2" in u for u in calls), calls


def _mane_only_handler(request: httpx.Request) -> httpx.Response:
    """VEP returns a single MYH7 consequence whose MANE Select is NM_000257.4."""
    u = str(request.url)
    if "/vep/human/hgvs" in u:
        return httpx.Response(200, json=[{
            "assembly_name": "GRCh38", "seq_region_name": "14",
            "start": 23424081, "end": 23424081, "allele_string": "G/A", "strand": -1,
            "most_severe_consequence": "missense_variant",
            "transcript_consequences": [{
                "transcript_id": "ENST00000355349", "gene_id": "ENSG00000092054",
                "gene_symbol": "MYH7", "biotype": "protein_coding",
                "mane_select": "NM_000257.4",
                "hgvsc": "ENST00000355349.4:c.1208G>A",
                "hgvsp": "ENSP00000347507.3:p.Arg403Gln",
                "consequence_terms": ["missense_variant"], "exon": "13/40", "impact": "MODERATE",
            }],
        }])
    if "/lookup/symbol" in u:
        return httpx.Response(200, json={"Transcript": []})
    return httpx.Response(404, json={})


def test_fetch_vep_requested_transcript_honored():
    res = _run_with_handler(
        _mane_only_handler, lambda: vep.fetch_vep("MYH7", "c.1208G>A", "NM_000257.4")
    )
    assert res["ok"] is True, res
    assert res["requested_transcript"] == "NM_000257.4"
    assert res["requested_transcript_honored"] is True


def test_fetch_vep_requested_transcript_not_honored():
    res = _run_with_handler(
        _mane_only_handler, lambda: vep.fetch_vep("MYH7", "c.1208G>A", "NM_999999.9")
    )
    assert res["ok"] is True, res
    assert res["requested_transcript"] == "NM_999999.9"
    assert res["requested_transcript_honored"] is False


def test_fetch_vep_no_requested_transcript_is_none():
    res = _run_with_handler(
        _mane_only_handler, lambda: vep.fetch_vep("MYH7", "c.1208G>A")
    )
    assert res["ok"] is True, res
    assert res["requested_transcript"] is None
    assert res["requested_transcript_honored"] is None


def test_recode_rsid_to_coords_primary_assembly():
    def handler(request: httpx.Request) -> httpx.Response:
        if "/variant_recoder/human" in str(request.url):
            return httpx.Response(200, json=[{"T": {
                "spdi": ["NC_000014.9:23429277:C:T"],
                "vcf_string": ["14-23429278-C-T"],
            }}])
        return httpx.Response(404, json={})
    res = _run_with_handler(handler, lambda: vep.recode_rsid_to_coords("rs727503113"))
    assert res == {"chrom": "14", "pos": 23429278, "ref": "C", "alt": "T"}


def test_recode_rsid_to_coords_skips_alt_contig():
    def handler(request: httpx.Request) -> httpx.Response:
        if "/variant_recoder/human" in str(request.url):
            return httpx.Response(200, json=[{"T": {
                "spdi": ["NW_009646201.1:100:C:T"],
                "vcf_string": ["KI270722.1-100-C-T"],
            }}])
        return httpx.Response(404, json={})
    assert _run_with_handler(handler, lambda: vep.recode_rsid_to_coords("rs999")) is None


def test_recode_rsid_to_coords_not_found_returns_none():
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(404, json={})
    assert _run_with_handler(handler, lambda: vep.recode_rsid_to_coords("rs1")) is None


def test_fetch_vep_by_coordinates_success():
    def handler(request: httpx.Request) -> httpx.Response:
        u = str(request.url)
        if "/vep/human/region" in u:
            assert request.method == "POST"
            return httpx.Response(200, json=[{
                "assembly_name": "GRCh38",
                "seq_region_name": "7",
                "start": 117548628,
                "end": 117548628,
                "allele_string": "C/T",
                "strand": 1,
                "most_severe_consequence": "missense_variant",
                "transcript_consequences": [{
                    "transcript_id": "ENST00000003084",
                    "gene_id": "ENSG00000001626",
                    "gene_symbol": "CFTR",
                    "biotype": "protein_coding",
                    "mane_select": "NM_000492.4",
                    "hgvsc": "ENST00000003084.11:c.1521_1523delCTT",
                    "hgvsp": "ENSP00000003084.6:p.Phe508del",
                    "consequence_terms": ["missense_variant"],
                    "exon": "11/27",
                    "impact": "MODERATE",
                    "revel_score": ".",
                    "cadd_phred": 22.0,
                }],
            }])
        if "/lookup/id" in u:
            return httpx.Response(404, json={})
        return httpx.Response(404, json={})

    res = _run_with_handler(
        handler,
        lambda: vep.fetch_vep_by_coordinates("7", 117548628, "C", "T"),
    )
    assert res["ok"] is True, res
    assert res["input_format"] == "coordinates"
    assert res["input_coords"] == "7-117548628-C-T"
    assert gnomad_variant_id(res) == "7-117548628-C-T"


def test_fetch_vep_by_coordinates_http_failure_preserves_input():
    def handler(request: httpx.Request) -> httpx.Response:
        if "/vep/human/region" in str(request.url):
            return httpx.Response(400, text="bad region input")
        return httpx.Response(404, json={})

    res = _run_with_handler(
        handler,
        lambda: vep.fetch_vep_by_coordinates("7", 117548628, "C", "T",
                                             build="GRCh37"),
    )
    assert res["ok"] is False, res
    assert res["input_format"] == "coordinates"
    assert res["input_build"] == "GRCh37"
    assert res["input_coords"] == "7-117548628-C-T"
    assert res.get("error"), res


def test_fetch_vep_by_coordinates_empty_payload_failure():
    def handler(request: httpx.Request) -> httpx.Response:
        if "/vep/human/region" in str(request.url):
            return httpx.Response(200, json=[])
        return httpx.Response(404, json={})

    res = _run_with_handler(
        handler,
        lambda: vep.fetch_vep_by_coordinates("7", 117548628, "C", "T"),
    )
    assert res["ok"] is False, res
    assert res["input_format"] == "coordinates"
    assert res["input_coords"] == "7-117548628-C-T"


def test_build_transcript_table_sort_flags_and_uniform_impact():
    consequences = [
        {"transcript_id": "ENST00000222222", "biotype": "protein_coding",
         "consequence_terms": ["missense_variant"], "impact": "MODERATE",
         "hgvsc": "ENST00000222222.1:c.10G>A", "hgvsp": "ENSP1:p.Gly4Ser"},
        {"transcript_id": "ENST00000355349", "biotype": "protein_coding",
         "consequence_terms": ["missense_variant"], "impact": "MODERATE",
         "mane_select": "NM_000257.4",
         "hgvsc": "ENST00000355349.4:c.1208G>A", "hgvsp": "ENSP2:p.Arg403Gln"},
        {"transcript_id": "ENST00000333333", "biotype": "protein_coding",
         "consequence_terms": ["missense_variant"], "impact": "MODERATE",
         "mane_plus_clinical": "NM_111111.1"},
    ]
    picked = consequences[1]
    rows, differs, tset = _build_transcript_table(consequences, picked)
    assert rows[0]["transcript_id"] == "ENST00000355349"
    assert rows[0]["is_mane_select"] is True
    assert rows[0]["is_picked"] is True
    assert rows[0]["source"] == "GENCODE"
    assert rows[1]["transcript_id"] == "ENST00000333333"
    assert rows[1]["is_mane_plus_clinical"] is True
    assert tset == "GENCODE"
    assert differs is False
    assert "_rank" not in rows[0]


def test_build_transcript_table_flags_lof_vs_missense():
    consequences = [
        {"transcript_id": "ENST00000000001", "biotype": "protein_coding",
         "consequence_terms": ["missense_variant"], "impact": "MODERATE"},
        {"transcript_id": "ENST00000000002", "biotype": "protein_coding",
         "consequence_terms": ["stop_gained"], "impact": "HIGH"},
    ]
    _rows, differs, _tset = _build_transcript_table(consequences, consequences[0])
    assert differs is True


def test_build_transcript_table_refseq_source_and_impact_fallback():
    consequences = [
        {"transcript_id": "NM_000257.4", "biotype": "protein_coding",
         "consequence_terms": ["synonymous_variant"]},
        {"transcript_id": "NM_999999.1", "biotype": "protein_coding",
         "consequence_terms": ["frameshift_variant"]},
    ]
    rows, differs, tset = _build_transcript_table(consequences, consequences[0])
    assert tset == "RefSeq"
    assert all(r["source"] == "RefSeq" for r in rows)
    assert differs is True


def test_build_transcript_table_single_transcript_never_differs():
    consequences = [
        {"transcript_id": "ENST00000355349", "biotype": "protein_coding",
         "consequence_terms": ["missense_variant"], "impact": "MODERATE"},
    ]
    rows, differs, _tset = _build_transcript_table(consequences, consequences[0])
    assert len(rows) == 1
    assert differs is False


def test_build_transcript_table_empty():
    rows, differs, tset = _build_transcript_table([], {})
    assert rows == []
    assert differs is False
    assert tset is None


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
