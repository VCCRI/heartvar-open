"""Tests for the LOCAL-FIRST ProtVar client.

Covers:
  * Local path against the REAL local caches (uniprot.db + AlphaMissense
    tabix table) — asserts the full return shape, conservation_score/
    pocket_score = None, feature_types sourced from UniProt, and that the
    live EBI fallback is NEVER invoked (httpx is hard-blocked).
  * Local path with monkeypatched canned local data — deterministic shape
    + am_* population, EBI not called.
  * Synthetic AlphaMissense tabix roundtrip — builds a tiny bgzipped +
    tabix-indexed TSV, reads it via the alphamissense client AND drives the
    protvar local path off it, proving the parse/read roundtrip without the
    multi-GB real download.
  * Live fallback — both local sources absent → the original EBI behaviour
    runs against a mocked HTTP transport.

Sync tests drive the coroutine with ``asyncio.run`` (matching test_uniprot.py).

Run:
    .venv/bin/python -m pytest backend/tests/test_protvar.py -q
"""

from __future__ import annotations

import asyncio
from pathlib import Path

import httpx

from backend.clients import protvar

_ROOT = Path(__file__).resolve().parent.parent.parent
_UNIPROT_DB = _ROOT / "data" / "uniprot.db"
_ALPHAMISSENSE = _ROOT / "data" / "AlphaMissense_hg38.tsv.gz"

_FOUND_KEYS = {
    "ok", "applicable", "found", "input", "gene", "uniprot",
    "protein_position", "ref_aa", "alt_aa", "consequence",
    "conservation_score", "feature_types", "function_summary",
    "pocket_score", "colocated_variants", "colocated_total",
    "am_pathogenicity", "am_class", "url",
}


def _myh7_v39m_vep() -> dict:
    """A realistic VEP result for MYH7 p.Val39Met (a residue that carries a
    pathogenic UniProt natural variant + sits inside the Myosin motor domain).

    MYH7 is on chr14, the minus strand — the strand=-1 flag exercises the
    reverse-complement branch in the shared coordinate derivation. The
    genomic ref/alt here are the TRANSCRIPT-strand alleles VEP reports for
    HGVS-c input; the client complements them for the forward-strand input
    string / url.
    """
    return {
        "ok": True,
        "seq_region_name": "14",
        "start": 23433561,
        "allele_string": "C/T",
        "strand": -1,
        "most_severe_consequence": "missense_variant",
        "gene_symbol": "MYH7",
        "hgvsp": "ENSP00000347507.3:p.Val39Met",
        "hgvsc": "ENST00000355349.4:c.115G>A",
    }


def _no_network_transport() -> httpx.MockTransport:
    def _handler(request: httpx.Request) -> httpx.Response:  # pragma: no cover
        raise AssertionError(
            f"unexpected live HTTP call to {request.url} — local path must not "
            "touch the network"
        )
    return httpx.MockTransport(_handler)


def _patch_async_client(monkeypatch, transport: httpx.MockTransport) -> None:
    """Force any ``httpx.AsyncClient()`` constructed inside protvar to ride the
    given MockTransport, so the live fallback (or an unexpected live call) is
    intercepted deterministically."""
    orig = httpx.AsyncClient

    def _factory(*args, **kwargs):
        kwargs["transport"] = transport
        return orig(*args, **kwargs)

    monkeypatch.setattr(protvar.httpx, "AsyncClient", _factory)


def test_local_path_real_caches_no_network(monkeypatch):
    if not _UNIPROT_DB.exists() or not _ALPHAMISSENSE.exists():
        import pytest
        pytest.skip("real local caches (uniprot.db / AlphaMissense) not present")

    monkeypatch.setenv("UNIPROT_DB_PATH", str(_UNIPROT_DB))
    monkeypatch.setenv("ALPHAMISSENSE_PATH", str(_ALPHAMISSENSE))
    _patch_async_client(monkeypatch, _no_network_transport())
    called = {"live": False}
    orig_live = protvar._fetch_protvar_live

    async def _spy(vep):  # pragma: no cover - only runs on failure
        called["live"] = True
        return await orig_live(vep)

    monkeypatch.setattr(protvar, "_fetch_protvar_live", _spy)

    out = asyncio.run(protvar.fetch_protvar(_myh7_v39m_vep()))

    assert called["live"] is False, "live EBI fallback must not run on the local path"
    assert out["ok"] is True
    assert out["applicable"] is True
    assert out["found"] is True
    assert set(out.keys()) == _FOUND_KEYS
    assert out["gene"] == "MYH7"
    assert out["uniprot"] == "P12883"
    assert out["protein_position"] == 39
    assert out["ref_aa"] == "Val"
    assert out["alt_aa"] == "Met"
    assert out["consequence"] == "missense_variant"
    assert out["conservation_score"] is None
    assert out["pocket_score"] is None
    assert "Domain" in out["feature_types"]
    assert out["function_summary"] is None
    assert out["colocated_total"] >= 1
    assert all(
        set(v.keys()) == {"variant", "clinical_significance", "source_db"}
        for v in out["colocated_variants"]
    )
    assert any(v["source_db"] == "UniProt" for v in out["colocated_variants"])
    assert out["url"] == (
        "https://www.ebi.ac.uk/ProtVar/query?search=14%2023433561%20G%20A"
    )
    assert out["input"] == "14 23433561 G A"


def test_local_path_canned_no_network(monkeypatch):
    _patch_async_client(monkeypatch, _no_network_transport())

    async def _canned_am(vep):
        return {
            "available": True,
            "score": 0.912,
            "classification": "likely_pathogenic",
            "protein_variant": "R502W",
            "not_applicable": False,
        }

    async def _canned_uniprot(gene):
        assert gene == "MYH7"
        return {
            "ok": True,
            "found": True,
            "accession": "P12883",
            "features": [
                {"type": "Domain", "start": 85, "end": 778,
                 "description": "Myosin motor"},
                {"type": "Region", "start": 655, "end": 677,
                 "description": "Actin-binding"},
            ],
            "domains": [{"name": "Myosin motor", "start": 85, "end": 778}],
            "natural_variants": [
                {"position": 502, "original_aa": "R", "variant_aa": "Q",
                 "description": "in CMH1; dbSNP:rs1", "clinical_significance": "pathogenic",
                 "dbsnp": "rs1", "clinvar": None, "uniprot_var_id": "VAR_X"},
                {"position": 502, "original_aa": "R", "variant_aa": "W",
                 "description": "in CMH1; dbSNP:rs2", "clinical_significance": "pathogenic",
                 "dbsnp": "rs2", "clinvar": None, "uniprot_var_id": "VAR_Y"},
            ],
        }

    monkeypatch.setattr(protvar, "fetch_alphamissense", _canned_am)
    monkeypatch.setattr(protvar, "fetch_uniprot", _canned_uniprot)

    vep = {
        "ok": True,
        "seq_region_name": "14",
        "start": 23425678,
        "allele_string": "C/T",
        "strand": -1,
        "most_severe_consequence": "missense_variant",
        "gene_symbol": "MYH7",
        "hgvsp": "ENSP00000347507.3:p.Arg502Trp",
    }
    out = asyncio.run(protvar.fetch_protvar(vep))

    assert set(out.keys()) == _FOUND_KEYS
    assert out["found"] is True
    assert out["am_pathogenicity"] == 0.912
    assert out["am_class"] == "likely_pathogenic"
    assert out["uniprot"] == "P12883"
    assert out["protein_position"] == 502
    assert out["ref_aa"] == "Arg"
    assert out["alt_aa"] == "Trp"
    assert out["conservation_score"] is None
    assert out["pocket_score"] is None
    assert "Domain" in out["feature_types"]
    assert out["colocated_total"] == 2
    assert {v["variant"] for v in out["colocated_variants"]} == {"R->Q", "R->W"}
    assert all(v["source_db"] == "UniProt" for v in out["colocated_variants"])


def _write_synthetic_alphamissense(tmp_path: Path) -> Path:
    """Write a tiny AlphaMissense-format TSV (10 columns), bgzip + tabix it,
    and return the .gz path — proving the parse/read roundtrip the deploy-time
    slice relies on, without the multi-GB real download."""
    import pysam

    raw = tmp_path / "am_fixture.tsv"
    rows = [
        "chr14\t23425678\tG\tA\thg38\tP12883\tENST00000355349.4\tR502W\t0.912\tlikely_pathogenic",
        "chr14\t23425679\tT\tC\thg38\tP12883\tENST00000355349.4\tR502Q\t0.150\tlikely_benign",
        "chr7\t100000\tA\tG\thg38\tQ00000\tENST00000000001.1\tA10V\t0.500\tambiguous",
    ]
    raw.write_text("\n".join(rows) + "\n")
    gz = pysam.tabix_index(
        str(raw), seq_col=0, start_col=1, end_col=1, force=True, zerobased=False,
    )
    return Path(gz)


def test_synthetic_alphamissense_roundtrip(monkeypatch, tmp_path):
    from backend.clients import alphamissense

    gz = _write_synthetic_alphamissense(tmp_path)
    assert gz.exists()
    assert (tmp_path / (gz.name + ".tbi")).exists()
    monkeypatch.setenv("ALPHAMISSENSE_PATH", str(gz))

    vep = {
        "ok": True,
        "seq_region_name": "14",
        "start": 23425678,
        "allele_string": "C/T",
        "strand": -1,
        "most_severe_consequence": "missense_variant",
        "gene_symbol": "MYH7",
        "hgvsp": "ENSP00000347507.3:p.Arg502Trp",
    }

    am = asyncio.run(alphamissense.fetch_alphamissense(vep))
    assert am["available"] is True
    assert am["score"] == 0.912
    assert am["classification"] == "likely_pathogenic"
    assert am["protein_variant"] == "R502W"

    _patch_async_client(monkeypatch, _no_network_transport())

    async def _no_uniprot(gene):
        return {"ok": True, "found": False, "gene": gene}

    monkeypatch.setattr(protvar, "fetch_uniprot", _no_uniprot)

    out = asyncio.run(protvar.fetch_protvar(vep))
    assert set(out.keys()) == _FOUND_KEYS
    assert out["found"] is True
    assert out["am_pathogenicity"] == 0.912
    assert out["am_class"] == "likely_pathogenic"
    assert out["uniprot"] is None
    assert out["protein_position"] == 502
    assert out["ref_aa"] == "Arg"
    assert out["alt_aa"] == "Trp"
    assert out["conservation_score"] is None
    assert out["pocket_score"] is None
    assert out["feature_types"] == []
    assert out["colocated_variants"] == []
    assert out["colocated_total"] == 0


_LIVE_MAPPING_PAYLOAD = {
    "content": {
        "inputs": [{
            "derivedGenomicVariants": [{
                "chromosome": "14", "position": 23425678,
                "refBase": "G", "altBase": "A",
                "genes": [{
                    "geneName": "MYH7",
                    "refAllele": "G", "altAllele": "A",
                    "isoforms": [{
                        "canonical": True,
                        "accession": "P12883",
                        "isoformPosition": 502,
                        "refAA": "Arg",
                        "variantAA": "Trp",
                        "consequences": "missense",
                        "amScore": {"amPathogenicity": 0.95, "amClass": "likely_pathogenic"},
                    }],
                }],
            }],
        }],
    },
}


def test_live_fallback_when_local_absent(monkeypatch):
    monkeypatch.delenv("ALPHAMISSENSE_PATH", raising=False)

    async def _no_uniprot(gene):
        return {"ok": True, "found": False, "gene": gene}

    monkeypatch.setattr(protvar, "fetch_uniprot", _no_uniprot)

    seen_urls: list[str] = []

    def _handler(request: httpx.Request) -> httpx.Response:
        seen_urls.append(str(request.url))
        path = request.url.path
        if path.endswith("/mapping"):
            return httpx.Response(200, json=_LIVE_MAPPING_PAYLOAD)
        if "/function/" in path:
            return httpx.Response(200, json={})
        if "/population/" in path:
            return httpx.Response(200, json={"variants": []})
        return httpx.Response(404, json={})  # pragma: no cover

    _patch_async_client(monkeypatch, httpx.MockTransport(_handler))

    ran = {"live": False}
    orig_live = protvar._fetch_protvar_live

    async def _spy(vep):
        ran["live"] = True
        return await orig_live(vep)

    monkeypatch.setattr(protvar, "_fetch_protvar_live", _spy)

    out = asyncio.run(protvar.fetch_protvar(_myh7_v39m_vep()))

    assert ran["live"] is True, "live fallback should run when both local sources absent"
    assert any("/mapping" in u for u in seen_urls)
    assert out["ok"] is True
    assert out["applicable"] is True
    assert out["found"] is True
    assert out["gene"] == "MYH7"
    assert out["uniprot"] == "P12883"
    assert out["protein_position"] == 502
    assert out["ref_aa"] == "Arg"
    assert out["alt_aa"] == "Trp"
    assert out["am_pathogenicity"] is None
    assert out["am_class"] is None
    assert set(out.keys()) == _FOUND_KEYS


def test_non_missense_short_circuits(monkeypatch):
    _patch_async_client(monkeypatch, _no_network_transport())
    out = asyncio.run(protvar.fetch_protvar({
        "ok": True,
        "most_severe_consequence": "stop_gained",
        "seq_region_name": "14", "start": 1, "allele_string": "C/T", "strand": 1,
    }))
    assert out == {
        "ok": True,
        "applicable": False,
        "reason": "non-missense variant",
        "consequence": "stop_gained",
    }


def test_vep_unavailable_returns_error():
    out = asyncio.run(protvar.fetch_protvar({"ok": False}))
    assert out["ok"] is False
    assert "VEP lookup unavailable" in out["error"]
