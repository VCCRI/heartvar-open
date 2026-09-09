"""Unit tests for the phyloP100way conservation backfill (Rec B).

Covers the contig-name mapping, the None guards, the BP7-eligibility gate, and
the VEP-result enrichment that backfills phyloP for synonymous / non-coding
variants where dbNSFP returned nothing. Fully offline — the actual bigWig read
(``fetch_phylop100way``) is stubbed; no pyBigWig / network needed.

Runnable directly (``python -m backend.tests.test_conservation``).
"""
from __future__ import annotations

import asyncio

import backend.clients.ensembl_vep as vep
from backend.clients.conservation import _ucsc_chrom, fetch_phylop100way
from backend.clients.ensembl_vep import _needs_phylop_backfill


def test_ucsc_chrom_mapping():
    assert _ucsc_chrom("14") == "chr14"
    assert _ucsc_chrom("X") == "chrX"
    assert _ucsc_chrom("MT") == "chrM"
    assert _ucsc_chrom("M") == "chrM"
    assert _ucsc_chrom("chr7") == "chr7"


def test_fetch_phylop_none_guards():
    assert asyncio.run(fetch_phylop100way(None, None)) is None
    assert asyncio.run(fetch_phylop100way("14", None)) is None
    assert asyncio.run(fetch_phylop100way(None, 100)) is None


def test_needs_phylop_backfill_gate():
    assert _needs_phylop_backfill("synonymous_variant") is True
    assert _needs_phylop_backfill("intron_variant") is True
    assert _needs_phylop_backfill("splice_region_variant&intron_variant") is True
    assert _needs_phylop_backfill("3_prime_utr_variant") is True
    assert _needs_phylop_backfill("missense_variant") is False
    assert _needs_phylop_backfill("stop_gained") is False
    assert _needs_phylop_backfill("") is False


def _with_stub(return_value, fn):
    """Swap ensembl_vep.fetch_phylop100way for a stub, run fn(), restore."""
    async def _stub(chrom, pos):
        return return_value
    orig = vep.fetch_phylop100way
    vep.fetch_phylop100way = _stub
    try:
        return fn()
    finally:
        vep.fetch_phylop100way = orig


def test_enrich_backfills_synonymous_when_phylop_absent():
    r = _with_stub(1.5, lambda: asyncio.run(vep._enrich_phylop_conservation({
        "ok": True, "most_severe_consequence": "synonymous_variant",
        "seq_region_name": "11", "start": 100, "phylop100way": None,
        "allele_string": "C/T",
    })))
    assert r["phylop100way"] == 1.5


def test_enrich_skips_non_snv_indel():
    r = _with_stub(1.5, lambda: asyncio.run(vep._enrich_phylop_conservation({
        "ok": True, "most_severe_consequence": "intron_variant",
        "seq_region_name": "11", "start": 100, "phylop100way": None,
        "allele_string": "C/-",
    })))
    assert r.get("phylop100way") is None


def test_enrich_skips_missense():
    r = _with_stub(1.5, lambda: asyncio.run(vep._enrich_phylop_conservation({
        "ok": True, "most_severe_consequence": "missense_variant",
        "seq_region_name": "11", "start": 100, "phylop100way": None,
    })))
    assert r.get("phylop100way") is None


def test_enrich_skips_when_phylop_already_present():
    r = _with_stub(9.9, lambda: asyncio.run(vep._enrich_phylop_conservation({
        "ok": True, "most_severe_consequence": "synonymous_variant",
        "seq_region_name": "11", "start": 100, "phylop100way": 0.3,
    })))
    assert r["phylop100way"] == 0.3


def test_enrich_noop_on_failed_vep():
    r = _with_stub(1.5, lambda: asyncio.run(vep._enrich_phylop_conservation({
        "ok": False,
    })))
    assert r.get("phylop100way") is None


def test_enrich_handles_lookup_miss():
    r = _with_stub(None, lambda: asyncio.run(vep._enrich_phylop_conservation({
        "ok": True, "most_severe_consequence": "synonymous_variant",
        "seq_region_name": "11", "start": 100, "phylop100way": None,
        "allele_string": "C/T",
    })))
    assert r.get("phylop100way") is None


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"ok  {fn.__name__}")
    print(f"\n{len(fns)} tests passed")
