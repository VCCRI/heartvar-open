"""Audit §7a/§6b: a VEP failure must ABORT the curation with a clear, kind-aware
error — never silently degrade to a 'VUS'. The dangerous path this guards: a
VEP blip → variant_id=None → fetch_gnomad(None) returns variant_found:False →
PM2 fires 'absent' on a transient network error. Aborting before the dependent
spawn removes that path entirely.

Exercises the REAL _gather_evidence_sse with every network fetcher stubbed."""

from __future__ import annotations

import asyncio
import json

import pytest

import backend.app as app
import backend.evidence as evidence


def _ok_stub(**extra):
    async def _f(*a, **k):
        return {"ok": True, "found": False, **extra}
    return _f


@pytest.fixture
def stub_fetchers(monkeypatch):
    """Replace every network fetcher reached by the gather with a fast stub so
    the test is offline and deterministic. The VEP stub is set per-test.

    The fetchers + gather machinery were extracted into backend.evidence, so
    the stubs are installed there (where _gather_evidence_sse now resolves its
    module globals) — not on backend.app, which only re-exports them."""
    names = [
        "fetch_clinvar", "fetch_uniprot", "fetch_gtex", "fetch_fetal_heart",
        "fetch_panelapp", "fetch_gencc", "fetch_medgen", "fetch_mgi",
        "fetch_biogrid", "fetch_gnomad", "fetch_spliceai",
        "fetch_protvar", "fetch_alphamissense", "fetch_opentargets",
        "fetch_erepo", "get_pm5_evidence", "get_gene_variant_landscape",
        "_build_domain_plp_evidence", "fetch_pmc_fulltext",
    ]
    for n in names:
        monkeypatch.setattr(evidence, n, _ok_stub(), raising=False)

    async def _lit(*a, **k):
        return {"ok": True, "variant_papers": [], "papers": []}
    monkeypatch.setattr(evidence, "_fetch_literature_with_vep_aa", _lit, raising=False)
    monkeypatch.setattr(evidence, "get_gene_disease_literature", _lit, raising=False)

    async def _strings(*a, **k):
        return []
    monkeypatch.setattr(evidence, "get_gene_phenotype_strings", _strings, raising=False)
    monkeypatch.setattr(evidence, "build_cardiac_keywords", lambda *a, **k: [], raising=False)
    monkeypatch.setattr(evidence, "chdgene_lookup", lambda *a, **k: {"ok": True, "found": False}, raising=False)
    monkeypatch.setattr(evidence, "resolved_token_map", lambda *a, **k: {}, raising=False)


def _gather(gene="MYH7", hgvs_c="c.1988G>A"):
    req = app.CurationRequest(gene=gene, hgvs_c=hgvs_c, ai_mode="none", enable_erepo=False)
    state: dict = {}
    events = []

    async def _run():
        async for ev in app._gather_evidence_sse(req, state):
            events.append(ev)
    asyncio.run(_run())
    return state, _parse(events)


def _parse(events):
    out = []
    for block in events:
        ev, data = None, None
        for line in block.splitlines():
            if line.startswith("event: "):
                ev = line[len("event: "):].strip()
            elif line.startswith("data: "):
                data = json.loads(line[len("data: "):])
        if ev:
            out.append((ev, data or {}))
    return out


def _vep(kind, err):
    async def _f(*a, **k):
        return {"ok": False, "error": err, "failure_kind": kind}
    return _f


def test_vep_not_found_aborts_with_clear_error(stub_fetchers, monkeypatch):
    monkeypatch.setattr(evidence, "fetch_vep", _vep("not_found", "400: ref allele mismatch"))
    state, events = _gather()
    kinds = [d.get("kind") for e, d in events if e == "error"]
    assert "vep_unresolved" in kinds, events
    assert "evidence" not in state
    assert not any(e in ("db_only_complete", "stage2_complete") for e, _ in events)


def test_vep_reference_allele_mismatch_blames_gene(stub_fetchers, monkeypatch):
    monkeypatch.setattr(evidence, "fetch_vep", _vep(
        "not_found",
        "400: {\"error\":\"Unable to parse HGVS notation "
        "'ENST00000423902:c.1988G>A': : Reference allele extracted from "
        "ENST00000423902:60781322-60781322 (A) does not match reference allele "
        "given by HGVS notation ENST00000423902:c.1988G>A\"}",
    ))
    state, events = _gather()
    err = [d for e, d in events if e == "error"]
    assert err and err[0].get("kind") == "vep_gene_mismatch", events
    msg = err[0].get("message", "")
    assert "gene" in msg.lower(), msg
    assert "400:" not in msg and "HGVS notation" not in msg, msg
    assert "evidence" not in state


def test_vep_transient_aborts_with_retry_message(stub_fetchers, monkeypatch):
    monkeypatch.setattr(evidence, "fetch_vep", _vep("transient", "transient Ensembl failure after retries"))
    state, events = _gather()
    err = [d for e, d in events if e == "error"]
    assert err and err[0].get("kind") == "vep_transient", events
    assert "retry" in err[0].get("message", "").lower()
    assert "evidence" not in state


def test_vep_missing_failure_kind_defaults_to_unresolved(stub_fetchers, monkeypatch):
    async def _vep_bare(*a, **k):
        return {"ok": False, "error": "VEP lookup failed"}
    monkeypatch.setattr(evidence, "fetch_vep", _vep_bare)
    state, events = _gather()
    kinds = [d.get("kind") for e, d in events if e == "error"]
    assert kinds == ["vep_unresolved"], events


def test_vep_failure_does_not_spawn_gnomad(stub_fetchers, monkeypatch):
    """The whole point: dependent variant-level sources must NOT be queried with
    variant_id=None (which would make gnomAD read 'absent' and fire PM2)."""
    called = {"gnomad": False}

    async def _gnomad_spy(*a, **k):
        called["gnomad"] = True
        return {"ok": True, "variant_found": False}
    monkeypatch.setattr(evidence, "fetch_gnomad", _gnomad_spy)
    monkeypatch.setattr(evidence, "fetch_vep", _vep("not_found", "400"))
    _gather()
    assert called["gnomad"] is False


def _vep_ok():
    async def _f(*a, **k):
        return {
            "ok": True, "assembly_name": "GRCh38", "seq_region_name": "14",
            "start": 23429116, "end": 23429116, "allele_string": "C/T", "strand": -1,
            "hgvsc": "c.1988G>A", "hgvsp": "p.Arg663His",
            "gene_id": "ENSG00000092054", "transcript_id": "ENST00000355349",
            "most_severe_consequence": "missense_variant",
            "protein_start": 663, "protein_end": 663,
        }
    return _f


def test_deadline_cuts_slow_source_but_completes(stub_fetchers, monkeypatch):
    """A slow non-critical source is cut at the deadline, but because VEP
    resolved the curation still completes (partial result, not an abort)."""
    monkeypatch.setattr(evidence, "fetch_vep", _vep_ok())

    async def _slow(*a, **k):
        await asyncio.sleep(5)
        return {"ok": True}
    monkeypatch.setattr(evidence, "fetch_gtex", _slow)
    monkeypatch.setattr(evidence, "_CURATION_DEADLINE_S", 0.5)
    state, events = _gather()
    deadline_marked = [
        (e, d) for e, d in events
        if e == "db_done" and (d.get("data") or {}).get("failure_kind") == "deadline"
    ]
    assert deadline_marked, events
    assert "evidence" in state
    assert not any(e == "error" for e, _ in events)


def test_deadline_cutting_gnomad_aborts(stub_fetchers, monkeypatch):
    """If a TIER-CRITICAL source (gnomAD) is cut at the deadline, abort + retry
    rather than score on incomplete frequency evidence (review finding)."""
    monkeypatch.setattr(evidence, "fetch_vep", _vep_ok())

    async def _slow(*a, **k):
        await asyncio.sleep(5)
        return {"ok": True, "variant_found": False}
    monkeypatch.setattr(evidence, "fetch_gnomad", _slow)
    monkeypatch.setattr(evidence, "_CURATION_DEADLINE_S", 0.5)
    state, events = _gather()
    err = [d for e, d in events if e == "error"]
    assert err and err[0].get("kind") == "vep_transient", events
    assert "gnomad" in err[0].get("message", "").lower()
    assert "evidence" not in state


def test_deadline_with_vep_hang_aborts(stub_fetchers, monkeypatch):
    """If VEP itself hasn't resolved by the deadline, abort like a transient
    VEP failure rather than score a truncated result."""
    async def _slow_vep(*a, **k):
        await asyncio.sleep(5)
        return {"ok": True}
    monkeypatch.setattr(evidence, "fetch_vep", _slow_vep)
    monkeypatch.setattr(evidence, "_CURATION_DEADLINE_S", 0.3)
    state, events = _gather()
    err = [d for e, d in events if e == "error"]
    assert err and err[0].get("kind") == "vep_transient", events
    assert "evidence" not in state


def test_protein_input_aborts_with_clear_message(stub_fetchers):
    state, events = _gather(hgvs_c="p.Arg403Gln")
    err = [d for e, d in events if e == "error"]
    assert err and err[0].get("kind") == "input_protein", events
    assert "evidence" not in state


def test_rsid_unresolved_aborts(stub_fetchers, monkeypatch):
    async def _no_recode(*a, **k):
        return None
    monkeypatch.setattr(evidence, "recode_rsid_to_coords", _no_recode)
    state, events = _gather(hgvs_c="rs999999999")
    err = [d for e, d in events if e == "error"]
    assert err and err[0].get("kind") == "rsid_unresolved", events
    assert "evidence" not in state


def test_gene_alias_canonicalised_in_gather(stub_fetchers, monkeypatch):
    monkeypatch.setattr(evidence, "canonicalise_gene_symbol", lambda g: {
        "input": g, "approved": "MYH7", "is_alias": True,
        "recognized": True, "ambiguous": False,
    })
    monkeypatch.setattr(evidence, "fetch_vep", _vep_ok())
    state, _events = _gather(gene="CSX", hgvs_c="c.1208G>A")
    assert state.get("gene_canonicalised") == {"from": "CSX", "to": "MYH7"}
    assert state.get("gene") == "MYH7"
