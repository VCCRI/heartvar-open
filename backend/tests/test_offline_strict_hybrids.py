"""HEARTVAR_OFFLINE_STRICT gate on the offline-primary hybrid clients.

Each test forces a LOCAL miss (so control reaches the gated live-fallback),
sets the flag, and monkeypatches the client's live network function(s) to raise.
A passing test proves the client returns a graceful shape WITHOUT any network
call in strict mode. Default-off behaviour is covered by the rest of the suite
(which runs with the flag unset).

``asyncio.run`` drives the coroutines so these pass under the repo's
``PYTEST_DISABLE_PLUGIN_AUTOLOAD=1`` convention (no pytest-asyncio).
"""

from __future__ import annotations

import asyncio
import importlib

import pytest


async def _boom_async(*a, **k):
    raise AssertionError("live network call made in offline-strict mode")


def _boom_sync(*a, **k):
    raise AssertionError("live network client opened in offline-strict mode")


@pytest.fixture
def strict(monkeypatch):
    monkeypatch.setenv("HEARTVAR_OFFLINE_STRICT", "1")
    return monkeypatch


@pytest.mark.parametrize(
    "modname,fetch,livefn",
    [
        ("uniprot", "fetch_uniprot", "_fetch_uniprot_live"),
        ("gtex", "fetch_gtex", "_fetch_gtex_live"),
        ("medgen", "fetch_medgen", "_fetch_medgen_live"),
        ("mgi", "fetch_mgi", "_fetch_mgi_live"),
    ],
)
def test_gene_keyed_client_strict_no_live(strict, modname, fetch, livefn):
    mod = importlib.import_module(f"backend.clients.{modname}")
    strict.setattr(mod, "_query_local_sync", lambda *a, **k: None, raising=True)
    strict.setattr(mod, livefn, _boom_async, raising=True)
    res = asyncio.run(getattr(mod, fetch)("MYH7"))
    assert isinstance(res, dict)


def test_conservation_strict_no_ucsc(strict):
    from backend.clients import conservation

    strict.delenv("HEARTVAR_PHYLOP_PATH", raising=False)
    res = asyncio.run(conservation.fetch_phylop100way("1", 100000))
    assert res is None


def test_spliceai_strict_no_live(strict):
    from backend.clients import spliceai

    strict.setattr(spliceai, "_fetch_spliceai_live", _boom_async, raising=True)
    res = asyncio.run(spliceai.fetch_spliceai("1-100000-A-G", "GRCh37"))
    assert isinstance(res, dict)


def test_gnomad_strict_no_live(strict):
    from backend.clients import gnomad

    strict.setattr(gnomad, "make_async_client", _boom_sync, raising=True)
    strict.setattr(gnomad, "request_with_retry", _boom_async, raising=True)
    res = asyncio.run(gnomad.fetch_gnomad("1-100000-A-G", "ZZZ0NOTAGENE"))
    assert isinstance(res, dict)


def test_gnomad_indel_rsid_strict_no_recoder(strict):
    """The in-panel indel miss path must not call the live Ensembl recoder in
    strict mode — it marks the frequency unresolved instead."""
    from backend.clients import ensembl_vep, gnomad

    strict.setattr(ensembl_vep, "fetch_variant_recoder_rsid", _boom_async, raising=True)
    res = asyncio.run(
        gnomad._resolve_indel_freq_by_rsid(
            {"variant_found": False}, "NM_000257.4:c.1del", False
        )
    )
    assert res.get("indel_unresolved") is True


def test_panelapp_strict_no_live(strict):
    from backend.clients import panelapp

    strict.setattr(panelapp, "_load_panelapp_snapshot", lambda: None, raising=True)
    strict.setattr(panelapp, "make_async_client", _boom_async, raising=False)
    strict.setattr(panelapp, "request_with_retry", _boom_async, raising=False)
    res = asyncio.run(panelapp.fetch_panelapp("MYH7", None))
    assert isinstance(res, dict)


def test_protvar_strict_no_live(strict):
    from backend.clients import protvar

    strict.setattr(protvar, "_fetch_protvar_live", _boom_async, raising=True)
    vep = {
        "gene_symbol": "ZZZ0NOTAGENE",
        "most_severe_consequence": "missense_variant",
        "amino_acids": "A/T",
        "protein_start": 10,
        "protein_end": 10,
    }
    res = asyncio.run(protvar.fetch_protvar(vep))
    assert isinstance(res, dict)
