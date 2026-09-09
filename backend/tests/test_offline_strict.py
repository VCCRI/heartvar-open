"""HEARTVAR_OFFLINE_STRICT: offline-primary clients must not make a live call
on a local miss when the flag is set. Default off keeps the live fallback.

Async is driven with ``asyncio.run`` inside sync tests so these pass under the
repo's ``PYTEST_DISABLE_PLUGIN_AUTOLOAD=1`` convention (no pytest-asyncio).
"""

from __future__ import annotations

import asyncio

import pytest

from backend.clients import _offline


def test_offline_strict_flag_parsing(monkeypatch):
    monkeypatch.delenv("HEARTVAR_OFFLINE_STRICT", raising=False)
    assert _offline.offline_strict() is False
    for truthy in ("1", "true", "yes", "on", "ON", "True", " yes "):
        monkeypatch.setenv("HEARTVAR_OFFLINE_STRICT", truthy)
        assert _offline.offline_strict() is True, truthy
    for falsy in ("0", "", "no", "off", "false", "disabled"):
        monkeypatch.setenv("HEARTVAR_OFFLINE_STRICT", falsy)
        assert _offline.offline_strict() is False, falsy


def test_opentargets_strict_skips_live(monkeypatch):
    """In strict mode Open Targets must return None (unavailable) without ever
    opening a live client."""
    from backend.clients import opentargets_client as ot

    monkeypatch.setenv("HEARTVAR_OFFLINE_STRICT", "1")
    monkeypatch.setenv("OPENTARGETS_DB_PATH", "/nonexistent/heartvar-test-ot.db")

    def _boom(*a, **k):
        raise AssertionError("make_async_client called in offline-strict mode")

    monkeypatch.setattr(ot, "make_async_client", _boom)
    result = asyncio.run(
        ot.OpenTargetsClient().get_gene_disease_evidence("ENSG00000106804")
    )
    assert result is None


def test_opentargets_reaches_live_when_not_strict(monkeypatch):
    """Sanity: with the flag off, the code path DOES reach make_async_client
    (proven by a sentinel — no real network is performed)."""
    from backend.clients import opentargets_client as ot

    monkeypatch.delenv("HEARTVAR_OFFLINE_STRICT", raising=False)
    monkeypatch.setenv("OPENTARGETS_DB_PATH", "/nonexistent/heartvar-test-ot.db")

    class _Sentinel(Exception):
        pass

    def _boom(*a, **k):
        raise _Sentinel()

    monkeypatch.setattr(ot, "make_async_client", _boom)
    with pytest.raises(_Sentinel):
        asyncio.run(
            ot.OpenTargetsClient().get_gene_disease_evidence("ENSG00000106804")
        )


def test_gencc_strict_no_live_refresh(monkeypatch):
    """In strict mode GenCC serves the local snapshot (or a clean 'unavailable')
    and never calls the live thegencc.org refresh."""
    from backend.clients import gencc

    monkeypatch.setenv("HEARTVAR_OFFLINE_STRICT", "1")
    monkeypatch.setattr(gencc, "_INDEX", None, raising=False)
    monkeypatch.setattr(gencc, "_META", None, raising=False)

    async def _boom():
        raise AssertionError(
            "_refresh_cache (live thegencc.org) called in offline-strict mode"
        )

    monkeypatch.setattr(gencc, "_refresh_cache", _boom)

    res = asyncio.run(gencc.fetch_gencc("MYH7"))
    assert isinstance(res, dict) and "ok" in res


def test_opentargets_local_snapshot_used(monkeypatch, tmp_path):
    """With a local snapshot present, Open Targets reads it and makes NO live
    call — the offline builder's output path (data/opentargets.db)."""
    import json
    import sqlite3

    from backend.clients import opentargets_client as ot

    db = tmp_path / "opentargets.db"
    con = sqlite3.connect(db)
    con.execute(
        "CREATE TABLE opentargets (ensembl_id TEXT PRIMARY KEY, symbol TEXT, payload TEXT)"
    )
    target = {
        "approvedSymbol": "MYH7",
        "associatedDiseases": {
            "rows": [
                {
                    "disease": {"id": "EFO_0000407", "name": "dilated cardiomyopathy"},
                    "score": 0.88,
                    "datatypeScores": [{"id": "genetic_association", "score": 0.9}],
                }
            ]
        },
    }
    con.execute(
        "INSERT INTO opentargets VALUES (?, ?, ?)",
        ("ENSG00000092054", "MYH7", json.dumps(target)),
    )
    con.commit()
    con.close()

    monkeypatch.setenv("OPENTARGETS_DB_PATH", str(db))
    monkeypatch.delenv("HEARTVAR_OFFLINE_STRICT", raising=False)

    def _boom(*a, **k):
        raise AssertionError("live GraphQL called despite a local snapshot")

    monkeypatch.setattr(ot, "make_async_client", _boom)
    res = asyncio.run(
        ot.OpenTargetsClient().get_gene_disease_evidence("ENSG00000092054")
    )
    assert isinstance(res, dict)
    assert res.get("ensembl_id") == "ENSG00000092054"
    names = [str(d.get("name", "")).lower() for d in (res.get("top_diseases") or [])]
    assert any("cardiomyopathy" in n for n in names)
