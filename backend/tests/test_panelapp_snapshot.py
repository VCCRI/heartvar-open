"""Unit tests for PanelApp's local snapshot lookup + live fallback.

fetch_panelapp prefers a local {GENE: [record,…]} snapshot
(scripts/build_panelapp_snapshot.py) and falls back to the live
panelapp-aus.org API when no snapshot is built. These tests prove (1) a
snapshot hit produces the canonical result shape with NO network, and (2)
the live fallback still works when the snapshot is absent.

Runnable with pytest or directly
(python -m backend.tests.test_panelapp_snapshot).
"""
from __future__ import annotations

import asyncio
import json

import httpx
import pytest

import backend.clients.panelapp as panelapp

_REAL_ASYNC_CLIENT = httpx.AsyncClient


def _reset_snapshot(monkeypatch, path: str | None):
    panelapp._panelapp_snapshot = None
    panelapp._panelapp_snapshot_loaded = False
    panelapp.CARDIAC_HPO_DESCENDANTS.clear()
    if path is None:
        monkeypatch.delenv("PANELAPP_SNAPSHOT_PATH", raising=False)
    else:
        monkeypatch.setenv("PANELAPP_SNAPSHOT_PATH", path)


def test_snapshot_hit_no_network(monkeypatch, tmp_path):
    snap = tmp_path / "panelapp_aus_snapshot.json"
    snap.write_text(json.dumps({"genes": {"MYH7": [{
        "panel": {"id": 1, "name": "Hypertrophic cardiomyopathy", "version": "1.5"},
        "confidence_level": "3",
        "mode_of_inheritance": "BIALLELIC, autosomal or pseudoautosomal",
        "phenotypes": ["Hypertrophic cardiomyopathy"],
        "entity_status": "green",
    }]}}))
    _reset_snapshot(monkeypatch, str(snap))

    def _boom(*a, **k):  # pragma: no cover - must never run
        raise AssertionError("snapshot hit must not touch the network")

    monkeypatch.setattr(panelapp.httpx, "AsyncClient", _boom)
    res = asyncio.run(panelapp.fetch_panelapp("MYH7", "HP:0001639"))
    assert res["ok"] is True
    assert res["on_cardiovascular_panel"] is True
    assert res["green_panels"] == 1
    assert res["panels_found"][0]["panel_name"] == "Hypertrophic cardiomyopathy"
    assert res["panels_found"][0]["confidence"] == "green"
    assert res["panels_found"][0]["contributes_to_pp4"] is True


def test_snapshot_gene_absent_returns_empty_no_network(monkeypatch, tmp_path):
    snap = tmp_path / "panelapp_aus_snapshot.json"
    snap.write_text(json.dumps({"genes": {"MYH7": []}}))
    _reset_snapshot(monkeypatch, str(snap))
    monkeypatch.setattr(panelapp.httpx, "AsyncClient",
                        lambda *a, **k: (_ for _ in ()).throw(
                            AssertionError("must not touch network")))
    res = asyncio.run(panelapp.fetch_panelapp("BRCA1", "HP:0001639"))
    assert res["ok"] is True
    assert res["on_cardiovascular_panel"] is False
    assert res["panels_found"] == []


def test_live_fallback_when_snapshot_absent(monkeypatch, tmp_path):
    _reset_snapshot(monkeypatch, str(tmp_path / "nope.json"))

    def handler(request: httpx.Request) -> httpx.Response:
        assert "/genes/" in str(request.url)
        return httpx.Response(200, json={"results": [{
            "panel": {"id": 2, "name": "Dilated Cardiomyopathy", "version": "2.0"},
            "confidence_level": "3",
            "mode_of_inheritance": "MONOALLELIC",
            "phenotypes": ["Dilated cardiomyopathy"],
            "entity_status": "green",
        }]})

    def _mock_client(*a, **k):
        return _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler),
                                  follow_redirects=True)

    monkeypatch.setattr(panelapp.httpx, "AsyncClient", _mock_client)
    res = asyncio.run(panelapp.fetch_panelapp("LMNA", "HP:0001644"))
    assert res["ok"] is True
    assert res["on_cardiovascular_panel"] is True
    assert res["panels_found"][0]["panel_name"] == "Dilated Cardiomyopathy"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
