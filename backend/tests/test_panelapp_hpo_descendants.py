"""Unit tests for PanelApp's offline cardiac HPO descendant closure.

The closure (backend/data/cardiac_hpo_descendants.json, built by
scripts/build_hpo_labels.py from hp.obo) replaces the 5 sequential live JAX
/descendants calls that caused the ~9 s post-restart cold-start. These tests
confirm the local closure loads, drives check_hpo_relevance, and that the
live JAX path is NOT touched when the closure is present.

Runnable with pytest or directly (python -m backend.tests.test_panelapp_hpo_descendants).
"""
from __future__ import annotations

import asyncio
from pathlib import Path

import pytest

import backend.clients.panelapp as panelapp

_CLOSURE = Path(__file__).resolve().parent.parent / "data" / "cardiac_hpo_descendants.json"


def _reset():
    panelapp.CARDIAC_HPO_DESCENDANTS.clear()


@pytest.mark.skipif(not _CLOSURE.is_file(), reason="closure JSON not built")
def test_local_closure_loads_and_categorises():
    _reset()
    assert panelapp._load_local_descendants() is True
    assert len(panelapp.CARDIAC_HPO_DESCENDANTS) > 800
    hcm = panelapp.CARDIAC_HPO_DESCENDANTS.get("HP:0001639", set())
    assert "hcm" in hcm
    assert "Congenital heart disease" in hcm


@pytest.mark.skipif(not _CLOSURE.is_file(), reason="closure JSON not built")
def test_ensure_cache_uses_local_not_jax(monkeypatch):
    """When the closure is present, no live JAX /descendants call is made."""
    _reset()

    async def _boom(_hpo_id):  # pragma: no cover - must never run
        raise AssertionError("live JAX must not be called when closure exists")

    monkeypatch.setattr(panelapp, "fetch_hpo_descendants", _boom)
    asyncio.run(panelapp._ensure_descendant_cache())
    assert len(panelapp.CARDIAC_HPO_DESCENDANTS) > 800
    assert panelapp.check_hpo_relevance("HP:0001639", "hcm") is True


def test_ensure_cache_falls_back_to_jax_when_closure_absent(monkeypatch):
    """With no closure file, the cache is built from the live JAX fetch."""
    _reset()
    monkeypatch.setattr(panelapp, "_DESCENDANTS_JSON_PATH", Path("/nonexistent/closure.json"))

    calls: list[str] = []

    async def _fake(hpo_id):
        calls.append(hpo_id)
        return {"HP:0001639"}

    monkeypatch.setattr(panelapp, "fetch_hpo_descendants", _fake)
    asyncio.run(panelapp._ensure_descendant_cache())
    assert calls, "live JAX fetch should have been used"
    assert "HP:0001639" in panelapp.CARDIAC_HPO_DESCENDANTS


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
