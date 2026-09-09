"""Endpoint test for GET /api/gene/validate — the landing form's warn-early
HGNC gene-symbol check.

The endpoint maps ``canonicalise_gene_symbol()``'s result to a
``{symbol, status, approved, message}`` shape the frontend renders. The
canonicaliser is monkeypatched to exercise each branch deterministically (the
local HGNC map's contents vary across builds), plus one real-symbol smoke test
that hits the actual map to confirm the wiring.

No network: a plain ``TestClient(app_module.app)`` (no ``with``) skips the lifespan
warm-up, and the endpoint only does a local dict lookup.

Runnable with pytest or directly
(``python -m backend.tests.test_gene_validate_endpoint``).
"""
from __future__ import annotations

import types

import pytest
from starlette.testclient import TestClient

import backend.app as app_module


@pytest.fixture(autouse=True)
def _restore_canon():
    """Each test may monkeypatch the module-level canonicaliser; snapshot +
    restore it around every test so the fake never leaks."""
    original = app_module.canonicalise_gene_symbol
    yield
    app_module.canonicalise_gene_symbol = original


def _get(symbol: str) -> dict:
    client = TestClient(app_module.app)
    resp = client.get("/api/gene/validate", params={"symbol": symbol})
    assert resp.status_code == 200, resp.text
    return resp.json()


def test_gene_validate_empty():
    d = _get("")
    assert d["status"] == "empty"
    assert d["approved"] is None
    assert d["message"] == ""


def test_gene_validate_approved():
    app_module.canonicalise_gene_symbol = lambda g: {
        "input": g, "approved": "MYH7", "is_alias": False,
        "recognized": True, "ambiguous": False,
    }
    d = _get("MYH7")
    assert d["status"] == "approved"
    assert d["approved"] == "MYH7"
    assert d["message"] == ""


def test_gene_validate_alias_suggests_approved():
    app_module.canonicalise_gene_symbol = lambda g: {
        "input": g, "approved": "MYH7", "is_alias": True,
        "recognized": True, "ambiguous": False,
    }
    d = _get("CMH1")
    assert d["status"] == "alias"
    assert d["approved"] == "MYH7"
    assert "MYH7" in d["message"]


def test_gene_validate_unrecognized():
    app_module.canonicalise_gene_symbol = lambda g: {
        "input": g, "approved": None, "is_alias": False,
        "recognized": False, "ambiguous": False,
    }
    d = _get("NOTAREALGENE123")
    assert d["status"] == "unrecognized"
    assert d["approved"] is None
    assert "not a recognised" in d["message"].lower()


def test_gene_validate_real_symbol_smoke():
    d = _get("MYH7")
    assert d["status"] in ("approved", "alias")
    assert d["approved"]


def _run_all():
    """Manual runner (pytest fixtures don't fire under direct execution), so we
    restore the canonicaliser between tests ourselves."""
    original = app_module.canonicalise_gene_symbol
    fns = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and isinstance(v, types.FunctionType)]
    try:
        for fn in fns:
            fn()
            app_module.canonicalise_gene_symbol = original
            print(f"  ok  {fn.__name__}")
    finally:
        app_module.canonicalise_gene_symbol = original
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
