"""Unit tests for backend.clients.hgnc_alias.canonicalise_gene_symbol.

Covers the LOCAL-FIRST resolver against a tiny temporary map (so the asserts
hold whether the shipped asset is the real HGNC download or a hand-written
placeholder), plus a fail-open check when the map is absent.

The resolver caches the loaded map at module level, so every test rebinds
``hgnc_alias.MAP_PATH`` to a controlled fixture and resets the module's
load-state caches before exercising the function.

Self-contained: no network is touched.
"""
from __future__ import annotations

import json

import pytest

import backend.clients.hgnc_alias as hgnc_alias

_FIXTURE = {
    "approved": ["MYH7", "MYBPC3", "NKX2-5", "KCNQ1"],
    "aliases": {
        "CSX": "NKX2-5",
        "CMH1": "MYH7",
        "KCNA9": "KCNQ1",
    },
}

_RETURN_KEYS = {"input", "approved", "is_alias", "recognized", "ambiguous"}


def _reset_module_state() -> None:
    """Clear the module-level lazy caches so the next call reloads MAP_PATH."""
    hgnc_alias._APPROVED_UPPER_TO_CANONICAL = None
    hgnc_alias._ALIASES = None
    hgnc_alias._LOADED = False
    hgnc_alias._WARNED_MISSING = False


@pytest.fixture()
def map_path(tmp_path, monkeypatch):
    """Write the fixture map to a temp file and point the module at it."""
    p = tmp_path / "hgnc_alias_map.json"
    p.write_text(json.dumps(_FIXTURE), encoding="utf-8")
    monkeypatch.setattr(hgnc_alias, "MAP_PATH", p)
    _reset_module_state()
    yield p
    _reset_module_state()


def test_approved_symbol(map_path):
    out = hgnc_alias.canonicalise_gene_symbol("MYH7")
    assert set(out) == _RETURN_KEYS
    assert out["recognized"] is True
    assert out["is_alias"] is False
    assert out["approved"] == "MYH7"
    assert out["ambiguous"] is False
    assert out["input"] == "MYH7"


def test_known_alias(map_path):
    out = hgnc_alias.canonicalise_gene_symbol("CSX")
    assert out["recognized"] is True
    assert out["is_alias"] is True
    assert out["approved"] == "NKX2-5"


def test_unknown_symbol(map_path):
    out = hgnc_alias.canonicalise_gene_symbol("NOT_A_GENE_XYZ")
    assert out["recognized"] is False
    assert out["is_alias"] is False
    assert out["approved"] is None


def test_case_insensitive_approved(map_path):
    out = hgnc_alias.canonicalise_gene_symbol("nkx2-5")
    assert out["recognized"] is True
    assert out["is_alias"] is False
    assert out["approved"] == "NKX2-5"


def test_case_insensitive_alias(map_path):
    out = hgnc_alias.canonicalise_gene_symbol("csx")
    assert out["recognized"] is True
    assert out["is_alias"] is True
    assert out["approved"] == "NKX2-5"


def test_fail_open_when_map_absent(tmp_path, monkeypatch):
    """A missing map must not block a curation: pass the input through
    unchanged, flagged recognized=True so callers proceed."""
    missing = tmp_path / "does_not_exist.json"
    monkeypatch.setattr(hgnc_alias, "MAP_PATH", missing)
    _reset_module_state()
    try:
        out = hgnc_alias.canonicalise_gene_symbol("MYH7")
        assert out["recognized"] is True
        assert out["is_alias"] is False
        assert out["approved"] == "MYH7"
    finally:
        _reset_module_state()


def test_shipped_asset_known_cardiac_aliases():
    """If the shipped asset is present (real download or placeholder), the
    well-known cardiac aliases used in this repo must resolve. Skips when the
    asset hasn't been built so the suite passes in a clean checkout."""
    _reset_module_state()
    try:
        if not hgnc_alias.MAP_PATH.exists():
            pytest.skip("hgnc_alias_map.json not built; run build_hgnc_alias_db.py")
        out = hgnc_alias.canonicalise_gene_symbol("CSX")
        assert out["approved"] == "NKX2-5"
        assert out["is_alias"] is True
        approved = hgnc_alias.canonicalise_gene_symbol("MYH7")
        assert approved["approved"] == "MYH7"
        assert approved["is_alias"] is False
    finally:
        _reset_module_state()
