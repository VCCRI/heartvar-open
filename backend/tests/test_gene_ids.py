"""Unit tests for backend.clients._gene_ids.

The load-bearing test here is
``test_resolves_with_no_hgnc_set_as_in_production``: the 17 MB
``hgnc_complete_set.txt`` is BOTH gitignored and listed in ``.dockerignore``, so
it does not exist in the deployed container. A resolver that only reads it works
on a dev machine and silently fails-open in production — which would leave Open
Targets "unavailable" for every RefSeq curation while every local check passed.
That test pins the shipped-map path so the trap cannot come back.

Run:
    .venv/bin/python -m pytest backend/tests/test_gene_ids.py -q
"""
from __future__ import annotations

from pathlib import Path

import pytest

import backend.clients._gene_ids as gi

_ELN = "ENSG00000049540"
_MYH7 = "ENSG00000092054"
_TTN = "ENSG00000155657"


@pytest.fixture(autouse=True)
def _reset_caches():
    """The module caches for the process lifetime and ``_LOADED`` is sticky, so
    every test must start from an unloaded state to exercise the real load path.
    """
    gi._SYMBOL_TO_ENSG = None
    gi._ENTREZ_TO_ENSG = None
    gi._LOADED = False
    gi._WARNED_MISSING = False
    yield
    gi._SYMBOL_TO_ENSG = None
    gi._ENTREZ_TO_ENSG = None
    gi._LOADED = False
    gi._WARNED_MISSING = False


def test_shipped_map_is_present_and_tracked():
    """The compact map must exist in the tree — it is the only gene-ID source
    that reaches production. Rebuild with scripts/build_gene_id_map.py."""
    assert gi.GENE_ID_MAP_PATH.exists(), (
        f"{gi.GENE_ID_MAP_PATH} missing — run scripts/build_gene_id_map.py"
    )


def test_resolves_with_no_hgnc_set_as_in_production(monkeypatch):
    """THE production scenario: no hgnc_complete_set.txt anywhere, because it is
    excluded from the runtime image. Resolution must still work from the shipped
    map alone, or the Open Targets fix is inert once deployed."""
    monkeypatch.setattr(gi, "HGNC_TSV_PATH", Path("/nonexistent/hgnc_complete_set.txt"))
    assert gi.resolve_ensembl_gene_id("2006", "ELN") == _ELN
    assert gi.resolve_ensembl_gene_id("2006", None) == _ELN
    assert gi.resolve_ensembl_gene_id(None, "ELN") == _ELN
    assert gi.resolve_ensembl_gene_id(None, "TTN") == _TTN
    assert gi.resolve_ensembl_gene_id("4625", "MYH7") == _MYH7


def test_entrez_wins_over_symbol_when_they_disagree(monkeypatch):
    """VEP's Entrez id describes the variant actually being curated, while the
    symbol is user-supplied and may be an alias, so Entrez is authoritative."""
    monkeypatch.setattr(gi, "HGNC_TSV_PATH", Path("/nonexistent/x.txt"))
    assert gi.resolve_ensembl_gene_id("7273", "ELN") == _TTN


def test_ensg_passes_through_with_version_stripped(monkeypatch):
    """An Ensembl id needs no lookup at all — and must not touch the map."""
    monkeypatch.setattr(gi, "GENE_ID_MAP_PATH", Path("/nonexistent/map.json.gz"))
    monkeypatch.setattr(gi, "HGNC_TSV_PATH", Path("/nonexistent/x.txt"))
    assert gi.resolve_ensembl_gene_id(_ELN, None) == _ELN
    assert gi.resolve_ensembl_gene_id(_ELN + ".16", None) == _ELN


def test_fails_open_when_no_source_exists(monkeypatch, caplog):
    """Neither source present → resolution disabled, warned once, never raising.
    A curation must not break because a reference file is missing."""
    monkeypatch.setattr(gi, "GENE_ID_MAP_PATH", Path("/nonexistent/map.json.gz"))
    monkeypatch.setattr(gi, "HGNC_TSV_PATH", Path("/nonexistent/hgnc.txt"))
    with caplog.at_level("WARNING"):
        assert gi.resolve_ensembl_gene_id("2006", "ELN") is None
    assert "resolution disabled" in caplog.text
    assert "build_gene_id_map.py" in caplog.text, "the warning must say how to fix it"


def test_unreadable_map_falls_back_to_hgnc_set(monkeypatch, tmp_path):
    """A corrupt shipped map must not disable resolution when the full HGNC set
    is available (dev machines, the builder image)."""
    bad = tmp_path / "gene_id_map.json.gz"
    bad.write_bytes(b"not gzip at all")
    monkeypatch.setattr(gi, "GENE_ID_MAP_PATH", bad)
    if not gi.HGNC_TSV_PATH.exists():
        pytest.skip("hgnc_complete_set.txt not present in this checkout")
    assert gi.resolve_ensembl_gene_id("2006", "ELN") == _ELN


def test_junk_never_guesses(monkeypatch):
    monkeypatch.setattr(gi, "HGNC_TSV_PATH", Path("/nonexistent/x.txt"))
    assert gi.resolve_ensembl_gene_id("NOT_AN_ENSG", None) is None
    assert gi.resolve_ensembl_gene_id(None, "NOT_A_REAL_GENE_XYZ") is None
    assert gi.resolve_ensembl_gene_id(None, None) is None
    assert gi.resolve_ensembl_gene_id("", "") is None
