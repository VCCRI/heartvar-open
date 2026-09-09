"""Unit tests for backend.clients.erepo_client.

Covers the local-dump-first lookup (against the real on-disk
``backend/data/erepo_all.tsv`` when present), transcript-agnostic matching,
a clean miss that does NOT touch the network, and the live-API fallback when
the dump is absent.

Self-contained: the live path is exercised with an httpx.MockTransport so no
network is touched. Runnable with pytest
(``pytest backend/tests/test_erepo_client.py``) or directly
(``python -m backend.tests.test_erepo_client``).
"""
from __future__ import annotations

import asyncio
from pathlib import Path

import httpx
import pytest

import backend.clients.erepo_client as erepo

_REAL_ASYNC_CLIENT = httpx.AsyncClient
_TSV = Path(__file__).resolve().parent.parent / "data" / "erepo_all.tsv"


def _reset_state(monkeypatch, tsv_path: str | None):
    """Clear the module's lazy caches and point the TSV path somewhere."""
    erepo._tsv_index = None
    erepo._tsv_index_loaded = False
    erepo._cache.clear()
    if tsv_path is None:
        monkeypatch.delenv("HEARTVAR_EREPO_TSV", raising=False)
    else:
        monkeypatch.setenv("HEARTVAR_EREPO_TSV", tsv_path)


@pytest.mark.skipif(not _TSV.is_file(), reason="erepo_all.tsv not built locally")
def test_local_dump_hit(monkeypatch):
    """A known PAH record resolves from the local dump with full detail."""
    _reset_state(monkeypatch, str(_TSV))
    res = asyncio.run(erepo.fetch_erepo("PAH", "NM_000277.2:c.1A>G"))
    assert res["ok"] is True
    assert res["found"] is True
    assert res["classification"] == "Pathogenic"
    assert "PS3" in res["criteria"] and "PM2" in res["criteria"]
    assert res["vcep"] == "Phenylketonuria VCEP"
    assert res["url"].startswith("http")


@pytest.mark.skipif(not _TSV.is_file(), reason="erepo_all.tsv not built locally")
def test_local_dump_transcript_agnostic(monkeypatch):
    """Matching ignores the transcript prefix (bare c. form)."""
    _reset_state(monkeypatch, str(_TSV))
    res = asyncio.run(erepo.fetch_erepo("PAH", "NM_999999.9(PAH):c.1A>G (p.Met1Val)"))
    assert res["found"] is True
    assert res["classification"] == "Pathogenic"


@pytest.mark.skipif(not _TSV.is_file(), reason="erepo_all.tsv not built locally")
def test_local_dump_miss_no_network(monkeypatch):
    """A variant absent from the dump returns found=False WITHOUT any HTTP.

    We sabotage the live client so any network attempt would raise — proving
    the local path is authoritative on a miss."""
    _reset_state(monkeypatch, str(_TSV))

    def _boom(*a, **k):  # pragma: no cover - must never be called
        raise AssertionError("network must not be touched on a local miss")

    monkeypatch.setattr(erepo.httpx, "AsyncClient", _boom)
    res = asyncio.run(erepo.fetch_erepo("PAH", "c.999999999A>G"))
    assert res == {"ok": True, "found": False}


def test_live_fallback_when_dump_absent(monkeypatch, tmp_path):
    """When the dump is missing, the live API is queried and parsed."""
    _reset_state(monkeypatch, str(tmp_path / "nope.tsv"))

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json={"variantInterpretations": [{
            "@id": "/api/interpretation/CA123/MONDO:0000001/006",
            "caid": "CA123",
            "hgvs": ["NM_000257.4:c.1A>G"],
            "guidelines": [{
                "outcome": {"label": "Likely Pathogenic"},
                "agents": [{
                    "affiliation": "Cardiomyopathy VCEP",
                    "evidenceCodes": [
                        {"status": "Met", "label": "PM2"},
                        {"status": "Not Met", "label": "PP3"},
                        {"status": "Met", "label": "PS4_Moderate"},
                    ],
                }],
            }],
        }]})

    def _mock_client(*a, **k):
        k.pop("transport", None)
        return _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler),
                                  follow_redirects=True,
                                  headers=k.get("headers"))

    monkeypatch.setattr(erepo.httpx, "AsyncClient", _mock_client)
    res = asyncio.run(erepo.fetch_erepo("MYH7", "NM_000257.4:c.1A>G"))
    assert res["ok"] is True and res["found"] is True
    assert res["classification"] == "Likely Pathogenic"
    assert res["criteria"] == ["PM2", "PS4_Moderate"]
    assert res["vcep"] == "Cardiomyopathy VCEP"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))


def test_variation_gene_is_parsed_from_the_records_own_hgvs():
    assert erepo._variation_gene(
        {"#Variation": "NM_005343.3(HRAS):c.520C>T (p.Pro174Ser)"}) == "HRAS"
    assert erepo._variation_gene({"#Variation": "no gene here"}) == ""
    assert erepo._variation_gene({}) == ""


def _mislabelled_tsv(tmp_path):
    """One row exactly as eRepo ships it: labelled LRRC56, HGVS on HRAS."""
    cols = ["#Variation", "ClinVar Variation Id", "Allele Registry Id",
            "HGVS Expressions", "HGNC Gene Symbol", "Disease", "Mondo Id",
            "Mode of Inheritance", "Assertion", "Applied Evidence Codes (Met)",
            "Applied Evidence Codes (Not Met)", "Summary of interpretation",
            "PubMed Articles", "Expert Panel", "Guideline", "Approval Date",
            "Published Date", "Retracted", "Evidence Repo Link", "Uuid"]
    row = {c: "" for c in cols}
    row.update({
        "#Variation": "NM_005343.3(HRAS):c.520C>T (p.Pro174Ser)",
        "HGVS Expressions": "NM_005343.3:c.520C>T, NM_005343.2:c.520C>T",
        "HGNC Gene Symbol": "LRRC56",
        "Disease": "RASopathy",
        "Assertion": "Benign",
        "Applied Evidence Codes (Met)": "BA1",
        "Expert Panel": "RASopathy VCEP",
        "Retracted": "false",
        "Evidence Repo Link": "https://erepo.genome.network/x",
    })
    path = tmp_path / "erepo_mislabelled.tsv"
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        fh.write("\t".join(row[c] for c in cols) + "\n")
    return path


def test_a_mislabelled_record_is_found_under_the_gene_it_is_actually_in(
        monkeypatch, tmp_path):
    """The bug: curating HRAS c.520C>T could not find eRepo's own RASopathy-VCEP
    call on it, because the record was filed under a neighbouring gene."""
    _reset_state(monkeypatch, str(_mislabelled_tsv(tmp_path)))
    res = asyncio.run(erepo.fetch_erepo("HRAS", "NM_005343.3:c.520C>T"))
    assert res["found"] is True, res
    assert res["classification"] == "Benign"
    assert res["vcep"] == "RASopathy VCEP"


def test_the_label_still_resolves_too(monkeypatch, tmp_path):
    """Indexed under BOTH symbols, not replaced. Most of the other 55
    disagreements are mitochondrial or overlapping-gene rows where the label is
    plausibly the curated choice, so adding an alias must not cost a lookup that
    works today."""
    _reset_state(monkeypatch, str(_mislabelled_tsv(tmp_path)))
    res = asyncio.run(erepo.fetch_erepo("LRRC56", "NM_005343.3:c.520C>T"))
    assert res["found"] is True, res


def test_an_unrelated_gene_still_misses(monkeypatch, tmp_path):
    """Aliasing must not make the index promiscuous."""
    _reset_state(monkeypatch, str(_mislabelled_tsv(tmp_path)))
    res = asyncio.run(erepo.fetch_erepo("MYH7", "NM_005343.3:c.520C>T"))
    assert res["found"] is False, res
