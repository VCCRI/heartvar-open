"""Unit tests for backend.clients.pmcoa.

Self-contained: mocks the NCBI ID-converter, PMC OA, and BioC endpoints
with an httpx.MockTransport so no network is touched. Runnable either with
pytest (``pytest backend/tests/test_pmcoa.py``) or directly
(``python -m backend.tests.test_pmcoa``) since the repo has no pytest dep.
"""
from __future__ import annotations

import asyncio

import httpx

import backend.clients.pmcoa as pmcoa
from backend.clients._cache import EXTERNAL_CACHE

_REAL_ASYNC_CLIENT = httpx.AsyncClient

_BIG_RESULTS = "The c.1504C>T variant was identified in 39 probands. " * 200
_PMCID = "PMC9999999"

_OA_FCGI_HITS: list[str] = []


def _make_client(oa_ok: bool = True, in_pmc: bool = True):
    """Return an AsyncClient wired to a MockTransport that emulates the
    three NCBI endpoints pmcoa uses."""
    def handler(request: httpx.Request) -> httpx.Response:
        u = str(request.url)
        if "idconv" in u:
            if not in_pmc:
                return httpx.Response(200, json={"records": [
                    {"pmid": 28615295, "status": "error",
                     "errmsg": "Identifier not found in PMC"}]})
            return httpx.Response(200, json={"records": [
                {"pmid": 28615295, "pmcid": _PMCID}]})
        if "oa.fcgi" in u:
            _OA_FCGI_HITS.append(u)
            return httpx.Response(404, text="Not Found")
        if "BioC_json" in u:
            if not oa_ok:
                return httpx.Response(200, text=(
                    "[Error] : No result can be found. <BR><HR><B> - "
                    "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/"
                    "</B><BR>"))
            doc = {"documents": [{
                "infons": {"year": "2017", "journal": "Circulation"},
                "passages": [
                    {"infons": {"section_type": "TITLE"}, "text": "MYBPC3 HCM study"},
                    {"infons": {"section_type": "RESULTS"}, "text": _BIG_RESULTS},
                    {"infons": {"section_type": "METHODS"}, "text": "Sanger sequencing of 100 probands."},
                ]}]}
            return httpx.Response(200, json=[doc])
        return httpx.Response(404)
    return _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler), follow_redirects=True)


def _patch(monkeypatch_target, oa_ok=True, in_pmc=True):
    """Replace pmcoa.httpx.AsyncClient with a MockTransport-backed one."""
    pmcoa.httpx.AsyncClient = lambda *a, **k: _make_client(oa_ok=oa_ok, in_pmc=in_pmc)


def _restore(orig):
    pmcoa.httpx.AsyncClient = orig


def test_pmcid_resolution_and_truncation():
    EXTERNAL_CACHE.clear()
    orig = pmcoa.httpx.AsyncClient
    _patch(orig, oa_ok=True, in_pmc=True)
    try:
        res = asyncio.run(pmcoa.fetch_pmc_fulltext(["28615295"]))
    finally:
        _restore(orig)
    rec = res["28615295"]
    assert rec["available"] is True, rec
    assert rec["pmcid"] == _PMCID, rec
    assert rec["results_text"], "results_text should be populated"
    combined = len(rec["results_text"]) + len(rec["methods_text"])
    assert combined <= pmcoa.MAX_TEXT_CHARS, combined
    assert rec["results_text"].endswith("[…truncated]"), "long results should be truncated"
    assert rec["year"] == "2017" and rec["journal"] == "Circulation", rec
    print("test_pmcid_resolution_and_truncation: PASS")


def test_not_open_access_returns_unavailable():
    """Non-OA is signalled by BioC's own "[Error]" body, not by oa.fcgi."""
    EXTERNAL_CACHE.clear()
    orig = pmcoa.httpx.AsyncClient
    _patch(orig, oa_ok=False, in_pmc=True)
    try:
        res = asyncio.run(pmcoa.fetch_pmc_fulltext(["28615295"]))
    finally:
        _restore(orig)
    rec = res["28615295"]
    assert rec["available"] is False, rec
    assert rec.get("pmcid") == _PMCID, "PMCID resolved but article not OA"
    assert "error" not in rec or rec.get("available") is False
    print("test_not_open_access_returns_unavailable: PASS")


def test_oa_fcgi_is_never_consulted():
    """The retired oa.fcgi endpoint must not be on the path at all.

    It 404s for every PMCID now, and the old `_is_open_access` pre-check read
    any non-200 as "not open access" — so the entire PMC full-text channel
    returned available:False for every article, starving PS3 and PS4 of body
    text. This test fails if a pre-flight to that host is reintroduced.
    """
    EXTERNAL_CACHE.clear()
    _OA_FCGI_HITS.clear()
    orig = pmcoa.httpx.AsyncClient
    _patch(orig, oa_ok=True, in_pmc=True)
    try:
        res = asyncio.run(pmcoa.fetch_pmc_fulltext(["28615295"]))
    finally:
        _restore(orig)
    assert _OA_FCGI_HITS == [], f"oa.fcgi was queried: {_OA_FCGI_HITS}"
    assert res["28615295"]["available"] is True, res
    print("test_oa_fcgi_is_never_consulted: PASS")


def test_not_in_pmc_returns_unavailable():
    EXTERNAL_CACHE.clear()
    orig = pmcoa.httpx.AsyncClient
    _patch(orig, in_pmc=False)
    try:
        res = asyncio.run(pmcoa.fetch_pmc_fulltext(["28615295"]))
    finally:
        _restore(orig)
    rec = res["28615295"]
    assert rec["available"] is False and rec.get("pmcid") is None, rec
    print("test_not_in_pmc_returns_unavailable: PASS")


def test_empty_input():
    assert asyncio.run(pmcoa.fetch_pmc_fulltext([])) == {}
    print("test_empty_input: PASS")


def _make_batch_client():
    """MockTransport where idconv PARSES the comma-separated ``ids`` and returns
    records in SCRAMBLED order, with integer pmids, one error record, and the
    OA/BioC endpoints keyed per-PMCID. Guards the batch map-back-by-echoed-pmid."""
    def handler(request: httpx.Request) -> httpx.Response:
        u = str(request.url)
        if "idconv" in u:
            ids = dict(request.url.params).get("ids", "").split(",")
            assert len(ids) == 3, f"expected one batched idconv call, got ids={ids!r}"
            recs = [
                {"pmid": 30000002, "pmcid": "PMC1000002"},
                {"pmid": 30000003, "status": "error", "errmsg": "not in PMC"},
                {"pmid": 30000001, "pmcid": "PMC1000001"},
            ]
            return httpx.Response(200, json={"records": recs})
        if "oa.fcgi" in u:
            pmcid = dict(request.url.params).get("id", "")
            return httpx.Response(200, text=(
                f'<OA><records returned-count="1"><record id="{pmcid}"/>'
                f'</records></OA>'))
        if "BioC_json" in u:
            pmcid = u.split("BioC_json/")[1].split("/")[0]
            doc = {"documents": [{
                "infons": {"year": "2020", "journal": "J Test"},
                "passages": [
                    {"infons": {"section_type": "TITLE"}, "text": f"paper {pmcid}"},
                    {"infons": {"section_type": "RESULTS"}, "text": f"results for {pmcid}"},
                ]}]}
            return httpx.Response(200, json=[doc])
        return httpx.Response(404)
    return _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler), follow_redirects=True)


def test_batch_idconv_maps_back_by_pmid():
    orig = pmcoa.httpx.AsyncClient
    pmcoa.httpx.AsyncClient = lambda *a, **k: _make_batch_client()
    try:
        res = asyncio.run(pmcoa.fetch_pmc_fulltext(["30000001", "30000002", "30000003"]))
    finally:
        _restore(orig)
    assert res["30000001"]["pmcid"] == "PMC1000001", res["30000001"]
    assert res["30000002"]["pmcid"] == "PMC1000002", res["30000002"]
    assert "PMC1000001" in res["30000001"]["title"]
    assert "PMC1000002" in res["30000002"]["title"]
    assert res["30000003"]["available"] is False
    assert res["30000003"].get("pmcid") is None, res["30000003"]
    print("test_batch_idconv_maps_back_by_pmid: PASS")


if __name__ == "__main__":
    test_pmcid_resolution_and_truncation()
    test_not_open_access_returns_unavailable()
    test_not_in_pmc_returns_unavailable()
    test_empty_input()
    test_batch_idconv_maps_back_by_pmid()
    print("ALL PMCOA TESTS PASS")
