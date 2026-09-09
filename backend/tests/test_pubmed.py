"""Unit tests for backend.clients.pubmed.

Exercises the two public coroutines with the network fully mocked (no real
NCBI E-utilities traffic):

  * ``fetch_pubmed``            — variant-specific esearch + efetch
  * ``get_gene_disease_literature`` — gene-disease esearch + efetch

Mocking style follows test_medgen.py: the client builds its own
``httpx.AsyncClient`` inside the function body, so we patch
``pubmed.httpx.AsyncClient`` to ride an ``httpx.MockTransport`` that emulates
the esearch.fcgi (JSON) and efetch.fcgi (XML) endpoints. The error path
monkeypatches the throttled-GET helpers to raise ``NCBIError`` /
``httpx.TimeoutException`` directly, so no retry-backoff ``asyncio.sleep``
fires and the test stays fast and deterministic.

Both public functions are cached wrappers; conftest's autouse fixture clears
EXTERNAL_CACHE before each test, and the ``__main__`` runner clears it inline.

Run:
    .venv/bin/python -m pytest backend/tests/test_pubmed.py -q
"""
from __future__ import annotations

import asyncio

import httpx

import backend.clients.pubmed as pubmed

_REAL_ASYNC_CLIENT = httpx.AsyncClient


# structured AbstractText (Label attribute), an ArticleDate year, an
_EFETCH_XML = """<?xml version="1.0"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>11111111</PMID>
      <Article>
        <Journal>
          <ISOAbbreviation>Circulation</ISOAbbreviation>
          <JournalIssue><PubDate><Year>2019</Year></PubDate></JournalIssue>
        </Journal>
        <ArticleTitle>MYH7 c.1208G&gt;A in hypertrophic cardiomyopathy</ArticleTitle>
        <Abstract>
          <AbstractText Label="BACKGROUND">A founder variant was studied.</AbstractText>
          <AbstractText Label="RESULTS">It segregated with disease.</AbstractText>
        </Abstract>
        <AuthorList>
          <Author><LastName>Smith</LastName><Initials>AB</Initials></Author>
          <Author><LastName>Jones</LastName><Initials>C</Initials></Author>
        </AuthorList>
        <ArticleDate><Year>2019</Year><Month>03</Month><Day>01</Day></ArticleDate>
      </Article>
    </MedlineCitation>
    <PubmedData>
      <ArticleIdList>
        <ArticleId IdType="pubmed">11111111</ArticleId>
        <ArticleId IdType="doi">10.1161/CIRC.0001</ArticleId>
      </ArticleIdList>
    </PubmedData>
  </PubmedArticle>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>22222222</PMID>
      <Article>
        <Journal><Title>Journal of Cardiology</Title>
          <JournalIssue><PubDate><Year>2021</Year></PubDate></JournalIssue>
        </Journal>
        <ArticleTitle>Functional analysis of a cardiomyopathy allele</ArticleTitle>
        <Abstract>
          <AbstractText>An unstructured abstract paragraph.</AbstractText>
        </Abstract>
        <AuthorList>
          <Author><CollectiveName>The HCM Study Group</CollectiveName></Author>
        </AuthorList>
      </Article>
    </MedlineCitation>
    <PubmedData>
      <ArticleIdList>
        <ArticleId IdType="pubmed">22222222</ArticleId>
      </ArticleIdList>
    </PubmedData>
  </PubmedArticle>
</PubmedArticleSet>
"""


def _make_client(*, idlist, count, xml=_EFETCH_XML):
    """Build a MockTransport-backed AsyncClient that emulates esearch.fcgi
    (JSON, returning ``idlist``/``count``) and efetch.fcgi (XML, ``xml``)."""

    def handler(request: httpx.Request) -> httpx.Response:
        u = str(request.url)
        if "esearch.fcgi" in u:
            return httpx.Response(200, json={
                "esearchresult": {"idlist": list(idlist), "count": str(count)}
            })
        if "efetch.fcgi" in u:
            return httpx.Response(
                200, text=xml, headers={"Content-Type": "text/xml"}
            )
        return httpx.Response(404)  # pragma: no cover

    return _REAL_ASYNC_CLIENT(
        transport=httpx.MockTransport(handler), follow_redirects=True
    )


def _patch_client(monkeypatch, client):
    """Make every ``httpx.AsyncClient(...)`` built inside pubmed return
    ``client``. The client ignores the headers kwarg the production code
    passes — MockTransport short-circuits the transport layer."""
    monkeypatch.setattr(
        pubmed.httpx, "AsyncClient", lambda *a, **k: client
    )


def test_fetch_pubmed_successful_parse(monkeypatch):
    """A representative esearch+efetch round trip parses both articles and
    carries the true (uncapped) hit count."""
    _patch_client(
        monkeypatch,
        _make_client(idlist=["11111111", "22222222"], count=14),
    )
    res = asyncio.run(
        pubmed.fetch_pubmed("MYH7", "NM_000257.4:c.1208G>A", "Arg403Gln")
    )

    assert res["ok"] is True, res
    assert res["gene"] == "MYH7"
    assert res["hgvs_c"] == "NM_000257.4:c.1208G>A"
    assert res["amino_acid"] == "Arg403Gln"
    assert '"c.1208G>A"[tiab]' in res["variant_query"]
    assert '"Arg403Gln"[tiab]' in res["variant_query"]
    assert "MYH7[tiab]" in res["variant_query"]
    assert res["variant_total"] == 14
    assert res["variant_search_url"].startswith(
        "https://pubmed.ncbi.nlm.nih.gov/?term="
    )

    papers = res["variant_papers"]
    assert len(papers) == 2, papers
    p0, p1 = papers
    assert p0["pmid"] == "11111111"
    assert p0["title"] == "MYH7 c.1208G>A in hypertrophic cardiomyopathy"
    assert p0["first_author"] == "Smith AB"
    assert p0["year"] == "2019"
    assert p0["journal"] == "Circulation"
    assert p0["doi"] == "10.1161/CIRC.0001"
    assert p0["url"] == "https://pubmed.ncbi.nlm.nih.gov/11111111/"
    assert "BACKGROUND: A founder variant" in p0["abstract"]
    assert "RESULTS: It segregated" in p0["abstract"]
    assert p1["first_author"] == "The HCM Study Group"
    assert p1["journal"] == "Journal of Cardiology"
    assert p1["doi"] == ""
    assert p1["year"] == "2021"


def test_fetch_pubmed_empty_results(monkeypatch):
    """esearch returns zero hits → no efetch, empty paper list, ok stays True."""
    _patch_client(monkeypatch, _make_client(idlist=[], count=0))
    res = asyncio.run(pubmed.fetch_pubmed("MYH7", "c.9999G>A", "Arg9999Gln"))

    assert res["ok"] is True, res
    assert res["variant_papers"] == []
    assert res["variant_total"] == 0
    assert res["variant_query"]
    assert res["variant_search_url"]


def test_fetch_pubmed_no_variant_query_short_circuits(monkeypatch):
    """No HGVS-c and no amino acid → no variant-specific query to run; returns
    the empty result without touching the network (transport asserts loudly)."""

    def _boom(request):  # pragma: no cover - must never be hit
        raise AssertionError(f"unexpected HTTP call: {request.url}")

    _patch_client(
        monkeypatch,
        _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(_boom)),
    )
    res = asyncio.run(pubmed.fetch_pubmed("MYH7", "", None))

    assert res["ok"] is True, res
    assert res["variant_query"] == ""
    assert res["variant_papers"] == []
    assert res["variant_total"] == 0


def test_fetch_pubmed_timeout_path(monkeypatch):
    """A timeout raised inside the esearch surfaces as ok=False + error, never
    propagating. Patch _esearch directly so no retry-backoff sleeps fire."""

    async def _raise_timeout(client, query):
        raise httpx.TimeoutException("read timed out")

    monkeypatch.setattr(pubmed, "_esearch", _raise_timeout)
    _patch_client(monkeypatch, _make_client(idlist=["1"], count=1))

    res = asyncio.run(pubmed.fetch_pubmed("MYH7", "c.1208G>A", "Arg403Gln"))
    assert res["ok"] is False, res
    assert "error" in res
    assert "TimeoutException" in res["error"]
    assert res["variant_papers"] == []


def test_gene_disease_literature_successful_parse(monkeypatch):
    """Gene-disease search parses papers and adds the 300-char abstract
    snippet the gene-tab prompt quotes."""
    _patch_client(
        monkeypatch,
        _make_client(idlist=["11111111", "22222222"], count=2),
    )
    res = asyncio.run(
        pubmed.get_gene_disease_literature(
            "MYH7", ["hypertrophic cardiomyopathy"]
        )
    )

    assert res["ok"] is True, res
    assert res["gene"] == "MYH7"
    assert res["phenotype_keywords"] == ["hypertrophic cardiomyopathy"]
    assert "MYH7[tiab]" in res["query"]
    assert '"hypertrophic cardiomyopathy"[tiab]' in res["query"]
    assert res["search_url"].startswith("https://pubmed.ncbi.nlm.nih.gov/?term=")

    papers = res["papers"]
    assert len(papers) == 2, papers
    assert papers[0]["pmid"] == "11111111"
    assert "abstract_snippet" in papers[0]
    assert papers[0]["abstract_snippet"]
    assert len(papers[0]["abstract_snippet"]) <= pubmed.GENE_LIT_SNIPPET_CHARS + 1


def test_gene_disease_literature_empty_phenotype_uses_fallback(monkeypatch):
    """No phenotype keywords → falls back to the broad cardiac vocabulary, and
    an empty esearch yields an empty (but ok) result."""
    _patch_client(monkeypatch, _make_client(idlist=[], count=0))
    res = asyncio.run(pubmed.get_gene_disease_literature("MYH7", []))

    assert res["ok"] is True, res
    assert res["papers"] == []
    assert res["phenotype_keywords"] == list(pubmed._DEFAULT_CARDIAC_KEYWORDS)
    assert '"congenital heart disease"[tiab]' in res["query"]


def test_gene_disease_literature_missing_gene_error():
    """Empty gene symbol short-circuits with the documented error shape, no
    network."""
    res = asyncio.run(pubmed.get_gene_disease_literature("", ["x"]))
    assert res["ok"] is False, res
    assert res["error"] == "gene_symbol is required"
    assert res["papers"] == []
    assert res["query"] == ""


def test_gene_disease_literature_esearch_error(monkeypatch):
    """An NCBIError from the gene-literature esearch surfaces as ok=False with
    an ``esearch ...`` error, without raising. Patch the helper so no retry
    sleeps run."""
    from backend.clients._ncbi_throttle import NCBIError

    async def _raise(client, url, params, timeout):
        raise NCBIError("HTTP 429")

    monkeypatch.setattr(pubmed, "eutils_get_json", _raise)
    _patch_client(monkeypatch, _make_client(idlist=["1"], count=1))

    res = asyncio.run(
        pubmed.get_gene_disease_literature("MYH7", ["cardiomyopathy"])
    )
    assert res["ok"] is False, res
    assert res["error"].startswith("esearch "), res
    assert res["papers"] == []


if __name__ == "__main__":
    from backend.clients._cache import EXTERNAL_CACHE

    class _MP:
        """Minimal monkeypatch shim for the direct-run path."""

        def __init__(self):
            self._attr = []

        def setattr(self, obj, name, val):
            self._attr.append((obj, name, getattr(obj, name)))
            setattr(obj, name, val)

        def undo(self):
            for obj, name, old in reversed(self._attr):
                setattr(obj, name, old)

    tests = [
        test_fetch_pubmed_successful_parse,
        test_fetch_pubmed_empty_results,
        test_fetch_pubmed_no_variant_query_short_circuits,
        test_fetch_pubmed_timeout_path,
        test_gene_disease_literature_successful_parse,
        test_gene_disease_literature_empty_phenotype_uses_fallback,
        test_gene_disease_literature_esearch_error,
    ]
    for fn in tests:
        EXTERNAL_CACHE.clear()
        mp = _MP()
        try:
            fn(mp)
        finally:
            mp.undo()
        print(f"  ok  {fn.__name__}")
    EXTERNAL_CACHE.clear()
    test_gene_disease_literature_missing_gene_error()
    print("  ok  test_gene_disease_literature_missing_gene_error")
    print("\nALL PUBMED TESTS PASS")
