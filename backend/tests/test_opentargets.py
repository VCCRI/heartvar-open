"""Unit tests for backend.clients.opentargets_client.

Drives the public ``fetch_opentargets`` coroutine with the network fully
mocked (no real Open Targets GraphQL traffic). The client builds its async
client via ``_http_retry.make_async_client`` (imported into the module
namespace), so we patch ``opentargets_client.make_async_client`` to return an
``httpx.AsyncClient`` riding an ``httpx.MockTransport`` that answers the single
``POST /graphql`` the client issues.

Covered:
  * Successful HPO-matched parse (HCM → matched disease + datatype scores).
  * Fallback parse when no HPO maps (top-5 diseases, matched_disease null).
  * Empty-result path (gene has zero associated diseases).
  * Error path — a non-200 GraphQL response is terminal in
    ``request_with_retry``, so ``_post`` returns None and ``fetch_opentargets``
    yields the ``ok: False`` shape WITHOUT any retry-backoff sleeps.

The error path uses a 404 (rather than a transport raise or a 5xx) precisely so
``request_with_retry`` returns immediately — keeping the test fast/offline.

Run:
    .venv/bin/python -m pytest backend/tests/test_opentargets.py -q
"""
from __future__ import annotations

import asyncio
import json

import httpx

import backend.clients.opentargets_client as ot
from backend.clients._gene_ids import resolve_ensembl_gene_id

_REAL_ASYNC_CLIENT = httpx.AsyncClient

_ENSG_MYH7 = "ENSG00000092054"
_ENSG_ELN = "ENSG00000049540"


def _graphql_response(rows, approved_symbol="MYH7"):
    """Wrap ``rows`` in the GraphQL envelope the client parses:
    data.target.associatedDiseases.rows[]."""
    return {
        "data": {
            "target": {
                "approvedSymbol": approved_symbol,
                "associatedDiseases": {"rows": rows},
            }
        }
    }


def _patch_post(monkeypatch, *, status=200, json_body=None, text=None):
    """Patch ot.make_async_client so the single POST the client issues is
    served by a MockTransport returning the given status/body."""

    def handler(request: httpx.Request) -> httpx.Response:
        assert request.method == "POST", request.method
        if text is not None:
            return httpx.Response(status, text=text)
        return httpx.Response(status, json=json_body if json_body is not None else {})

    def factory(*a, **k):
        return _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler))

    monkeypatch.setattr(ot, "make_async_client", factory)


def test_fetch_opentargets_hpo_matched(monkeypatch):
    """HP:0001639 maps to HCM (EFO_0000408); the matching row drives the
    matched-disease branch with overall score + datatype breakdown."""
    rows = [
        {
            "disease": {"id": "EFO_0000408", "name": "hypertrophic cardiomyopathy"},
            "score": 0.876543,
            "datatypeScores": [
                {"id": "genetic_association", "score": 0.9},
                {"id": "literature", "score": 0.4},
            ],
        },
        {
            "disease": {"id": "EFO_9999999", "name": "some other disease"},
            "score": 0.5,
            "datatypeScores": [{"id": "animal_model", "score": 0.2}],
        },
    ]
    _patch_post(monkeypatch, json_body=_graphql_response(rows))

    res = asyncio.run(ot.fetch_opentargets(_ENSG_MYH7, ["HP:0001639"]))

    assert res["ok"] is True, res
    assert res["gene_symbol"] == "MYH7"
    assert res["ensembl_id"] == _ENSG_MYH7
    assert res["matched_hpo"] == "HP:0001639"
    assert res["matched_disease"] == {
        "id": "EFO_0000408", "name": "hypertrophic cardiomyopathy",
    }
    assert res["overall_association_score"] == 0.877
    assert res["top_diseases"] == []
    dts = res["datatype_scores"]
    assert dts["genetic_association"] == 0.9
    assert dts["literature"] == 0.4
    for dt in ot._DATATYPES:
        assert dt in dts
    assert dts["somatic_mutation"] == 0.0
    assert res["evidence_summary"].startswith("Strong")
    assert "hypertrophic cardiomyopathy" in res["evidence_summary"]
    assert res["url"] == f"https://platform.opentargets.org/target/{_ENSG_MYH7}"


def test_fetch_opentargets_fallback_top_diseases(monkeypatch):
    """No matching HPO (none supplied) → fallback to the top-5 diseases for the
    gene, matched_disease null, top_diseases populated."""
    rows = [
        {
            "disease": {"id": "EFO_0000408", "name": "hypertrophic cardiomyopathy"},
            "score": 0.6,
            "datatypeScores": [{"id": "genetic_association", "score": 0.55}],
        },
        {
            "disease": {"id": "EFO_0000400", "name": "dilated cardiomyopathy"},
            "score": 0.42,
            "datatypeScores": [{"id": "literature", "score": 0.3}],
        },
    ]
    _patch_post(monkeypatch, json_body=_graphql_response(rows))

    res = asyncio.run(ot.fetch_opentargets(_ENSG_MYH7, None))

    assert res["ok"] is True, res
    assert res["matched_hpo"] is None
    assert res["matched_disease"] is None
    assert res["overall_association_score"] == 0.6
    assert len(res["top_diseases"]) == 2
    assert res["top_diseases"][0] == {
        "id": "EFO_0000408", "name": "hypertrophic cardiomyopathy", "score": 0.6,
    }
    assert res["evidence_summary"].startswith("Moderate")
    assert "hypertrophic cardiomyopathy" in res["evidence_summary"]


def test_fetch_opentargets_empty_rows(monkeypatch):
    """Target exists but no associated diseases → ok=True with the
    'no associations indexed' summary and null scores."""
    _patch_post(monkeypatch, json_body=_graphql_response([]))

    res = asyncio.run(ot.fetch_opentargets(_ENSG_MYH7, ["HP:0001639"]))

    assert res["ok"] is True, res
    assert res["matched_disease"] is None
    assert res["overall_association_score"] is None
    assert res["top_diseases"] == []
    assert "No Open Targets associations indexed" in res["evidence_summary"]
    for dt in ot._DATATYPES:
        assert res["datatype_scores"][dt] == 0.0


def test_fetch_opentargets_no_target(monkeypatch):
    """A GraphQL data envelope with target=null → no data → ok=False shape."""
    _patch_post(monkeypatch, json_body={"data": {"target": None}})

    res = asyncio.run(ot.fetch_opentargets(_ENSG_MYH7, None))
    assert res["ok"] is False, res
    assert res["error"] == f"Open Targets has no target record for {_ENSG_MYH7}"


def test_fetch_opentargets_http_error(monkeypatch):
    """A 404 is terminal in request_with_retry → _post returns None → the
    module helper returns the 'unavailable' error shape, never raising."""
    _patch_post(monkeypatch, status=404, text="not found")

    res = asyncio.run(ot.fetch_opentargets(_ENSG_MYH7, ["HP:0001639"]))
    assert res["ok"] is False, res
    assert res["error"] == "Open Targets returned HTTP 404"


def test_fetch_opentargets_graphql_errors_field(monkeypatch):
    """A 200 response carrying a populated GraphQL ``errors`` array is treated
    as no data → ok=False (no exception)."""
    body = {"errors": [{"message": "boom"}], "data": None}
    _patch_post(monkeypatch, json_body=body)

    res = asyncio.run(ot.fetch_opentargets(_ENSG_MYH7, None))
    assert res["ok"] is False, res
    assert res["error"].startswith("Open Targets GraphQL error:")


def test_fetch_opentargets_bad_ensembl_id_no_network(monkeypatch):
    """An unresolvable gene id with no symbol to fall back on short-circuits to
    None → ok=False, without issuing any HTTP (transport asserts loudly)."""

    def _boom(request):  # pragma: no cover - must never be hit
        raise AssertionError(f"unexpected HTTP call: {request.url}")

    monkeypatch.setattr(
        ot, "make_async_client",
        lambda *a, **k: _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(_boom)),
    )
    res = asyncio.run(ot.fetch_opentargets("NOT_AN_ENSG", None))
    assert res["ok"] is False, res
    assert "No Ensembl gene ID resolved" in res["error"], res
    assert "NOT_AN_ENSG" in res["error"], res


def test_entrez_gene_id_from_refseq_vep_resolves_and_queries(monkeypatch):
    """gene_id='2006' (Entrez, ELN) + symbol resolves to ENSG00000049540 and
    the GraphQL query is actually issued against the resolved Ensembl id."""
    seen: dict = {}

    def handler(request: httpx.Request) -> httpx.Response:
        seen["ensemblId"] = json.loads(request.content)["variables"]["ensemblId"]
        return httpx.Response(200, json=_graphql_response(
            [{
                "disease": {"id": "MONDO_0008504",
                            "name": "supravalvular aortic stenosis"},
                "score": 0.783,
                "datatypeScores": [{"id": "genetic_association", "score": 0.878}],
            }],
            approved_symbol="ELN",
        ))

    monkeypatch.setattr(
        ot, "make_async_client",
        lambda *a, **k: _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler)),
    )

    res = asyncio.run(ot.fetch_opentargets("2006", None, gene_symbol="ELN"))
    assert res["ok"] is True, res
    assert seen["ensemblId"] == _ENSG_ELN, seen
    assert res["ensembl_id"] == _ENSG_ELN
    assert res["url"].endswith(_ENSG_ELN)
    assert res["top_diseases"][0]["name"] == "supravalvular aortic stenosis"


def test_missing_gene_id_falls_back_to_symbol(monkeypatch):
    """VEP supplying no gene_id at all must still resolve via the symbol."""
    seen: dict = {}

    def handler(request: httpx.Request) -> httpx.Response:
        seen["ensemblId"] = json.loads(request.content)["variables"]["ensemblId"]
        return httpx.Response(200, json=_graphql_response([], approved_symbol="ELN"))

    monkeypatch.setattr(
        ot, "make_async_client",
        lambda *a, **k: _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler)),
    )

    res = asyncio.run(ot.fetch_opentargets(None, None, gene_symbol="ELN"))
    assert res["ok"] is True, res
    assert seen["ensemblId"] == _ENSG_ELN, seen


def test_resolver_maps_entrez_and_symbol_offline():
    """The resolver itself is offline-only (HGNC complete set), so it works
    under HEARTVAR_OFFLINE_STRICT. Direct unit check of the three input forms."""
    assert resolve_ensembl_gene_id("2006", "ELN") == _ENSG_ELN
    assert resolve_ensembl_gene_id("2006", None) == _ENSG_ELN
    assert resolve_ensembl_gene_id(None, "ELN") == _ENSG_ELN
    assert resolve_ensembl_gene_id("4625", "MYH7") == _ENSG_MYH7
    assert resolve_ensembl_gene_id(_ENSG_ELN + ".16", None) == _ENSG_ELN
    assert resolve_ensembl_gene_id("NOT_AN_ENSG", None) is None
    assert resolve_ensembl_gene_id(None, "NOT_A_REAL_GENE_XYZ") is None


if __name__ == "__main__":

    class _MP:
        def __init__(self):
            self._attr = []

        def setattr(self, obj, name, val):
            self._attr.append((obj, name, getattr(obj, name)))
            setattr(obj, name, val)

        def undo(self):
            for obj, name, old in reversed(self._attr):
                setattr(obj, name, old)

    tests = [
        test_fetch_opentargets_hpo_matched,
        test_fetch_opentargets_fallback_top_diseases,
        test_fetch_opentargets_empty_rows,
        test_fetch_opentargets_no_target,
        test_fetch_opentargets_http_error,
        test_fetch_opentargets_graphql_errors_field,
        test_fetch_opentargets_bad_ensembl_id_no_network,
        test_entrez_gene_id_from_refseq_vep_resolves_and_queries,
        test_missing_gene_id_falls_back_to_symbol,
        test_resolver_maps_entrez_and_symbol_offline,
    ]
    for fn in tests:
        mp = _MP()
        try:
            fn(mp) if fn.__code__.co_argcount else fn()
        finally:
            mp.undo()
        print(f"  ok  {fn.__name__}")
    print("\nALL OPENTARGETS TESTS PASS")
