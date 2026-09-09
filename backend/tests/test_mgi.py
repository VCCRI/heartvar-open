"""Unit tests for backend.clients.mgi.

Covers both paths of the local-lookup conversion:
  (a) the LOCAL path — builds a tiny temporary SQLite DB with the same
      schema as scripts/build_mgi_db.py, points MGI_DB_PATH at it (via
      rebinding mgi.DB_PATH), and asserts fetch_mgi returns the canonical
      shape from the local lookup.
  (b) the live FALLBACK — points DB_PATH at a nonexistent file and mocks
      the three live hosts (genenames.org, Alliance orthologs, Alliance
      phenotypes) with an httpx.MockTransport, asserting the live chain
      still produces the identical dict shape.

Self-contained: no network is touched. Runnable with
``python -m pytest backend/tests/test_mgi.py`` or directly
(``python -m backend.tests.test_mgi``).
"""
from __future__ import annotations

import asyncio
import sqlite3
import tempfile
from pathlib import Path

import httpx

import backend.clients.mgi as mgi

_REAL_ASYNC_CLIENT = httpx.AsyncClient

_FOUND_KEYS = {
    "ok", "found", "gene", "hgnc_id", "mouse_symbol", "mgi_id",
    "ortholog_confidence", "is_best_score", "phenotype_count",
    "cardiac_phenotypes", "phenotype_sample", "pubmed_count", "url",
    "alliance_url",
}


def _build_fixture_db(path: Path) -> None:
    """Build a tiny mgi.db with the same schema as scripts/build_mgi_db.py."""
    conn = sqlite3.connect(path)
    conn.execute(
        "CREATE TABLE gene_ortholog ("
        "human_symbol TEXT PRIMARY KEY, hgnc_id TEXT, mgi_id TEXT, "
        "mouse_symbol TEXT)"
    )
    conn.execute(
        "CREATE TABLE gene_phenotype ("
        "mgi_id TEXT NOT NULL, mp_id TEXT, phenotype_statement TEXT, "
        "pubmed_id TEXT)"
    )
    conn.execute(
        "INSERT INTO gene_ortholog VALUES (?,?,?,?)",
        ("MYH7", "HGNC:7577", "MGI:2155600", "Myh7"),
    )
    conn.executemany(
        "INSERT INTO gene_phenotype VALUES (?,?,?,?)",
        [
            ("MGI:2155600", "MP:0000266", "abnormal heart morphology", "11111111"),
            ("MGI:2155600", "MP:0000266", "abnormal heart morphology", "22222222"),
            ("MGI:2155600", "MP:0002834", "decreased cardiac muscle contractility", "33333333"),
            ("MGI:2155600", "MP:0001262", "decreased body weight", "44444444"),
        ],
    )
    conn.commit()
    conn.close()


def test_local_path_returns_canonical_shape():
    orig_db = mgi.DB_PATH
    with tempfile.TemporaryDirectory() as tmp:
        db = Path(tmp) / "mgi.db"
        _build_fixture_db(db)
        mgi.DB_PATH = db
        try:
            res = asyncio.run(mgi.fetch_mgi("MYH7"))
        finally:
            mgi.DB_PATH = orig_db

    assert set(res.keys()) == _FOUND_KEYS, res
    assert res["ok"] is True and res["found"] is True, res
    assert res["gene"] == "MYH7"
    assert res["hgnc_id"] == "HGNC:7577"
    assert res["mgi_id"] == "MGI:2155600"
    assert res["mouse_symbol"] == "Myh7"
    assert res["ortholog_confidence"] is None
    assert res["is_best_score"] is None
    assert res["phenotype_count"] == 3, res
    assert res["cardiac_phenotypes"] == [
        "abnormal heart morphology",
        "decreased cardiac muscle contractility",
    ], res
    assert res["pubmed_count"] == 4, res
    assert res["url"] == "https://www.informatics.jax.org/marker/MGI:2155600"
    assert res["alliance_url"] == "https://www.alliancegenome.org/gene/MGI:2155600"
    print("test_local_path_returns_canonical_shape: PASS")


def test_local_miss_falls_back_to_live():
    """Gene absent from gene_ortholog → local miss → live fallback fires.
    The fixture DB exists but has no row for the queried gene, so the code
    must reach the (mocked) live chain."""
    orig_db = mgi.DB_PATH
    orig_client = mgi.httpx.AsyncClient
    with tempfile.TemporaryDirectory() as tmp:
        db = Path(tmp) / "mgi.db"
        _build_fixture_db(db)
        mgi.DB_PATH = db
        mgi.httpx.AsyncClient = lambda *a, **k: _make_live_client()
        try:
            res = asyncio.run(mgi.fetch_mgi("TTN"))
        finally:
            mgi.DB_PATH = orig_db
            mgi.httpx.AsyncClient = orig_client
    assert res["found"] is True, res
    assert res["mgi_id"] == "MGI:98864", res
    assert res["ortholog_confidence"] == "good", res
    print("test_local_miss_falls_back_to_live: PASS")


def _make_live_client():
    """AsyncClient wired to a MockTransport emulating the three live MGI
    hosts: genenames.org HGNC fetch, Alliance orthologs, Alliance
    phenotypes."""
    def handler(request: httpx.Request) -> httpx.Response:
        u = str(request.url)
        if request.url.host == "rest.genenames.org":
            return httpx.Response(200, json={
                "response": {"docs": [{"hgnc_id": "HGNC:12403"}]}
            })
        if u.endswith("/orthologs"):
            return httpx.Response(200, json={"results": [
                {"geneToGeneOrthologyGenerated": {
                    "objectGene": {
                        "primaryExternalId": "MGI:98864",
                        "geneSymbol": {"displayText": "Ttn"},
                        "taxon": {"curie": "NCBITaxon:10090"},
                    },
                    "confidence": {"name": "good"},
                    "isBestScore": {"name": "Yes"},
                }}
            ]})
        if "/phenotypes" in u:
            return httpx.Response(200, json={"results": [
                {"phenotypeStatement": "abnormal cardiac muscle morphology",
                 "pubmedPubModIDs": ["PMID:12345678"]},
                {"phenotypeStatement": "decreased body length",
                 "pubmedPubModIDs": ["87654321"]},
            ]})
        return httpx.Response(404)
    return _REAL_ASYNC_CLIENT(
        transport=httpx.MockTransport(handler), follow_redirects=True
    )


def test_live_fallback_when_db_absent():
    """DB_PATH points at a nonexistent file → local path short-circuits to
    a miss → the (mocked) live chain produces the canonical shape."""
    orig_db = mgi.DB_PATH
    orig_client = mgi.httpx.AsyncClient
    mgi._HGNC_CACHE.clear()
    mgi._ORTHO_CACHE.clear()
    mgi.DB_PATH = Path("/nonexistent/path/to/mgi.db")
    mgi.httpx.AsyncClient = lambda *a, **k: _make_live_client()
    try:
        res = asyncio.run(mgi.fetch_mgi("TTN"))
    finally:
        mgi.DB_PATH = orig_db
        mgi.httpx.AsyncClient = orig_client
        mgi._HGNC_CACHE.clear()
        mgi._ORTHO_CACHE.clear()

    assert set(res.keys()) == _FOUND_KEYS, res
    assert res["ok"] is True and res["found"] is True, res
    assert res["gene"] == "TTN"
    assert res["hgnc_id"] == "HGNC:12403"
    assert res["mgi_id"] == "MGI:98864"
    assert res["mouse_symbol"] == "Ttn"
    assert res["ortholog_confidence"] == "good", res
    assert res["is_best_score"] is True, res
    assert res["phenotype_count"] == 2, res
    assert res["cardiac_phenotypes"] == ["abnormal cardiac muscle morphology"], res
    assert res["pubmed_count"] == 2, res
    assert res["url"] == "https://www.informatics.jax.org/marker/MGI:98864"
    print("test_live_fallback_when_db_absent: PASS")


if __name__ == "__main__":
    test_local_path_returns_canonical_shape()
    test_local_miss_falls_back_to_live()
    test_live_fallback_when_db_absent()
    print("ALL MGI TESTS PASS")
