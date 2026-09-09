"""Unit tests for backend.clients.gtex.

Covers both paths of the local-first GTEx client:

  (a) LOCAL: build a tiny temporary SQLite DB with the production schema,
      point the module DB_PATH at it, and assert fetch_gtex returns the
      byte-identical success shape from the local lookup.

  (b) LIVE FALLBACK: point DB_PATH at a nonexistent file and mock the two
      live GTEx Portal endpoints (/reference/gene then
      /expression/medianGeneExpression) with an httpx.MockTransport so no
      network is touched — confirming the fallback still works when the DB
      is absent, and a gene GTEx doesn't have falls through to a live miss.

Runnable with pytest (``python -m pytest backend/tests/test_gtex.py``) or
directly (``python -m backend.tests.test_gtex``).
"""
from __future__ import annotations

import asyncio
import sqlite3
import tempfile
from pathlib import Path

import httpx

import backend.clients.gtex as gtex

_REAL_ASYNC_CLIENT = httpx.AsyncClient

_DDL = """
CREATE TABLE gtex_cardiac_expression (
    gene_symbol   TEXT,
    gencode_id    TEXT,
    tissue_id     TEXT,
    median_tpm    REAL
)
"""


def _make_temp_db(rows: list[tuple]) -> Path:
    """Write a tiny gtex.db fixture with the production schema. Returns its
    path (caller is responsible for cleanup via the parent tempdir)."""
    tmpdir = tempfile.mkdtemp(prefix="gtex_test_")
    db_path = Path(tmpdir) / "gtex.db"
    conn = sqlite3.connect(db_path)
    conn.execute(_DDL)
    conn.executemany(
        "INSERT INTO gtex_cardiac_expression VALUES (?,?,?,?)", rows
    )
    conn.execute(
        "CREATE INDEX idx_gtex_gene_upper "
        "ON gtex_cardiac_expression (UPPER(gene_symbol))"
    )
    conn.commit()
    conn.close()
    return db_path


def test_local_lookup_returns_correct_shape():
    """A gene present in the local DB returns the full success dict from the
    DB without touching the network."""
    db_path = _make_temp_db([
        ("MYH7", "ENSG00000092054.13", "Heart_Left_Ventricle", 4527.28),
        ("MYH7", "ENSG00000092054.13", "Heart_Atrial_Appendage", 558.177),
        ("MYH7", "ENSG00000092054.13", "Artery_Aorta", 12.34),
        ("MYH7", "ENSG00000092054.13", "Artery_Coronary", 8.91),
    ])
    orig_db = gtex.DB_PATH
    orig_live = gtex._fetch_gtex_live

    async def _boom(gene):  # noqa: ANN001
        raise AssertionError("live fallback must NOT run on a local hit")

    gtex.DB_PATH = db_path
    gtex._fetch_gtex_live = _boom
    try:
        res = asyncio.run(gtex.fetch_gtex("MYH7"))
    finally:
        gtex.DB_PATH = orig_db
        gtex._fetch_gtex_live = orig_live

    assert res["ok"] is True, res
    assert res["found"] is True, res
    assert res["gene"] == "MYH7", res
    assert res["gencode_id"] == "ENSG00000092054.13", res
    assert res["dataset"] == "gtex_v10", res
    assert res["url"] == "https://gtexportal.org/home/gene/ENSG00000092054.13", res
    assert res["tissues"] == [
        {"tissue_id": "Heart_Left_Ventricle",
         "tissue_label": "Heart Left Ventricle",
         "median_tpm": 4527.28},
        {"tissue_id": "Heart_Atrial_Appendage",
         "tissue_label": "Heart Atrial Appendage",
         "median_tpm": 558.177},
        {"tissue_id": "Artery_Aorta",
         "tissue_label": "Artery Aorta",
         "median_tpm": 12.34},
        {"tissue_id": "Artery_Coronary",
         "tissue_label": "Artery Coronary",
         "median_tpm": 8.91},
    ], res
    assert set(res) == {
        "ok", "found", "gene", "gencode_id", "dataset", "tissues", "url"
    }, res
    print("test_local_lookup_returns_correct_shape: PASS")


def test_local_lookup_is_case_insensitive():
    """Lookup matches on UPPER(gene_symbol) — a lowercase query still hits."""
    db_path = _make_temp_db([
        ("JAG1", "ENSG00000101384.12", "Heart_Left_Ventricle", 11.3342),
        ("JAG1", "ENSG00000101384.12", "Heart_Atrial_Appendage", 18.0444),
    ])
    orig_db = gtex.DB_PATH
    gtex.DB_PATH = db_path
    try:
        res = asyncio.run(gtex.fetch_gtex("jag1"))
    finally:
        gtex.DB_PATH = orig_db
    assert res["found"] is True and res["gencode_id"] == "ENSG00000101384.12", res
    print("test_local_lookup_is_case_insensitive: PASS")


def _make_live_client(found: bool = True):
    """Return an AsyncClient wired to a MockTransport emulating the two live
    GTEx Portal endpoints fetch_gtex falls back to."""
    gencode_id = "ENSG00000175084.12"

    def handler(request: httpx.Request) -> httpx.Response:
        u = str(request.url)
        if "/reference/gene" in u:
            if not found:
                return httpx.Response(200, json={"data": []})
            return httpx.Response(200, json={"data": [
                {"geneSymbol": "DES", "gencodeId": gencode_id}]})
        if "/expression/medianGeneExpression" in u:
            return httpx.Response(200, json={"data": [
                {"tissueSiteDetailId": "Heart_Left_Ventricle", "median": 1234.5},
                {"tissueSiteDetailId": "Heart_Atrial_Appendage", "median": 678.9},
                {"tissueSiteDetailId": "Artery_Aorta", "median": 45.6},
                {"tissueSiteDetailId": "Artery_Coronary", "median": 23.4},
            ]})
        return httpx.Response(404)

    return _REAL_ASYNC_CLIENT(
        transport=httpx.MockTransport(handler), follow_redirects=True
    )


def test_live_fallback_when_db_absent():
    """With DB_PATH pointed at a nonexistent file, fetch_gtex falls back to
    the live API (here mocked) and returns the identical success shape."""
    orig_db = gtex.DB_PATH
    orig_client = gtex.httpx.AsyncClient
    gtex.DB_PATH = Path("/nonexistent/heartvar-test/gtex.db")
    gtex._GENCODE_CACHE.clear()
    gtex.httpx.AsyncClient = lambda *a, **k: _make_live_client(found=True)
    try:
        res = asyncio.run(gtex.fetch_gtex("DES"))
    finally:
        gtex.DB_PATH = orig_db
        gtex.httpx.AsyncClient = orig_client
        gtex._GENCODE_CACHE.clear()

    assert res["ok"] is True and res["found"] is True, res
    assert res["gene"] == "DES", res
    assert res["gencode_id"] == "ENSG00000175084.12", res
    assert res["dataset"] == "gtex_v10", res
    assert res["url"] == "https://gtexportal.org/home/gene/ENSG00000175084.12", res
    assert res["tissues"] == [
        {"tissue_id": "Heart_Left_Ventricle",
         "tissue_label": "Heart Left Ventricle", "median_tpm": 1234.5},
        {"tissue_id": "Heart_Atrial_Appendage",
         "tissue_label": "Heart Atrial Appendage", "median_tpm": 678.9},
        {"tissue_id": "Artery_Aorta",
         "tissue_label": "Artery Aorta", "median_tpm": 45.6},
        {"tissue_id": "Artery_Coronary",
         "tissue_label": "Artery Coronary", "median_tpm": 23.4},
    ], res
    print("test_live_fallback_when_db_absent: PASS")


def test_local_miss_falls_through_to_live():
    """A gene ABSENT from the local DB still falls through to the live API
    (so a partial cache can't permanently hide a gene GTEx actually has)."""
    db_path = _make_temp_db([
        ("MYH7", "ENSG00000092054.13", "Heart_Left_Ventricle", 4527.28),
        ("MYH7", "ENSG00000092054.13", "Heart_Atrial_Appendage", 558.177),
    ])
    orig_db = gtex.DB_PATH
    orig_client = gtex.httpx.AsyncClient
    gtex.DB_PATH = db_path
    gtex._GENCODE_CACHE.clear()
    gtex.httpx.AsyncClient = lambda *a, **k: _make_live_client(found=True)
    try:
        res = asyncio.run(gtex.fetch_gtex("DES"))
    finally:
        gtex.DB_PATH = orig_db
        gtex.httpx.AsyncClient = orig_client
        gtex._GENCODE_CACHE.clear()
    assert res["found"] is True and res["gencode_id"] == "ENSG00000175084.12", res
    print("test_local_miss_falls_through_to_live: PASS")


def test_live_fallback_gene_not_in_gtex():
    """When neither the (absent) DB nor the live API has the gene, the live
    'not found' shape is preserved (ok=True, found=False, gene, error)."""
    orig_db = gtex.DB_PATH
    orig_client = gtex.httpx.AsyncClient
    gtex.DB_PATH = Path("/nonexistent/heartvar-test/gtex.db")
    gtex._GENCODE_CACHE.clear()
    gtex.httpx.AsyncClient = lambda *a, **k: _make_live_client(found=False)
    try:
        res = asyncio.run(gtex.fetch_gtex("FAKEGENE"))
    finally:
        gtex.DB_PATH = orig_db
        gtex.httpx.AsyncClient = orig_client
        gtex._GENCODE_CACHE.clear()
    assert res["ok"] is True and res["found"] is False, res
    assert res["gene"] == "FAKEGENE", res
    assert "error" in res, res
    print("test_live_fallback_gene_not_in_gtex: PASS")


if __name__ == "__main__":
    test_local_lookup_returns_correct_shape()
    test_local_lookup_is_case_insensitive()
    test_live_fallback_when_db_absent()
    test_local_miss_falls_through_to_live()
    test_live_fallback_gene_not_in_gtex()
    print("ALL GTEX TESTS PASS")
