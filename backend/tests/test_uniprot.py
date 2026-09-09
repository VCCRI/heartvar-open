"""Unit tests for backend.clients.uniprot's local-DB + live-fallback paths.

Self-contained:
  - The local-path tests build a tiny temporary SQLite DB with the same
    schema ``scripts/build_uniprot_db.py`` produces, point UNIPROT_DB_PATH
    at it, and assert ``fetch_uniprot`` returns the correct shape from the
    DB — including the byte-for-byte payload that ``parse_entry`` would
    have produced from a streamed entry.
  - The live-fallback test points UNIPROT_DB_PATH at a nonexistent file and
    mocks the ``/uniprotkb/search`` HTTP call with an httpx.MockTransport,
    asserting the DB-vs-live results are identical for the same entry.

Runnable with pytest (``python -m pytest backend/tests/test_uniprot.py``)
or directly (``python -m backend.tests.test_uniprot``) since the repo's
clients have no hard pytest dependency.
"""
from __future__ import annotations

import asyncio
import json
import os
import sqlite3
import tempfile
from pathlib import Path

import httpx

import backend.clients.uniprot as uniprot

_REAL_ASYNC_CLIENT = httpx.AsyncClient


_ENTRY = {
    "primaryAccession": "P12883",
    "uniProtkbId": "MYH7_HUMAN",
    "proteinDescription": {
        "recommendedName": {"fullName": {"value": "Myosin-7"}}
    },
    "sequence": {"length": 1935},
    "genes": [
        {
            "geneName": {"value": "MYH7"},
            "synonyms": [{"value": "MYHCB"}],
        }
    ],
    "features": [
        {
            "type": "Domain",
            "location": {"start": {"value": 32}, "end": {"value": 81}},
            "description": "Myosin N-terminal SH3-like",
        },
        {
            "type": "Binding site",
            "location": {"start": {"value": 179}, "end": {"value": 186}},
            "description": "ATP",
        },
        {
            "type": "Natural variant",
            "location": {"start": {"value": 403}, "end": {"value": 403}},
            "description": "in CMH1; pathogenic; dbSNP:rs121913624",
            "alternativeSequence": {
                "originalSequence": "R",
                "alternativeSequences": ["Q"],
            },
            "featureId": "VAR_007552",
            "featureCrossReferences": [
                {"database": "dbSNP", "id": "rs121913624"},
            ],
        },
        {
            "type": "Natural variant",
            "location": {"start": {"value": 3}, "end": {"value": 3}},
            "description": "in dbSNP:rs3729993",
            "alternativeSequence": {
                "originalSequence": "D",
                "alternativeSequences": ["A"],
            },
            "featureId": "VAR_029430",
            "featureCrossReferences": [
                {"database": "dbSNP", "id": "rs3729993"},
            ],
        },
    ],
}


def _build_tiny_db(path: Path) -> None:
    """Create a one-entry uniprot.db with the same schema + index the build
    script produces, indexing the entry under primary + synonym (uppercased)
    and storing the parse_entry payload (without the 'gene' key)."""
    conn = sqlite3.connect(path)
    conn.execute(
        "CREATE TABLE uniprot_entry (gene_symbol TEXT, accession TEXT, payload TEXT)"
    )
    payload = uniprot.parse_entry(_ENTRY)
    payload.pop("gene", None)
    payload_json = json.dumps(payload, separators=(",", ":"))
    primary, synonyms = uniprot.gene_symbols(_ENTRY)
    seen: set[str] = set()
    for sym in ([primary] if primary else []) + synonyms:
        up = sym.upper()
        if up in seen:
            continue
        seen.add(up)
        conn.execute(
            "INSERT INTO uniprot_entry (gene_symbol, accession, payload) VALUES (?,?,?)",
            (up, payload.get("accession"), payload_json),
        )
    conn.execute("CREATE INDEX idx_uniprot_gene ON uniprot_entry (gene_symbol)")
    conn.commit()
    conn.close()


def test_local_db_returns_correct_shape():
    """fetch_uniprot reads the local DB and returns the full structured
    shape (ok/found/gene/accession/uniprot_id/.../natural_variants)."""
    with tempfile.TemporaryDirectory() as td:
        db = Path(td) / "uniprot.db"
        _build_tiny_db(db)
        prev = os.environ.get("UNIPROT_DB_PATH")
        os.environ["UNIPROT_DB_PATH"] = str(db)
        try:
            res = asyncio.run(uniprot.fetch_uniprot("MYH7"))
        finally:
            if prev is None:
                os.environ.pop("UNIPROT_DB_PATH", None)
            else:
                os.environ["UNIPROT_DB_PATH"] = prev

    assert res["ok"] is True and res["found"] is True, res
    assert res["gene"] == "MYH7", res
    assert res["accession"] == "P12883", res
    assert res["uniprot_id"] == "MYH7_HUMAN", res
    assert res["protein_name"] == "Myosin-7", res
    assert res["length"] == 1935, res
    assert res["url"] == "https://www.uniprot.org/uniprotkb/P12883", res
    assert res["domains"] == [
        {"name": "Myosin N-terminal SH3-like", "start": 32, "end": 81}
    ], res
    assert res["active_sites"] == [], res
    assert res["binding_sites"] == [
        {"start": 179, "end": 186, "description": "ATP"}
    ], res
    assert {f["type"] for f in res["features"]} == {"Domain", "Binding site"}, res
    assert res["natural_variant_count"] == 2, res
    by_pos = {v["position"]: v for v in res["natural_variants"]}
    assert by_pos[403]["clinical_significance"] == "pathogenic", by_pos[403]
    assert by_pos[403]["original_aa"] == "R" and by_pos[403]["variant_aa"] == "Q"
    assert by_pos[403]["dbsnp"] == "rs121913624", by_pos[403]
    assert by_pos[3]["clinical_significance"] == "tolerated", by_pos[3]
    assert set(res.keys()) == {
        "ok", "found", "gene", "accession", "uniprot_id", "protein_name",
        "length", "sequence", "url", "domains", "active_sites",
        "binding_sites", "features", "natural_variant_count",
        "natural_variants",
    }, sorted(res.keys())
    print("test_local_db_returns_correct_shape: PASS")


def test_synonym_lookup_hits_local_db_and_echoes_symbol():
    """A synonym (MYHCB) resolves to the same entry via the local DB, and
    the returned 'gene' echoes the REQUESTED symbol (matching the live
    gene_exact behaviour) while accession stays canonical."""
    with tempfile.TemporaryDirectory() as td:
        db = Path(td) / "uniprot.db"
        _build_tiny_db(db)
        prev = os.environ.get("UNIPROT_DB_PATH")
        os.environ["UNIPROT_DB_PATH"] = str(db)
        try:
            res = asyncio.run(uniprot.fetch_uniprot("myhcb"))
        finally:
            if prev is None:
                os.environ.pop("UNIPROT_DB_PATH", None)
            else:
                os.environ["UNIPROT_DB_PATH"] = prev
    assert res["found"] is True, res
    assert res["gene"] == "myhcb", res
    assert res["accession"] == "P12883", res
    print("test_synonym_lookup_hits_local_db_and_echoes_symbol: PASS")


def _mock_search_client() -> httpx.AsyncClient:
    """An AsyncClient wired to a MockTransport that emulates the
    /uniprotkb/search endpoint, returning the synthetic entry."""
    def handler(request: httpx.Request) -> httpx.Response:
        if "uniprotkb/search" in str(request.url):
            return httpx.Response(200, json={"results": [_ENTRY]})
        return httpx.Response(404)
    return _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler))


def test_live_fallback_when_db_absent_matches_local():
    """When the DB file is absent, fetch_uniprot falls back to the live
    /uniprotkb/search path. The fallback result is byte-identical to the
    local-DB result for the same entry (DB-vs-live parity)."""
    with tempfile.TemporaryDirectory() as td:
        db = Path(td) / "uniprot.db"
        _build_tiny_db(db)
        prev = os.environ.get("UNIPROT_DB_PATH")
        os.environ["UNIPROT_DB_PATH"] = str(db)
        try:
            local_res = asyncio.run(uniprot.fetch_uniprot("MYH7"))
        finally:
            if prev is None:
                os.environ.pop("UNIPROT_DB_PATH", None)
            else:
                os.environ["UNIPROT_DB_PATH"] = prev

    prev = os.environ.get("UNIPROT_DB_PATH")
    os.environ["UNIPROT_DB_PATH"] = str(Path(tempfile.gettempdir()) / "does_not_exist_uniprot.db")
    orig_client = uniprot.httpx.AsyncClient
    uniprot.httpx.AsyncClient = lambda *a, **k: _mock_search_client()
    try:
        live_res = asyncio.run(uniprot.fetch_uniprot("MYH7"))
    finally:
        uniprot.httpx.AsyncClient = orig_client
        if prev is None:
            os.environ.pop("UNIPROT_DB_PATH", None)
        else:
            os.environ["UNIPROT_DB_PATH"] = prev

    assert live_res["found"] is True, live_res
    assert live_res == local_res, (
        "live-fallback result must be byte-identical to the local-DB result\n"
        f"local={local_res}\nlive={live_res}"
    )
    print("test_live_fallback_when_db_absent_matches_local: PASS")


def test_live_fallback_not_found():
    """When the DB is absent and the live search returns no results,
    fetch_uniprot returns the {ok, found:False, gene} shape unchanged."""
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json={"results": []})

    prev = os.environ.get("UNIPROT_DB_PATH")
    os.environ["UNIPROT_DB_PATH"] = str(Path(tempfile.gettempdir()) / "does_not_exist_uniprot.db")
    orig_client = uniprot.httpx.AsyncClient
    uniprot.httpx.AsyncClient = lambda *a, **k: _REAL_ASYNC_CLIENT(
        transport=httpx.MockTransport(handler)
    )
    try:
        res = asyncio.run(uniprot.fetch_uniprot("ZZZ_NOT_A_GENE"))
    finally:
        uniprot.httpx.AsyncClient = orig_client
        if prev is None:
            os.environ.pop("UNIPROT_DB_PATH", None)
        else:
            os.environ["UNIPROT_DB_PATH"] = prev
    assert res == {"ok": True, "found": False, "gene": "ZZZ_NOT_A_GENE"}, res
    print("test_live_fallback_not_found: PASS")


if __name__ == "__main__":
    test_local_db_returns_correct_shape()
    test_synonym_lookup_hits_local_db_and_echoes_symbol()
    test_live_fallback_when_db_absent_matches_local()
    test_live_fallback_not_found()
    print("ALL UNIPROT TESTS PASS")
