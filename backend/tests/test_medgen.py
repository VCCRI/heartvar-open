"""Unit tests for backend.clients.medgen.

Covers BOTH lookup paths of the local-first MedGen client:

  (a) Local DB hit — build a tiny temporary SQLite DB with the same
      schema as scripts/build_medgen_db.py, point MEDGEN_DB_PATH at it,
      and assert fetch_medgen returns the canonical shape (cardiac-first
      ordering, MAX_CONDITIONS cap, total_found pre-cap, URL builders).

  (b) Live fallback — with the local DB absent (MEDGEN_DB_PATH pointed at
      a nonexistent file AND medgen.DB_PATH monkeypatched), the three
      NCBI E-utilities calls (esearch + elink + esummary) are mocked with
      an httpx.MockTransport and the same canonical shape is produced.

Self-contained: no network is touched. Runnable with pytest
(``python -m pytest backend/tests/test_medgen.py``) or directly as a
module from the repo root (``python -m backend.tests.test_medgen``).
"""
from __future__ import annotations

import asyncio
import os
import sqlite3
from pathlib import Path

import httpx

import backend.clients.medgen as medgen

_REAL_ASYNC_CLIENT = httpx.AsyncClient


def _build_fixture_db(path: Path) -> None:
    """Create a tiny medgen.db with the same schema the build script writes.

    Two cardiac and two non-cardiac conditions for gene FOO, inserted in
    an order that does NOT already have the cardiac entries first, so the
    cardiac-first stable sort is genuinely exercised. A fifth+sixth
    condition push the count past MAX_CONDITIONS (5) to exercise the cap.
    """
    conn = sqlite3.connect(path)
    conn.execute(
        """
        CREATE TABLE conditions (
            gene_symbol   TEXT,
            cui           TEXT,
            mim           TEXT,
            name          TEXT,
            source        TEXT,
            last_updated  TEXT
        )
        """
    )
    rows = [
        ("foo", "C0000001", "100001", "Distal skeletal myopathy", "MONDO", "Jan 1 2020"),
        ("foo", "C0000002", None, "Hypertrophic cardiomyopathy 1", "MONDO", "Jan 1 2020"),
        ("foo", "C0000003", "100003", "Sensorineural hearing loss", "NCBI curation", "Jan 1 2020"),
        ("foo", "C0000004", "100004", "Dilated cardiomyopathy 1S", "MONDO", "Jan 1 2020"),
        ("foo", "C0000005", "100005", "Retinitis pigmentosa", "MONDO", "Jan 1 2020"),
        ("foo", "C0000006", "100006", "Polydactyly syndrome", "MONDO", "Jan 1 2020"),
    ]
    conn.executemany("INSERT INTO conditions VALUES (?,?,?,?,?,?)", rows)
    conn.execute(
        "CREATE INDEX idx_conditions_gene ON conditions (UPPER(gene_symbol))"
    )
    conn.commit()
    conn.close()


def _make_live_client():
    """AsyncClient wired to a MockTransport emulating the three E-utilities
    endpoints the live MedGen fallback hits."""
    def handler(request: httpx.Request) -> httpx.Response:
        u = str(request.url)
        if "esearch.fcgi" in u:
            return httpx.Response(200, json={
                "esearchresult": {"idlist": ["4625"]}
            })
        if "elink.fcgi" in u:
            return httpx.Response(200, json={
                "linksets": [{
                    "linksetdbs": [{
                        "dbto": "medgen",
                        "links": ["111", "222", "333"],
                    }]
                }]
            })
        if "esummary.fcgi" in u:
            return httpx.Response(200, json={
                "result": {
                    "uids": ["111", "222", "333"],
                    "111": {
                        "conceptid": "C9000001",
                        "title": "Skeletal myopathy NOS",
                        "semantictype": {"value": "Disease or Syndrome"},
                        "conceptmeta": "",
                    },
                    "222": {
                        "conceptid": "C9000002",
                        "title": "Hypertrophic cardiomyopathy 1",
                        "semantictype": {"value": "Disease or Syndrome"},
                        "conceptmeta": (
                            '<ConceptMeta>'
                            '<Name SAB="OMIM" SDUI="192600" TTY="PT">'
                            'Hypertrophic cardiomyopathy 1</Name>'
                            '</ConceptMeta>'
                        ),
                    },
                    "333": {
                        "conceptid": "C9000003",
                        "title": "Some gene placeholder",
                        "semantictype": {"value": "Gene or Genome"},
                        "conceptmeta": "",
                    },
                }
            })
        return httpx.Response(404)

    return _REAL_ASYNC_CLIENT(
        transport=httpx.MockTransport(handler), follow_redirects=True
    )


def test_local_db_hit_shape_sort_and_cap(tmp_path):
    db = tmp_path / "medgen.db"
    _build_fixture_db(db)

    prev = os.environ.get("MEDGEN_DB_PATH")
    os.environ["MEDGEN_DB_PATH"] = str(db)
    try:
        res = asyncio.run(medgen.fetch_medgen("FOO"))
    finally:
        if prev is None:
            os.environ.pop("MEDGEN_DB_PATH", None)
        else:
            os.environ["MEDGEN_DB_PATH"] = prev

    assert res["ok"] is True, res
    assert res["gene"] == "FOO", res
    assert set(res.keys()) == {"ok", "gene", "conditions", "total_found", "url"}, res
    assert res["total_found"] == 6, res
    assert len(res["conditions"]) == medgen.MAX_CONDITIONS == 5, res
    assert res["url"] == "https://www.ncbi.nlm.nih.gov/medgen/?term=FOO[gene]", res

    cond_keys = {"name", "cui", "mim", "moi", "heart_related", "omim_url", "medgen_url"}
    for c in res["conditions"]:
        assert set(c.keys()) == cond_keys, c
        assert c["moi"] is None

    names = [c["name"] for c in res["conditions"]]
    assert names[0] == "Hypertrophic cardiomyopathy 1", names
    assert names[1] == "Dilated cardiomyopathy 1S", names
    assert res["conditions"][0]["heart_related"] is True
    assert res["conditions"][1]["heart_related"] is True
    assert res["conditions"][2]["heart_related"] is False

    hcm = res["conditions"][0]
    assert hcm["cui"] == "C0000002"
    assert hcm["mim"] is None
    assert hcm["omim_url"] is None
    assert hcm["medgen_url"] == "https://www.ncbi.nlm.nih.gov/medgen/C0000002"
    dcm = res["conditions"][1]
    assert dcm["mim"] == "100004"
    assert dcm["omim_url"] == "https://www.omim.org/entry/100004"
    print("test_local_db_hit_shape_sort_and_cap: PASS")


def test_local_db_case_insensitive_lookup(tmp_path):
    """Lowercased / mixed-case gene symbols resolve against the same rows."""
    db = tmp_path / "medgen.db"
    _build_fixture_db(db)
    prev = os.environ.get("MEDGEN_DB_PATH")
    os.environ["MEDGEN_DB_PATH"] = str(db)
    try:
        res = asyncio.run(medgen.fetch_medgen("foo"))
    finally:
        if prev is None:
            os.environ.pop("MEDGEN_DB_PATH", None)
        else:
            os.environ["MEDGEN_DB_PATH"] = prev
    assert res["ok"] is True and res["total_found"] == 6, res
    assert res["gene"] == "foo", res
    print("test_local_db_case_insensitive_lookup: PASS")


def test_live_fallback_when_db_absent(tmp_path, monkeypatch):
    missing = tmp_path / "does_not_exist.db"
    monkeypatch.setenv("MEDGEN_DB_PATH", str(missing))
    monkeypatch.setattr(medgen, "DB_PATH", missing)
    assert not missing.exists()

    orig = medgen.httpx.AsyncClient
    medgen.httpx.AsyncClient = lambda *a, **k: _make_live_client()
    try:
        res = asyncio.run(medgen.fetch_medgen("MYH7"))
    finally:
        medgen.httpx.AsyncClient = orig

    assert res["ok"] is True, res
    assert res["gene"] == "MYH7", res
    assert set(res.keys()) == {"ok", "gene", "conditions", "total_found", "url"}, res
    assert res["total_found"] == 2, res
    names = [c["name"] for c in res["conditions"]]
    assert "Some gene placeholder" not in names, "SemanticType filter failed"
    assert names[0] == "Hypertrophic cardiomyopathy 1", names
    hcm = res["conditions"][0]
    assert hcm["heart_related"] is True
    assert hcm["cui"] == "C9000002"
    assert hcm["mim"] == "192600", hcm
    assert hcm["omim_url"] == "https://www.omim.org/entry/192600", hcm
    assert hcm["medgen_url"] == "https://www.ncbi.nlm.nih.gov/medgen/C9000002", hcm
    cond_keys = {"name", "cui", "mim", "moi", "heart_related", "omim_url", "medgen_url"}
    assert set(hcm.keys()) == cond_keys, hcm
    print("test_live_fallback_when_db_absent: PASS")


def test_missing_gene_returns_error_shape():
    """Empty gene short-circuits with the legacy error shape on both paths."""
    res = asyncio.run(medgen.fetch_medgen(""))
    assert res["ok"] is False, res
    assert res["error"] == "missing gene symbol", res
    assert res["url"] == "https://www.ncbi.nlm.nih.gov/medgen/?term=[gene]", res
    print("test_missing_gene_returns_error_shape: PASS")


if __name__ == "__main__":
    import tempfile

    class _MP:
        """Minimal monkeypatch shim for the direct-run path."""
        def __init__(self):
            self._env = []
            self._attr = []
        def setenv(self, k, v):
            self._env.append((k, os.environ.get(k)))
            os.environ[k] = v
        def setattr(self, obj, name, val):
            self._attr.append((obj, name, getattr(obj, name)))
            setattr(obj, name, val)
        def undo(self):
            for k, old in reversed(self._env):
                if old is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = old
            for obj, name, old in reversed(self._attr):
                setattr(obj, name, old)

    with tempfile.TemporaryDirectory() as d:
        tp = Path(d)
        test_local_db_hit_shape_sort_and_cap(tp)
    with tempfile.TemporaryDirectory() as d:
        tp = Path(d)
        test_local_db_case_insensitive_lookup(tp)
    with tempfile.TemporaryDirectory() as d:
        tp = Path(d)
        mp = _MP()
        try:
            test_live_fallback_when_db_absent(tp, mp)
        finally:
            mp.undo()
    test_missing_gene_returns_error_shape()
    print("ALL MEDGEN TESTS PASS")
