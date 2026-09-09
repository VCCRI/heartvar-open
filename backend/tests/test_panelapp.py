"""Unit tests for backend/clients/panelapp.py public helpers + fetch_panelapp.

Complements (does NOT duplicate) the two existing PanelApp test files:
  - test_panelapp_snapshot.py    — the snapshot-hit / gene-absent / live-fallback
                                    acquisition paths of fetch_panelapp.
  - test_panelapp_hpo_descendants.py — the offline cardiac-HPO descendant closure
                                    and the live-JAX descendant fallback.

This file covers the OTHER public surface:
  - _label              confidence-level → green/amber/red/unknown mapping
  - _categorize_panel   name → cardiovascular category bucketing (specific
                        sub-buckets win over broad fallbacks; non-cardiac → None)
  - _normalise_phrase / _substring_phrase_match  phrase-key normalisation + the
                        ≥6-char longest-substring fallback
  - _parse_phenotype_input / _parse_hpo_terms    mixed HPO-ID + free-text parsing,
                        parent bubbling, and unrecognised-token surfacing
  - check_hpo_relevance the descendant / seed-map / parent-term / phrase match
                        ladder (incl. the seed-map fallback when the cache is empty)
  - _panelapp_results   the live error + transient-None (request_with_retry None)
                        acquisition paths
  - fetch_panelapp      amber/red confidence rollup, panel-id dedup, sorting,
                        per-category rollup, and the propagated error path.

All HTTP is mocked (httpx.MockTransport) or monkeypatched out — NO real network.

Runnable with pytest or directly (python -m backend.tests.test_panelapp).
"""
from __future__ import annotations

import asyncio

import httpx
import pytest

import backend.clients.panelapp as panelapp

_REAL_ASYNC_CLIENT = httpx.AsyncClient


def _reset_snapshot(monkeypatch, path: str | None) -> None:
    """Force the snapshot loader to re-resolve, point it at `path` (or unset),
    and clear the descendant cache so relevance checks start from a known state.
    Mirrors the helper in test_panelapp_snapshot.py."""
    panelapp._panelapp_snapshot = None
    panelapp._panelapp_snapshot_loaded = False
    panelapp.CARDIAC_HPO_DESCENDANTS.clear()
    if path is None:
        monkeypatch.delenv("PANELAPP_SNAPSHOT_PATH", raising=False)
    else:
        monkeypatch.setenv("PANELAPP_SNAPSHOT_PATH", path)


def _no_descendant_cache(monkeypatch) -> None:
    """Pin _ensure_descendant_cache to a no-op so fetch_panelapp never touches
    the (possibly absent) closure file or live JAX. Leaves CARDIAC_HPO_DESCENDANTS
    in whatever state the test set."""
    async def _noop() -> None:
        return None
    monkeypatch.setattr(panelapp, "_ensure_descendant_cache", _noop)


def test_label_maps_numeric_confidence_levels():
    assert panelapp._label("3") == "green"
    assert panelapp._label("2") == "amber"
    assert panelapp._label("1") == "red"
    assert panelapp._label(3) == "green"
    assert panelapp._label(2) == "amber"
    assert panelapp._label(1) == "red"


def test_label_none_is_unknown_and_unmapped_passes_through():
    assert panelapp._label(None) == "unknown"
    assert panelapp._label("0") == "0"
    assert panelapp._label("99") == "99"


def test_categorize_specific_subtypes_win_over_broad_buckets():
    assert panelapp._categorize_panel("Hypertrophic cardiomyopathy") == "hcm"
    assert panelapp._categorize_panel("Dilated Cardiomyopathy") == "dcm"
    assert panelapp._categorize_panel(
        "Arrhythmogenic right ventricular cardiomyopathy") == "cardiomyopathy_other"
    assert panelapp._categorize_panel("Long QT syndrome") == "channelopathy"
    assert panelapp._categorize_panel("Brugada syndrome") == "channelopathy"


def test_categorize_broad_fallback_buckets():
    assert panelapp._categorize_panel(
        "Paediatric or syndromic cardiomyopathy") == "Cardiomyopathy"
    assert panelapp._categorize_panel("Congenital heart disease") == "Congenital heart disease"
    assert panelapp._categorize_panel("Thoracic aortic aneurysm") == "Aortic / vascular"
    assert panelapp._categorize_panel("Marfan syndrome") == "Aortic / vascular"


def test_categorize_cardiovascular_but_no_subbucket_is_other_cardiac():
    assert panelapp._categorize_panel("Cardiac something unusual") == "Other cardiac"


def test_categorize_non_cardiac_returns_none():
    assert panelapp._categorize_panel("Intellectual disability") is None
    assert panelapp._categorize_panel("") is None
    assert panelapp._categorize_panel(None) is None


def test_normalise_phrase_lowercases_hyphens_whitespace_punct():
    assert panelapp._normalise_phrase(
        "  Non-Compaction  Cardiomyopathy. ") == "non compaction cardiomyopathy"
    assert panelapp._normalise_phrase("non-compaction") == "non compaction"
    assert panelapp._normalise_phrase("non compaction") == "non compaction"


def test_substring_phrase_match_returns_longest_key():
    assert panelapp._substring_phrase_match("hypertrophic cardiomyopathy 1") == ("hcm",)


def test_substring_phrase_match_short_keys_never_match():
    assert panelapp._SUBSTRING_MIN_KEY_LEN == 6
    assert panelapp._substring_phrase_match("the hcm gene") is None
    assert panelapp._substring_phrase_match("cadherin variant") is None
    assert panelapp._substring_phrase_match("") is None


def test_parse_phenotype_input_mixes_hpo_phrase_and_unrecognised():
    hpo_ids, phrase_cats, recognised, unrecognised = panelapp._parse_phenotype_input(
        "HP:0001639, dilated cardiomyopathy, gobbledygook"
    )
    assert hpo_ids == ["HP:0001639"]
    assert "dcm" in phrase_cats
    assert "Cardiomyopathy" in phrase_cats
    assert "HP:0001639" in recognised
    assert "dilated cardiomyopathy" in recognised
    assert unrecognised == ["gobbledygook"]


def test_parse_phenotype_input_accepts_list_and_semicolons_uppercases_hpo():
    hpo_ids, _, _, _ = panelapp._parse_phenotype_input(["hp:0001644", "DCM"])
    assert hpo_ids == ["HP:0001644"]
    hpo_ids2, _, _, _ = panelapp._parse_phenotype_input("HP:0001639; HP:0001644")
    assert hpo_ids2 == ["HP:0001639", "HP:0001644"]


def test_parse_phenotype_input_empty_inputs():
    assert panelapp._parse_phenotype_input(None) == ([], set(), [], [])
    assert panelapp._parse_phenotype_input("") == ([], set(), [], [])
    assert panelapp._parse_phenotype_input("   ,  ; ") == ([], set(), [], [])


def test_parse_hpo_terms_legacy_alias_returns_only_ids():
    assert panelapp._parse_hpo_terms("HP:0001639; DCM") == ["HP:0001639"]
    assert panelapp._parse_hpo_terms("dilated cardiomyopathy") == []


def test_relevance_parent_term_matches_every_category():
    panelapp.CARDIAC_HPO_DESCENDANTS.clear()
    assert panelapp.check_hpo_relevance("HP:0001627", "Aortic / vascular") is True
    assert panelapp.check_hpo_relevance("HP:0011675", "Congenital heart disease") is True


def test_relevance_free_text_phrase_match():
    panelapp.CARDIAC_HPO_DESCENDANTS.clear()
    assert panelapp.check_hpo_relevance("dilated cardiomyopathy", "dcm") is True
    assert panelapp.check_hpo_relevance("dilated cardiomyopathy", "channelopathy") is False


def test_relevance_no_terms_returns_false():
    panelapp.CARDIAC_HPO_DESCENDANTS.clear()
    assert panelapp.check_hpo_relevance("", "hcm") is False
    assert panelapp.check_hpo_relevance("totally unknown phenotype text", "hcm") is False


def test_relevance_descendant_cache_path():
    panelapp.CARDIAC_HPO_DESCENDANTS.clear()
    panelapp.CARDIAC_HPO_DESCENDANTS["HP:9999999"] = {"hcm"}
    try:
        assert panelapp.check_hpo_relevance("HP:9999999", "hcm") is True
        assert panelapp.check_hpo_relevance("HP:9999999", "dcm") is False
    finally:
        panelapp.CARDIAC_HPO_DESCENDANTS.clear()


def test_relevance_seed_map_fallback_when_cache_empty():
    panelapp.CARDIAC_HPO_DESCENDANTS.clear()
    assert "HP:0005157" in (panelapp._load_hpo_map().get("hcm") or [])
    assert panelapp.check_hpo_relevance("HP:0005157", "hcm") is True
    assert panelapp.check_hpo_relevance("HP:0005157", "Aortic / vascular") is False


def test_panelapp_results_live_non_200_returns_error(monkeypatch):
    _reset_snapshot(monkeypatch, "/nonexistent/snapshot.json")

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(404, text="not found")

    monkeypatch.setattr(
        panelapp.httpx, "AsyncClient",
        lambda *a, **k: _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler)),
    )
    results, error = asyncio.run(panelapp._panelapp_results("MYH7"))
    assert results is None
    assert error is not None and error.startswith("404")


def test_panelapp_results_transient_none_returns_error(monkeypatch):
    _reset_snapshot(monkeypatch, "/nonexistent/snapshot.json")

    async def _none(*a, **k):
        return None

    monkeypatch.setattr(panelapp, "request_with_retry", _none)
    monkeypatch.setattr(
        panelapp.httpx, "AsyncClient",
        lambda *a, **k: _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(
            lambda r: httpx.Response(200, json={"results": []}))),
    )
    results, error = asyncio.run(panelapp._panelapp_results("MYH7"))
    assert results is None
    assert error == "PanelApp transient failure after retries"


def test_panelapp_results_snapshot_never_errors_for_absent_gene(monkeypatch, tmp_path):
    snap = tmp_path / "snap.json"
    snap.write_text('{"genes": {"MYH7": [{"panel": {"id": 1, "name": "X"}}]}}')
    _reset_snapshot(monkeypatch, str(snap))
    monkeypatch.setattr(panelapp.httpx, "AsyncClient",
                        lambda *a, **k: (_ for _ in ()).throw(
                            AssertionError("snapshot path must not touch network")))
    results, error = asyncio.run(panelapp._panelapp_results("ABSENT"))
    assert results == []
    assert error is None


def test_fetch_panelapp_amber_red_rollup_dedup_and_sort(monkeypatch, tmp_path):
    """Snapshot with a green + amber + red panel plus a duplicate panel id.
    Verifies confidence counts, panel-id dedup, green-first sort ordering, and
    the per-category rollup flags."""
    snap = tmp_path / "snap.json"
    snap.write_text("""{"genes": {"MYH7": [
        {"panel": {"id": 10, "name": "Dilated Cardiomyopathy", "version": "1"},
         "confidence_level": "2", "mode_of_inheritance": "MONOALLELIC"},
        {"panel": {"id": 11, "name": "Long QT syndrome", "version": "1"},
         "confidence_level": "1", "mode_of_inheritance": "MONOALLELIC"},
        {"panel": {"id": 12, "name": "Hypertrophic cardiomyopathy", "version": "2"},
         "confidence_level": "3", "mode_of_inheritance": "MONOALLELIC"},
        {"panel": {"id": 12, "name": "Hypertrophic cardiomyopathy", "version": "2"},
         "confidence_level": "3", "mode_of_inheritance": "MONOALLELIC"},
        {"panel": {"id": 99, "name": "Intellectual disability", "version": "1"},
         "confidence_level": "3", "mode_of_inheritance": "MONOALLELIC"}
    ]}}""")
    _reset_snapshot(monkeypatch, str(snap))
    _no_descendant_cache(monkeypatch)
    monkeypatch.setattr(panelapp.httpx, "AsyncClient",
                        lambda *a, **k: (_ for _ in ()).throw(
                            AssertionError("snapshot path must not touch network")))

    res = asyncio.run(panelapp.fetch_panelapp("MYH7", hpo=None))
    assert res["ok"] is True
    assert res["total_panels"] == 3
    assert res["green_panels"] == 1
    assert res["amber_panels"] == 1
    assert res["red_panels"] == 1
    assert res["total_panels_overall"] == 5
    assert res["panels_found"][0]["confidence"] == "green"
    assert res["best_confidence_label"] == "green"
    confidences = [p["confidence"] for p in res["panels_found"]]
    assert confidences == ["green", "amber", "red"]
    assert res["any_green_with_hpo_match"] is False
    assert res["any_green_no_hpo_match"] is True
    assert all(p["contributes_to_pp4"] is False for p in res["panels_found"])
    assert res["categories"]["hcm"]["any_green"] is True
    assert res["categories"]["hcm"]["panel_count"] == 1
    assert res["categories"]["dcm"]["any_green"] is False


def test_fetch_panelapp_green_with_matching_phenotype_contributes_pp4(monkeypatch, tmp_path):
    snap = tmp_path / "snap.json"
    snap.write_text("""{"genes": {"MYH7": [
        {"panel": {"id": 20, "name": "Hypertrophic cardiomyopathy", "version": "1"},
         "confidence_level": "3", "mode_of_inheritance": "MONOALLELIC"}
    ]}}""")
    _reset_snapshot(monkeypatch, str(snap))
    _no_descendant_cache(monkeypatch)
    monkeypatch.setattr(panelapp.httpx, "AsyncClient",
                        lambda *a, **k: (_ for _ in ()).throw(
                            AssertionError("snapshot path must not touch network")))

    res = asyncio.run(panelapp.fetch_panelapp("MYH7", hpo="hypertrophic cardiomyopathy"))
    assert res["panels_found"][0]["hpo_match"] is True
    assert res["panels_found"][0]["contributes_to_pp4"] is True
    assert res["any_green_with_hpo_match"] is True
    assert res["any_green_no_hpo_match"] is False
    assert "hcm" in res["submitted_phrase_categories"]


def test_fetch_panelapp_propagates_error_from_failed_live_call(monkeypatch):
    _reset_snapshot(monkeypatch, "/nonexistent/snapshot.json")
    _no_descendant_cache(monkeypatch)

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(503, text="service unavailable")

    async def _err(*a, **k):
        return httpx.Response(503, text="service unavailable")

    monkeypatch.setattr(panelapp, "request_with_retry", _err)
    monkeypatch.setattr(
        panelapp.httpx, "AsyncClient",
        lambda *a, **k: _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler)),
    )
    res = asyncio.run(panelapp.fetch_panelapp("MYH7", hpo=None))
    assert res["ok"] is False
    assert res["error"].startswith("503")


def test_fetch_panelapp_live_success_and_unrecognised_token(monkeypatch):
    """Live fallback (no snapshot): a CHD-bucket green panel, plus an
    unrecognised phenotype token surfaced back to the caller."""
    _reset_snapshot(monkeypatch, "/nonexistent/snapshot.json")
    _no_descendant_cache(monkeypatch)

    def handler(request: httpx.Request) -> httpx.Response:
        assert "/genes/" in str(request.url)
        return httpx.Response(200, json={"results": [{
            "panel": {"id": 30, "name": "Congenital heart disease", "version": "3"},
            "confidence_level": "3",
            "mode_of_inheritance": "MONOALLELIC",
            "phenotypes": ["CHD"],
            "entity_status": "green",
        }]})

    monkeypatch.setattr(
        panelapp.httpx, "AsyncClient",
        lambda *a, **k: _REAL_ASYNC_CLIENT(transport=httpx.MockTransport(handler),
                                           follow_redirects=True),
    )
    res = asyncio.run(panelapp.fetch_panelapp("GATA4", hpo="tetralogy of fallot, mystery word"))
    assert res["ok"] is True
    assert res["on_chd_panel"] is True
    assert res["on_cardiovascular_panel"] is True
    assert res["green_panels"] == 1
    assert res["unrecognised_hpo_tokens"] == ["mystery word"]
    assert res["panels_found"][0]["contributes_to_pp4"] is True


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
