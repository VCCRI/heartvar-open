"""Unit tests for the live-source hardening changes (live-source plan §4):
the PubTator3 autocomplete-call cap and the NCBI good-citizen query params.

Offline — web_get_json is monkeypatched; no network.

Runnable directly
(``python -m backend.tests.test_live_source_hardening``).
"""
from __future__ import annotations

import asyncio
import time

from backend.clients import pmcoa, pubtator3
from backend.clients.pubmed import _eutils_params
from backend.clients import _ncbi_throttle


def test_resolve_entity_caps_autocomplete_calls():
    calls: list[str] = []

    async def fake_web_get_json(client, url, params, timeout):
        calls.append(params["query"])
        return []

    orig = pubtator3.web_get_json
    pubtator3.web_get_json = fake_web_get_json
    terms = ["rs1", "GENE A1B", "GENE Ala1Bcd", "GENE c.1A>T",
             "GENE 1A>T", "GENE p.A1B"]
    try:
        result = asyncio.run(pubtator3._resolve_entity(None, "GENE", terms))
    finally:
        pubtator3.web_get_json = orig

    assert result is None
    assert len(terms) > pubtator3.MAX_RESOLVE_TERMS
    assert calls == terms[:pubtator3.MAX_RESOLVE_TERMS]


def test_eutils_params_carry_tool_and_email():
    p = _eutils_params(db="pubmed", term="MYH7")
    assert p["tool"] == _ncbi_throttle.NCBI_TOOL == "heartvar"
    assert p["email"] == _ncbi_throttle.NCBI_CONTACT_EMAIL
    assert "@" in p["email"]
    assert p["db"] == "pubmed" and p["term"] == "MYH7"


def test_user_agent_identifies_operator():
    ua = _ncbi_throttle.NCBI_USER_AGENT
    assert "heartvar" in ua.lower()
    assert _ncbi_throttle.NCBI_CONTACT_EMAIL in ua


def test_web_bucket_is_independent_of_api_key():
    """The web bucket (PubTator3/PMC/BioC) must NOT inherit the keyed eutils
    rate — those services don't honour the key, so its interval is a fixed
    conservative constant regardless of whether NCBI_API_KEY is set."""
    assert _ncbi_throttle._WEB_MIN_INTERVAL >= 0.33


def test_web_throttle_spaces_consecutive_calls():
    """Two back-to-back web_throttle acquisitions are spaced by at least the
    web interval (lower bound — robust to CI jitter)."""
    async def go():
        await _ncbi_throttle.web_throttle()
        t0 = time.monotonic()
        await _ncbi_throttle.web_throttle()
        return time.monotonic() - t0
    assert asyncio.run(go()) >= _ncbi_throttle._WEB_MIN_INTERVAL * 0.9


def test_non_eutils_clients_use_the_web_bucket():
    """Wiring guard: the non-E-utilities clients route through the conservative
    bucket, not the key-aware eutils path. A future refactor that re-couples
    them (the original block bug) trips this."""
    assert pubtator3.web_get_json is _ncbi_throttle.web_get_json
    assert not hasattr(pubtator3, "eutils_get_json")
    assert pmcoa.web_throttle is _ncbi_throttle.web_throttle


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
