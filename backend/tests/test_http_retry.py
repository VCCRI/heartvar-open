"""Unit tests for the shared HTTP retry helper's Retry-After honoring
(live-source hardening plan §4.1 part A).

request_with_retry now waits at least as long as a server-supplied Retry-After
on a 429/5xx (capped), and falls through to the static backoff when the header is
absent — so the 8 hosts that don't send it are unaffected.

No network — httpx.MockTransport + a captured asyncio.sleep.
Runnable directly (``python -m backend.tests.test_http_retry``).
"""
from __future__ import annotations

import asyncio

import httpx

from backend.clients import _http_retry
from backend.clients._http_retry import _retry_after_seconds, request_with_retry, _RETRY_AFTER_CAP


def test_retry_after_delta_seconds():
    r = httpx.Response(429, headers={"Retry-After": "2"})
    assert _retry_after_seconds(r) == 2.0


def test_retry_after_capped():
    r = httpx.Response(429, headers={"Retry-After": "999"})
    assert _retry_after_seconds(r) == _RETRY_AFTER_CAP


def test_retry_after_absent_or_bad():
    assert _retry_after_seconds(httpx.Response(429)) is None
    assert _retry_after_seconds(httpx.Response(429, headers={"Retry-After": "soon"})) is None
    assert _retry_after_seconds(httpx.Response(429, headers={"Retry-After": "0"})) is None


def _run_with_captured_sleep(handler):
    """Drive request_with_retry against a MockTransport handler, capturing the
    backoff sleep durations. throttle_host is stubbed to a no-op so its own
    host-spacing sleeps don't land in the captured list (patching asyncio.sleep
    hits the shared module, so the throttle would otherwise leak into it)."""
    slept: list[float] = []

    async def fake_sleep(d):
        slept.append(d)

    async def no_throttle(_url):
        return None

    real_sleep = _http_retry.asyncio.sleep
    real_throttle = _http_retry.throttle_host
    _http_retry.asyncio.sleep = fake_sleep
    _http_retry.throttle_host = no_throttle
    client = httpx.AsyncClient(transport=httpx.MockTransport(handler))
    try:
        resp = asyncio.run(
            request_with_retry(client, "GET", "https://rest.ensembl.org/x",
                               timeout=5.0, name="test")
        )
    finally:
        _http_retry.asyncio.sleep = real_sleep
        _http_retry.throttle_host = real_throttle
        asyncio.run(client.aclose())
    return resp, slept


def test_backoff_waits_for_retry_after():
    state = {"n": 0}

    def handler(_req):
        state["n"] += 1
        if state["n"] == 1:
            return httpx.Response(429, headers={"Retry-After": "5"})
        return httpx.Response(200, text="ok")

    resp, slept = _run_with_captured_sleep(handler)
    assert resp is not None and resp.status_code == 200
    assert slept and slept[0] == 5.0, slept


def test_backoff_falls_back_without_header():
    state = {"n": 0}

    def handler(_req):
        state["n"] += 1
        if state["n"] == 1:
            return httpx.Response(503)
        return httpx.Response(200, text="ok")

    resp, slept = _run_with_captured_sleep(handler)
    assert resp is not None and resp.status_code == 200
    assert slept and 0.5 <= slept[0] < 0.8, slept


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


def _run_counting(handler, *, deadline=None, per_attempt=0.0):
    """Drive request_with_retry, counting attempts and measuring real elapsed.

    Backoff sleeps are stubbed out (they would dominate wall-clock) but the
    HANDLER consumes real time, so the deadline sees a clock that actually
    advances — which is the thing under test.
    """
    import time as _time
    calls = {"n": 0}

    async def fake_sleep(d):
        return None

    async def no_throttle(_url):
        return None

    def counting_handler(request):
        calls["n"] += 1
        if per_attempt:
            _time.sleep(per_attempt)
        return handler(request)

    real_sleep = _http_retry.asyncio.sleep
    real_throttle = _http_retry.throttle_host
    _http_retry.asyncio.sleep = fake_sleep
    _http_retry.throttle_host = no_throttle
    client = httpx.AsyncClient(transport=httpx.MockTransport(counting_handler))
    kw = {} if deadline is None else {"deadline": deadline}
    t0 = _time.perf_counter()
    try:
        resp = asyncio.run(request_with_retry(
            client, "GET", "https://rest.ensembl.org/x", timeout=5.0,
            name="test", **kw))
    finally:
        _http_retry.asyncio.sleep = real_sleep
        _http_retry.throttle_host = real_throttle
        asyncio.run(client.aclose())
    return resp, calls["n"], _time.perf_counter() - t0


def _boom(request):
    raise httpx.ConnectError("name resolution timed out", request=request)


def test_deadline_stops_retrying_before_max_attempts():
    """The production shape: every attempt fails slowly on a cold container."""
    resp, attempts, _ = _run_counting(_boom, deadline=0.25, per_attempt=0.12)
    assert resp is None, "an exhausted retry still returns None"
    assert attempts < 4, (
        f"made all {attempts} attempts despite the deadline — retries are "
        "still unbounded in total elapsed time"
    )
    assert attempts >= 1, "it must try at least once"


def test_a_generous_deadline_still_allows_every_attempt():
    """The deadline must not silently reduce resilience on fast failures."""
    resp, attempts, _ = _run_counting(_boom, deadline=60.0, per_attempt=0.0)
    assert resp is None
    assert attempts == 4


def test_success_on_the_first_attempt_is_unaffected():
    resp, attempts, _ = _run_counting(
        lambda request: httpx.Response(200, json={"ok": True}), deadline=0.25)
    assert resp is not None and resp.status_code == 200
    assert attempts == 1


def test_deadline_defaults_are_sane_and_tunable(monkeypatch):
    monkeypatch.delenv("HEARTVAR_HTTP_RETRY_DEADLINE", raising=False)
    default = _http_retry._retry_deadline()
    assert 2.0 <= default <= 12.0, (
        f"default deadline {default}s should bound the measured ~18 s schedule "
        "without truncating a legitimately slow upstream"
    )
    monkeypatch.setenv("HEARTVAR_HTTP_RETRY_DEADLINE", "3.5")
    assert _http_retry._retry_deadline() == 3.5


def test_deadline_falls_back_on_junk(monkeypatch):
    for raw in ("nonsense", "0", "-1", ""):
        monkeypatch.setenv("HEARTVAR_HTTP_RETRY_DEADLINE", raw)
        assert _http_retry._retry_deadline() == 8.0, f"{raw!r} should fall back"


def _captured_timeout(passed, **kw):
    """Return the timeout object request_with_retry hands to httpx."""
    seen = {}

    class _FakeClient:
        async def request(self, method, url, *, timeout, **_):
            seen["timeout"] = timeout
            return httpx.Response(200, request=httpx.Request(method, url))

    async def no_throttle(_url):
        return None

    real = _http_retry.throttle_host
    _http_retry.throttle_host = no_throttle
    try:
        asyncio.run(request_with_retry(
            _FakeClient(), "GET", "https://rest.ensembl.org/x",
            timeout=passed, name="t", **kw))
    finally:
        _http_retry.throttle_host = real
    return seen["timeout"]


def test_a_scalar_timeout_gets_a_short_connect_budget():
    t = _captured_timeout(30.0)
    assert isinstance(t, httpx.Timeout), (
        "a scalar timeout applies to CONNECT too — it must be normalised"
    )
    assert t.connect <= 5.0, f"connect budget {t.connect}s is still too long"
    assert t.read == 30.0, "the caller's read budget must be preserved"


def test_connect_is_never_longer_than_the_callers_own_timeout():
    """A caller asking for 1 s total must not get a 3 s connect."""
    t = _captured_timeout(1.0)
    assert t.connect <= 1.0


def test_an_explicit_httpx_timeout_is_passed_through_untouched():
    given = httpx.Timeout(9.0, connect=7.0)
    t = _captured_timeout(given)
    assert t is given, "an explicit Timeout means the caller has decided"


def test_connect_timeout_is_env_tunable(monkeypatch):
    monkeypatch.setenv("HEARTVAR_HTTP_CONNECT_TIMEOUT", "1.5")
    assert _http_retry._connect_timeout() == 1.5
    monkeypatch.setenv("HEARTVAR_HTTP_CONNECT_TIMEOUT", "junk")
    assert _http_retry._connect_timeout() == 3.0


if __name__ == "__main__":
    _run_all()
