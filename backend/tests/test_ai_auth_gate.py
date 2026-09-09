"""Tests for the AI-only sign-in gate and the per-user AI quota.

HeartVar is public: the site, every reference-database lookup and the
deterministic (AI-free) ACMG classification must stay reachable with no account.
Sign-in exists for one reason — the AI interpretation spends the deployment's own
Anthropic key.

The failure this file mainly guards against is a gate that exists only in the
browser. An earlier revision defined ``require_auth`` but never attached it to any
endpoint, so ``POST /api/curate/stream`` with ``ai_mode="server"`` spent the key
for anyone who skipped the page — client-side checks are UX, not access control.

Ordering matters as much as outcome: the 401 must be raised BEFORE
``/api/curate/stream`` returns its StreamingResponse, because after that the
status line is already on the wire and an exception inside the SSE generator can
only truncate a 200. The tests assert this by patching ``_check_daily_cap`` — the
step immediately after the gate — into a 418 sentinel, so "reached 418" proves the
request was allowed through and "got 401" proves it was not. No network, no
Anthropic calls, no real curation.
"""
from __future__ import annotations

import pytest
from fastapi import HTTPException
from starlette.testclient import TestClient

import backend.app as app_module
from backend import ratelimit

ACCOUNT = {
    "oid": "11111111-2222-3333-4444-555555555555",
    "preferred_username": "curator@hospital.org.au",
    "name": "Test Curator",
}


@pytest.fixture
def client():
    with TestClient(app_module.app, follow_redirects=False) as c:
        yield c


@pytest.fixture
def gate_on(monkeypatch):
    """Microsoft Entra ID configured → the gate is live. Uses the real predicate, which
    reads the env at call time, rather than stubbing it."""
    monkeypatch.setenv("CLIENT_ID", "test-client-id")
    monkeypatch.setenv("CLIENT_SECRET", "test-client-secret")
    monkeypatch.delenv("HEARTVAR_AI_REQUIRES_SIGNIN", raising=False)
    yield


@pytest.fixture
def gate_off(monkeypatch):
    for var in ("CLIENT_ID", "CLIENT_SECRET",
                "GOOGLE_CLIENT_ID", "GOOGLE_CLIENT_SECRET"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.delenv("HEARTVAR_AI_REQUIRES_SIGNIN", raising=False)
    yield


@pytest.fixture
def signed_in(monkeypatch):
    """Simulate a valid session.

    Patches ``session_account`` rather than minting a signed session cookie: the
    cookie is signed with a per-process random key when SESSION_SECRET_KEY is
    unset, so forging one in-test would assert on Starlette's signing rather than
    on our gate. Whether ``session_account`` reads the session correctly is
    covered separately by test_session_account_* below.
    """
    monkeypatch.setattr(app_module, "session_account", lambda request: dict(ACCOUNT))
    yield


@pytest.fixture(autouse=True)
def _fresh_quota(monkeypatch):
    monkeypatch.setattr(ratelimit, "_ai_quota", {"day": None, "counts": {}})
    yield


@pytest.fixture(autouse=True)
def _no_inbound_rate_limit(monkeypatch):
    """Exempt this file from the inbound per-IP limiter.

    Every test in the suite shares one limiter bucket, because TestClient always
    presents as "testclient" — so the 20/minute ceiling on /api/curate/stream is a
    budget shared across all files. This file makes ~16 calls to that endpoint,
    which is enough to exhaust it and make LATER files fail with 429s that look
    like unrelated breakage (it did exactly that on first run).

    These tests are about the auth gate, not rate limiting, so opt out rather than
    spend the shared allowance. Anything that genuinely exercises limiter behaviour
    should leave it enabled.
    """
    monkeypatch.setattr(ratelimit.limiter, "enabled", False)
    yield


@pytest.fixture
def past_gate(monkeypatch):
    """Turn the step *after* the gate into a 418 sentinel, so a test can prove the
    request was allowed through without running a real curation."""
    def _teapot() -> None:
        raise HTTPException(418, "past the gate")
    monkeypatch.setattr(app_module, "_check_daily_cap", _teapot)
    yield


def _curate(ai_mode: str) -> dict:
    return {"gene": "MYH7", "hgvs_c": "c.1208G>A", "ai_mode": ai_mode}


def test_anonymous_evidence_only_curation_is_allowed(client, gate_on, past_gate):
    """ai_mode="none" is the whole tool for an anonymous visitor: every database
    lookup plus the deterministic ACMG classification. 418 = past the gate."""
    assert client.post("/api/curate/stream", json=_curate("none")).status_code == 418


def test_health_and_lookups_need_no_session(client, gate_on):
    assert client.get("/health").status_code == 200
    assert client.get("/api/gene/validate?symbol=MYH7").status_code == 200


def test_landing_page_is_public(client, gate_on):
    assert client.get("/").status_code == 200


def test_data_status_is_public(client, gate_on):
    assert client.get("/api/data-status").status_code == 200


def test_ai_curation_requires_sign_in(client, gate_on, past_gate):
    resp = client.post("/api/curate/stream", json=_curate("server"))
    assert resp.status_code == 401
    assert resp.json()["detail"]["error"] == "signin_required"


def test_ai_curation_401_is_not_a_redirect(client, gate_on, past_gate):
    """/api/curate/stream is called with fetch(), so a 302 to /login would be
    followed transparently and hand the caller the login page's HTML instead of an
    error it can act on."""
    resp = client.post("/api/curate/stream", json=_curate("server"))
    assert resp.status_code == 401
    assert "location" not in {k.lower() for k in resp.headers}


def test_ai_curation_proceeds_when_signed_in(client, gate_on, signed_in, past_gate):
    assert client.post("/api/curate/stream", json=_curate("server")).status_code == 418


def test_chat_requires_sign_in(client, gate_on, monkeypatch):
    """The assistant is pure AI spend with no AI-free equivalent, so it is gated
    unconditionally — unlike curation, there is nothing to degrade to."""
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    resp = client.post("/api/chat", json={"context": "ctx", "question": "why?"})
    assert resp.status_code == 401
    assert resp.json()["detail"]["error"] == "signin_required"


def test_chat_reports_unconfigured_ai_before_asking_for_a_sign_in(client, gate_on,
                                                                 monkeypatch):
    """Deliberate precedence: with no Anthropic key the answer is 503, not 401.
    Asking someone to sign in for a feature that cannot work either way would send
    them through Microsoft to reach the same dead end."""
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", False)
    assert client.post(
        "/api/chat", json={"context": "ctx", "question": "why?"},
    ).status_code == 503


def test_gate_is_dormant_without_client_credentials(client, gate_off, past_gate):
    """Deploying this code before the app registration exists must leave AI
    behaviour exactly as it was, not lock the feature out of reach."""
    assert app_module.ai_auth_required() is False
    assert client.post("/api/curate/stream", json=_curate("server")).status_code == 418


def test_gate_can_be_forced_on(monkeypatch, client, gate_off, past_gate):
    """=1 insists on sign-in even with no credentials, so a misconfiguration fails
    loudly rather than silently serving AI to anonymous callers."""
    monkeypatch.setenv("HEARTVAR_AI_REQUIRES_SIGNIN", "1")
    assert client.post("/api/curate/stream", json=_curate("server")).status_code == 401


def test_gate_can_be_forced_off(monkeypatch, client, gate_on, past_gate):
    monkeypatch.setenv("HEARTVAR_AI_REQUIRES_SIGNIN", "0")
    assert client.post("/api/curate/stream", json=_curate("server")).status_code == 418


def test_account_key_prefers_the_immutable_object_id():
    """preferred_username is an email address and can be reassigned to another
    person; oid cannot."""
    assert app_module.account_key(ACCOUNT).endswith(ACCOUNT["oid"])
    assert app_module.account_key({"preferred_username": "a@b.org"}) == "u:preferred_username:a@b.org"
    assert app_module.account_key({"sub": "xyz"}) == "u:sub:xyz"


def test_account_key_never_returns_empty_for_an_odd_token():
    """An unrecognisable account keys to a shared bucket, not an unlimited one."""
    assert app_module.account_key({}) == "u:unidentified"


def test_account_label_prefers_an_address_then_a_name():
    assert app_module.account_label(ACCOUNT) == "curator@hospital.org.au"
    assert app_module.account_label({"name": "Only Name"}) == "Only Name"
    assert app_module.account_label({}) == "signed-in user"


def test_session_account_is_none_without_session_middleware():
    """Auth not being configured must degrade to "anonymous", never to a 500.
    request.session raises AssertionError when SessionMiddleware is absent."""
    class _NoSession:
        @property
        def session(self):
            raise AssertionError("SessionMiddleware must be installed")

    assert app_module.session_account(_NoSession()) is None


def test_session_account_reads_the_session():
    class _Req:
        session = {"account": dict(ACCOUNT)}

    assert app_module.session_account(_Req())["oid"] == ACCOUNT["oid"]
    assert app_module.session_account(type("R", (), {"session": {}})()) is None


def test_status_reports_the_gate_and_no_identity_when_anonymous(client, gate_on):
    body = client.get("/api/auth/status").json()
    assert body["ai_requires_signin"] is True
    assert body["authenticated"] is False
    assert body["label"] == ""
    assert [p["key"] for p in body["providers"]] == ["microsoft"]


def test_status_reports_the_signed_in_label(client, gate_on, signed_in):
    body = client.get("/api/auth/status").json()
    assert body["authenticated"] is True
    assert body["label"] == "curator@hospital.org.au"


def test_status_is_false_when_unconfigured(client, gate_off):
    """The frontend keys off this to render NO sign-in affordance — /login would
    raise a KeyError without CLIENT_ID."""
    assert client.get("/api/auth/status").json()["ai_requires_signin"] is False


def test_only_configured_providers_are_offered(client, gate_on, monkeypatch):
    """The chooser renders exactly this list, so an unconfigured provider must not
    appear — its login route raises KeyError on the missing settings."""
    monkeypatch.delenv("GOOGLE_CLIENT_ID", raising=False)
    monkeypatch.delenv("GOOGLE_CLIENT_SECRET", raising=False)
    assert client.get("/api/auth/status").json()["providers"] == [
        {"key": "microsoft", "label": "Microsoft"},
    ]


def test_google_appears_once_configured(client, gate_on, monkeypatch):
    monkeypatch.setenv("GOOGLE_CLIENT_ID", "g-client")
    monkeypatch.setenv("GOOGLE_CLIENT_SECRET", "g-secret")
    keys = [p["key"] for p in client.get("/api/auth/status").json()["providers"]]
    assert keys == ["microsoft", "google"]


def test_google_alone_still_enables_the_gate(client, gate_off, monkeypatch, past_gate):
    """A Google-only deployment must gate AI too — the predicate is "any provider",
    not "Microsoft"."""
    monkeypatch.setenv("GOOGLE_CLIENT_ID", "g-client")
    monkeypatch.setenv("GOOGLE_CLIENT_SECRET", "g-secret")
    assert app_module.ai_auth_required() is True
    assert client.post("/api/curate/stream", json=_curate("server")).status_code == 401


def test_google_login_404s_when_unconfigured(client, gate_off):
    """Rather than a 500 from os.environ["GOOGLE_CLIENT_ID"]."""
    assert client.get("/login/google").status_code == 404


def test_google_id_token_validation_rejects_each_failure(monkeypatch):
    import base64 as _b64
    import json as _json
    import time as _time
    monkeypatch.setenv("GOOGLE_CLIENT_ID", "g-client")

    def token(**over):
        claims = {"iss": "https://accounts.google.com", "aud": "g-client",
                  "sub": "123", "nonce": "n1", "exp": _time.time() + 300,
                  "email": "a@b.org"}
        claims.update(over)
        body = _b64.urlsafe_b64encode(_json.dumps(claims).encode()).rstrip(b"=").decode()
        return f"header.{body}.sig"

    assert app_module._claims_from_google_id_token(token(), "n1")["sub"] == "123"
    for bad in ({"iss": "https://evil.example"}, {"aud": "someone-else"},
                {"nonce": "replayed"}, {"exp": _time.time() - 3600}, {"sub": ""}):
        with pytest.raises(ValueError):
            app_module._claims_from_google_id_token(token(**bad), "n1")
    with pytest.raises(ValueError):
        app_module._claims_from_google_id_token("only.two", "n1")


def test_google_account_keys_on_sub():
    """Google issues no `oid`, so the key falls through to `sub` — still stable and
    still per-user, which is what the quota needs."""
    assert app_module.account_key(
        {"sub": "g-12345", "preferred_username": "a@b.org", "idp": "google"}
    ) == "u:sub:g-12345"


def test_microsoft_login_requests_no_graph_scopes(client, gate_on, monkeypatch):
    """The authorize request must ask for no Microsoft Graph permission.

    `scopes=["User.Read"]` bought a Graph access token that nothing ever read —
    the callback stores ID token claims only. Requesting any resource scope can
    require administrator consent, which replaces the login page with a consent
    prompt and blocks sign-in. Deleting the permission from the app registration
    does not help: what is evaluated is the scope we *request* here.

    Asserts on the call rather than on the redirect URL because constructing a real
    MSAL client performs instance discovery over the network.
    """
    seen: dict = {}

    class _FakeMsal:
        def initiate_auth_code_flow(self, scopes, redirect_uri):
            seen["scopes"] = list(scopes)
            return {"auth_uri": "https://login.microsoftonline.com/x/authorize?a=1"}

    monkeypatch.setattr(app_module, "_get_msal_app", lambda: _FakeMsal())
    monkeypatch.setenv("REDIRECT_URI", "https://heartvar.example/auth/callback")

    resp = client.get("/login?next=/")

    assert resp.status_code in (302, 307)
    assert seen["scopes"] == [], (
        "requesting any resource scope can require consent and block sign-in"
    )


def test_quota_blocks_the_user_after_the_cap(monkeypatch):
    monkeypatch.setitem(ratelimit._AI_QUOTA_CAPS, "curation", 3)
    for _ in range(3):
        ratelimit.check_ai_user_quota("u:oid:abc", "curation")
    with pytest.raises(HTTPException) as excinfo:
        ratelimit.check_ai_user_quota("u:oid:abc", "curation")
    assert excinfo.value.status_code == 429
    assert excinfo.value.detail["error"] == "ai_quota_exhausted"


def test_quota_is_per_user_and_per_kind(monkeypatch):
    monkeypatch.setitem(ratelimit._AI_QUOTA_CAPS, "curation", 1)
    monkeypatch.setitem(ratelimit._AI_QUOTA_CAPS, "chat", 1)
    ratelimit.check_ai_user_quota("u:oid:a", "curation")
    ratelimit.check_ai_user_quota("u:oid:b", "curation")
    ratelimit.check_ai_user_quota("u:oid:a", "chat")
    with pytest.raises(HTTPException):
        ratelimit.check_ai_user_quota("u:oid:a", "curation")


def test_quota_can_be_disabled(monkeypatch):
    monkeypatch.setitem(ratelimit._AI_QUOTA_CAPS, "curation", 0)
    for _ in range(40):
        ratelimit.check_ai_user_quota("u:oid:a", "curation")


def test_quota_enforced_through_the_endpoint(client, gate_on, signed_in, past_gate,
                                             monkeypatch):
    monkeypatch.setitem(ratelimit._AI_QUOTA_CAPS, "curation", 2)
    assert client.post("/api/curate/stream", json=_curate("server")).status_code == 418
    assert client.post("/api/curate/stream", json=_curate("server")).status_code == 418
    third = client.post("/api/curate/stream", json=_curate("server"))
    assert third.status_code == 429
    assert third.json()["detail"]["error"] == "ai_quota_exhausted"
    assert client.post("/api/curate/stream", json=_curate("none")).status_code == 418


def test_evidence_only_curations_do_not_consume_quota(client, gate_on, signed_in,
                                                      past_gate, monkeypatch):
    monkeypatch.setitem(ratelimit._AI_QUOTA_CAPS, "curation", 1)
    for _ in range(5):
        assert client.post("/api/curate/stream", json=_curate("none")).status_code == 418
    assert client.post("/api/curate/stream", json=_curate("server")).status_code == 418


def test_rate_limit_key_prefers_the_signed_in_user(monkeypatch):
    """A shared hospital egress otherwise throttles colleagues as a group, and an
    IP bucket is trivially rotated."""
    class _Req:
        session = {"account": dict(ACCOUNT)}
        client = type("c", (), {"host": "203.0.113.9"})()
        headers: dict = {}

    assert ratelimit._rate_limit_key(_Req()).endswith(ACCOUNT["oid"])


def test_rate_limit_key_falls_back_to_ip_when_anonymous(monkeypatch):
    monkeypatch.setattr(ratelimit, "_client_ip", lambda request: "203.0.113.9")

    class _Req:
        session: dict = {}
        headers: dict = {}

    assert ratelimit._rate_limit_key(_Req()) == "203.0.113.9"


def test_rate_limit_key_survives_missing_session_middleware(monkeypatch):
    monkeypatch.setattr(ratelimit, "_client_ip", lambda request: "198.51.100.7")

    class _Req:
        headers: dict = {}

        @property
        def session(self):
            raise AssertionError("SessionMiddleware must be installed")

    assert ratelimit._rate_limit_key(_Req()) == "198.51.100.7"


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
