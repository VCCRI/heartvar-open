"""Tests for POST /api/curate/interpret — adding AI to a finished curation.

The friction: forget to tick "Include AI interpretation", and getting it meant
editing the query and re-running the whole annotation — every reference lookup,
VEP call and literature fetch repeated to reach a step that needs none of them.
The AI half of the pipeline depends only on the gathered evidence, so the server
holds that evidence for half an hour (backend/interpret_cache.py) and this
endpoint resumes from it.

What these tests pin:

  * **The gather does NOT run again.** That is the entire point; the fake gather
    counts its calls, so a regression fails loudly rather than silently costing
    the curator another two minutes of annotation.
  * **One interpretive path.** A post-hoc interpretation goes through the same
    ``_interpret_sse`` as a ticked-box run, so the criteria, points and tier
    must be identical. If they can differ, the button produces a result the
    normal flow would not, and neither can be trusted.
  * **The token IS the gate.** No token is minted for an anonymous run, so the
    button cannot render and an anonymous visitor's only route to AI stays the
    sign-in pill (the user's explicit requirement). The endpoint re-checks
    sign-in anyway — a client-side gate is not a gate.
  * **A token is not a capability transfer.** The held evidence carries the
    curator's clinical context (HPO terms, family/segregation fields), so
    another account presenting the token gets the same 410 as an expired one —
    the response never reveals that the run exists.

No network and no Anthropic calls: the evidence gather and ``stream_claude`` are
both replaced with canned equivalents.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
from starlette.testclient import TestClient

import backend.app as app_module
from backend import interpret_cache, ratelimit

ACCOUNT = {
    "oid": "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee",
    "preferred_username": "curator@hospital.org.au",
    "name": "Test Curator",
}
OTHER_ACCOUNT = {
    "oid": "99999999-8888-7777-6666-555555555555",
    "preferred_username": "someone.else@hospital.org.au",
    "name": "Other Curator",
}

FAKE_MODEL_JSON = json.dumps({
    "classification": "ignored — the server recomputes this",
    "confidence": "moderate",
    "summary": "Canned interpretation for the test.",
    "criteria": [
        {
            "code": "PP2",
            "status": "met",
            "criteria_strength": "PP2_Supporting",
            "direction": "pathogenic",
            "evidence": "Missense-constrained gene.",
        },
    ],
    "gene_context": {},
})


@pytest.fixture(autouse=True)
def _clean_state(monkeypatch):
    """Restore the module-level fakes and empty the token store around each test.

    _gather_evidence_sse and stream_claude are reassigned on the module (the
    pattern test_no_ai_endpoint.py established), so a leak would silently break
    unrelated files.
    """
    original_gather = app_module._gather_evidence_sse
    original_stream = app_module.stream_claude
    interpret_cache.clear()
    monkeypatch.setattr(ratelimit, "_ai_quota", {"day": None, "counts": {}})
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    yield
    app_module._gather_evidence_sse = original_gather
    app_module.stream_claude = original_stream
    interpret_cache.clear()


@pytest.fixture(autouse=True)
def _no_inbound_rate_limit(monkeypatch):
    """Exempt this file from the inbound per-IP limiter.

    TestClient always presents as "testclient", so the 20/minute ceiling on
    /api/curate/stream is one budget shared across the whole suite. This file
    makes ~15 calls to it and would otherwise exhaust the allowance and make a
    LATER file fail with a 429 that looks like unrelated breakage (it did
    exactly that on first run). Same reasoning as test_ai_auth_gate.py.
    """
    monkeypatch.setattr(ratelimit.limiter, "enabled", False)
    yield


@pytest.fixture
def signed_in(monkeypatch):
    monkeypatch.setattr(app_module, "session_account", lambda request: dict(ACCOUNT))
    yield


@pytest.fixture
def anonymous(monkeypatch):
    """No session, and the sign-in gate live (Microsoft Entra ID configured)."""
    monkeypatch.setattr(app_module, "session_account", lambda request: None)
    monkeypatch.setenv("CLIENT_ID", "test-client-id")
    monkeypatch.setenv("CLIENT_SECRET", "test-client-secret")
    monkeypatch.delenv("HEARTVAR_AI_REQUIRES_SIGNIN", raising=False)
    yield


@pytest.fixture
def auth_unconfigured(monkeypatch):
    """No auth providers → ai_auth_required() is False and AI is already open."""
    monkeypatch.setattr(app_module, "session_account", lambda request: None)
    for var in ("CLIENT_ID", "CLIENT_SECRET",
                "GOOGLE_CLIENT_ID", "GOOGLE_CLIENT_SECRET"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.delenv("HEARTVAR_AI_REQUIRES_SIGNIN", raising=False)
    yield


class GatherSpy:
    """A canned evidence gather that records how many times it ran."""

    def __init__(self, evidence: dict | None = None, gene: str = "MYH7"):
        self.calls = 0
        self.evidence = evidence if evidence is not None else {}
        self.gene = gene

    def install(self):
        spy = self

        async def _fake_gather(req, state):
            spy.calls += 1
            state["evidence"] = dict(spy.evidence)
            state["variant_id"] = "14-23424115-G-A"
            state["gene"] = req.gene or spy.gene
            state["hgvs_c"] = req.hgvs_c
            state["timing"] = {}
            return
            yield

        app_module._gather_evidence_sse = _fake_gather
        return self


def install_fake_claude(text: str = FAKE_MODEL_JSON):
    """Replace the Claude stream with one canned chunk."""
    async def _fake_stream(prompt, max_tokens=None):
        yield text
    app_module.stream_claude = _fake_stream


def parse_sse(body: str) -> list[tuple[str, dict]]:
    events = []
    for block in body.split("\n\n"):
        event = data = None
        for line in block.splitlines():
            if line.startswith("event: "):
                event = line[len("event: "):].strip()
            elif line.startswith("data: "):
                data = line[len("data: "):]
        if event and data is not None:
            try:
                events.append((event, json.loads(data)))
            except json.JSONDecodeError:
                events.append((event, {}))
    return events


def only(events, kind) -> dict:
    payloads = [d for e, d in events if e == kind]
    assert len(payloads) == 1, [e for e, _ in events]
    return payloads[0]


def run_evidence_only(client, **overrides) -> dict:
    body = {"gene": "MYH7", "hgvs_c": "c.1988G>A", "ai_mode": "none"}
    body.update(overrides)
    resp = client.post("/api/curate/stream", json=body)
    assert resp.status_code == 200, resp.text
    return only(parse_sse(resp.text), "db_only_complete")


@pytest.fixture
def client():
    return TestClient(app_module.app)


def test_token_issued_for_a_signed_in_evidence_only_run(client, signed_in):
    GatherSpy().install()
    payload = run_evidence_only(client)
    assert payload["preliminary"] is True
    assert isinstance(payload["interpret_token"], str)
    assert len(payload["interpret_token"]) >= 16


def test_no_token_for_an_anonymous_run(client, anonymous):
    """The user's requirement: nothing offers AI after the fact to a visitor who
    is not signed in — the sign-in pill stays the only route."""
    GatherSpy().install()
    assert run_evidence_only(client)["interpret_token"] is None


def test_no_token_when_the_server_has_no_ai_key(client, signed_in, monkeypatch):
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", False)
    GatherSpy().install()
    assert run_evidence_only(client)["interpret_token"] is None


def test_token_issued_when_auth_is_unconfigured(client, auth_unconfigured):
    """Local / non-Azure deployment: AI needs no sign-in there, so withholding
    the token would make the feature unreachable rather than safe."""
    GatherSpy().install()
    assert isinstance(run_evidence_only(client)["interpret_token"], str)


def test_a_ticked_box_run_issues_no_token(client, signed_in):
    """AI already ran — there is nothing to add, so nothing is held."""
    GatherSpy().install()
    install_fake_claude()
    resp = client.post("/api/curate/stream", json={
        "gene": "MYH7", "hgvs_c": "c.1988G>A", "ai_mode": "server",
    })
    assert resp.status_code == 200, resp.text
    kinds = [e for e, _ in parse_sse(resp.text)]
    assert "stage2_complete" in kinds
    assert "db_only_complete" not in kinds
    assert interpret_cache.size() == 0


def test_interpret_resumes_without_re_gathering(client, signed_in):
    """The whole point of the feature."""
    spy = GatherSpy().install()
    token = run_evidence_only(client)["interpret_token"]
    assert spy.calls == 1

    install_fake_claude()
    resp = client.post("/api/curate/interpret", json={"token": token})
    assert resp.status_code == 200, resp.text
    stage2 = only(parse_sse(resp.text), "stage2_complete")
    assert stage2["classification"]
    assert stage2["summary"] == "Canned interpretation for the test."
    assert spy.calls == 1, "the evidence gather ran a second time"


def test_interpret_matches_a_ticked_box_run(client, signed_in):
    """Parity. A post-hoc interpretation must not be able to produce a result
    the normal flow would not — same generator, same criteria, same tier."""
    GatherSpy().install()
    install_fake_claude()

    direct = only(parse_sse(client.post("/api/curate/stream", json={
        "gene": "MYH7", "hgvs_c": "c.1988G>A", "ai_mode": "server",
    }).text), "stage2_complete")

    token = run_evidence_only(client)["interpret_token"]
    posthoc = only(parse_sse(
        client.post("/api/curate/interpret", json={"token": token}).text,
    ), "stage2_complete")

    assert posthoc["classification"] == direct["classification"]
    assert posthoc["points_total"] == direct["points_total"]
    assert posthoc["criteria"] == direct["criteria"]
    assert posthoc["variant_id"] == direct["variant_id"]


def test_interpret_can_be_retried_within_the_window(client, signed_in):
    """A dead stream or an unparseable response should be retryable, so the
    token is NOT burned on use."""
    GatherSpy().install()
    token = run_evidence_only(client)["interpret_token"]
    install_fake_claude()
    for _ in range(2):
        resp = client.post("/api/curate/interpret", json={"token": token})
        assert resp.status_code == 200, resp.text


def test_evidence_is_never_accepted_from_the_client(client, signed_in):
    """The request model carries a token and nothing else, so a caller cannot
    hand the prompt forged evidence."""
    GatherSpy().install()
    token = run_evidence_only(client)["interpret_token"]
    install_fake_claude()
    resp = client.post("/api/curate/interpret", json={
        "token": token,
        "db_evidence": {"gnomad": {"popmax_af": 0.5}},
    })
    assert resp.status_code == 200, resp.text
    assert only(parse_sse(resp.text), "stage2_complete")["classification"] != "Benign"


def test_unknown_token_is_410(client, signed_in):
    resp = client.post("/api/curate/interpret", json={"token": "n" * 43})
    assert resp.status_code == 410
    assert resp.json()["detail"]["error"] == "interpret_expired"


def test_expired_token_is_410(client, signed_in, monkeypatch):
    GatherSpy().install()
    token = run_evidence_only(client)["interpret_token"]
    monkeypatch.setattr(interpret_cache, "TTL_SECONDS", -1)
    assert client.post(
        "/api/curate/interpret", json={"token": token},
    ).status_code == 410


def test_another_accounts_token_is_410(client, signed_in, monkeypatch):
    """Indistinguishable from expired, deliberately: the held evidence contains
    the first curator's phenotype and family data, and the response must not
    confirm that their run exists."""
    GatherSpy().install()
    token = run_evidence_only(client)["interpret_token"]
    monkeypatch.setattr(
        app_module, "session_account", lambda request: dict(OTHER_ACCOUNT),
    )
    resp = client.post("/api/curate/interpret", json={"token": token})
    assert resp.status_code == 410
    assert resp.json()["detail"]["error"] == "interpret_expired"


def test_interpret_requires_sign_in(client, anonymous):
    """A client-side gate is UX, not access control — the endpoint spends the
    owner's key, so it enforces the session itself."""
    resp = client.post("/api/curate/interpret", json={"token": "n" * 43})
    assert resp.status_code == 401
    assert resp.json()["detail"]["error"] == "signin_required"


def test_interpret_refuses_when_the_daily_budget_is_spent(
    client, signed_in, monkeypatch,
):
    GatherSpy().install()
    token = run_evidence_only(client)["interpret_token"]
    monkeypatch.setattr(app_module.budget, "ai_within_budget", lambda: False)
    resp = client.post("/api/curate/interpret", json={"token": token})
    assert resp.status_code == 503
    assert resp.json()["detail"]["error"] == "budget"


def test_interpret_refuses_when_ai_is_unconfigured(client, signed_in, monkeypatch):
    GatherSpy().install()
    token = run_evidence_only(client)["interpret_token"]
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", False)
    resp = client.post("/api/curate/interpret", json={"token": token})
    assert resp.status_code == 503
    assert resp.json()["detail"]["error"] == "unconfigured"


def test_interpret_consumes_the_per_user_ai_quota(client, signed_in, monkeypatch):
    """It spends the key exactly like a ticked-box curation, so it counts."""
    seen = []
    monkeypatch.setattr(
        app_module, "check_ai_user_quota",
        lambda key, kind: seen.append((key, kind)),
    )
    GatherSpy().install()
    token = run_evidence_only(client)["interpret_token"]
    install_fake_claude()
    assert client.post(
        "/api/curate/interpret", json={"token": token},
    ).status_code == 200
    assert [kind for _, kind in seen] == ["curation"]


def test_a_dead_token_costs_no_quota(client, signed_in, monkeypatch):
    seen = []
    monkeypatch.setattr(
        app_module, "check_ai_user_quota",
        lambda key, kind: seen.append((key, kind)),
    )
    assert client.post(
        "/api/curate/interpret", json={"token": "n" * 43},
    ).status_code == 410
    assert seen == []


def test_a_short_token_is_rejected_by_the_model(client, signed_in):
    assert client.post(
        "/api/curate/interpret", json={"token": "abc"},
    ).status_code == 422


def test_store_evicts_the_oldest_beyond_the_cap():
    tokens = [
        interpret_cache.store({"n": i}, "acct")
        for i in range(interpret_cache.MAX_ENTRIES + 3)
    ]
    assert interpret_cache.size() == interpret_cache.MAX_ENTRIES
    assert interpret_cache.get(tokens[0], "acct") is None
    assert interpret_cache.get(tokens[-1], "acct") == {"n": len(tokens) - 1}


def test_store_expires_entries(monkeypatch):
    token = interpret_cache.store({"n": 1}, "acct")
    assert interpret_cache.get(token, "acct") == {"n": 1}
    monkeypatch.setattr(interpret_cache, "TTL_SECONDS", -1)
    assert interpret_cache.get(token, "acct") is None
    assert interpret_cache.size() == 0


def test_store_is_account_scoped():
    token = interpret_cache.store({"n": 1}, "acct-a")
    assert interpret_cache.get(token, "acct-b") is None
    assert interpret_cache.get(token, None) is None
    assert interpret_cache.get(token, "acct-a") == {"n": 1}


_JS_PATH = Path(__file__).resolve().parent.parent.parent / "static" / "heartvar.js"


def test_js_calls_the_interpret_endpoint():
    assert "/api/curate/interpret" in _JS_PATH.read_text()


def test_js_gates_the_button_on_the_token():
    """No token → no button. The token is minted server-side for signed-in
    curators only, so this guard is what keeps the post-hoc AI route off an
    anonymous visitor's screen."""
    src = _JS_PATH.read_text()
    guard = src.index("const _tok = payload.interpret_token")
    creation = src.index("btn.id = 'prelim-ai-btn'")
    assert guard < creation, (
        "the Add AI interpretation button is created before the "
        "interpret_token check — an anonymous run would render it"
    )
    assert "if (_tok) {" in src[guard:creation]


def test_js_sends_only_the_token():
    """Evidence must never round-trip through the browser: a client that can
    supply evidence can shape the prompt and the ACMG score."""
    src = _JS_PATH.read_text()
    assert "JSON.stringify({ token: HV_INTERPRET.token })" in src


def test_js_has_one_stage2_shaping_path():
    """The ticked-box stream and the post-hoc button must hand the renderer the
    same object. A second inline copy of the shaping would drift the moment a
    field is added to stage2_complete."""
    src = _JS_PATH.read_text()
    assert src.count("borderline_reasoning: (data && data.borderline_reasoning)") == 1
    assert src.count("function _shapeAiResult") == 1
