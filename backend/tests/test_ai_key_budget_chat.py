"""Tests for the server-side AI key + cost-budget + chatbot rewiring.

HeartVar moved from bring-your-own-key to a single OWNER Anthropic key with
abuse/cost controls. This file covers the new surfaces:

  1. backend/budget.py — per-model cost estimation (incl. cache-write/read
     multipliers + the unknown-model Opus fallback) and the per-UTC-day spend
     ledger that flips ``ai_within_budget()`` False once the daily USD cap is
     exceeded (and stays True when the cap is disabled with budget 0).
  2. POST /api/chat — the server-side "ask about this variant" assistant:
     {answer} on success, 503 when AI is unconfigured, 429 when the daily
     budget is exhausted, 422 on an empty question.
  3. The graceful curate-degrade path — ai_mode="server" but AI unavailable
     (budget exhausted / unconfigured) must emit a single ``db_only_complete``
     SSE event carrying the right ``ai_unavailable`` reason and NO
     ``stage2_complete``.
  4. ``_load_anthropic_key`` precedence — a present API.txt is authoritative;
     otherwise the ANTHROPIC_API_KEY env var is preserved.

No real network / Anthropic calls: ``backend.app.chat_reply`` is patched, the
evidence gather is stubbed exactly like test_no_ai_endpoint.py, and the budget
ledger is reset between tests so they are order-independent.

Runnable with pytest or directly
(``python -m backend.tests.test_ai_key_budget_chat``) — the __main__ runner
builds the same monkeypatch shim the fixtures rely on.
"""
from __future__ import annotations

import json

import pytest
from starlette.testclient import TestClient

import backend.app as app_module
from backend import budget


@pytest.fixture(autouse=True)
def _fresh_budget_ledger(monkeypatch):
    """Reset the process-wide spend ledger before every test so per-day tallies
    never leak across tests (the ledger is a module singleton)."""
    monkeypatch.setattr(budget, "_LEDGER", budget._DailyLedger())
    yield


@pytest.fixture(autouse=True)
def _restore_gather():
    """The curate tests reassign the module-level _gather_evidence_sse;
    snapshot + restore it so the fake never leaks to other test files."""
    original = app_module._gather_evidence_sse
    yield
    app_module._gather_evidence_sse = original


def _install_fake_gather(evidence: dict | None = None):
    """Replace _gather_evidence_sse with an async generator that populates
    `state` from a canned evidence dict and emits no SSE events. Mirrors
    test_no_ai_endpoint.py so the curate path is fully offline."""
    evidence = evidence or {}

    async def _fake_gather(req, state):
        state["evidence"] = dict(evidence)
        state["variant_id"] = "1-1000-A-T"
        state["gene"] = req.gene or "MYH7"
        state["hgvs_c"] = req.hgvs_c
        state["timing"] = {}
        return
        yield

    app_module._gather_evidence_sse = _fake_gather


def _parse_sse(text: str) -> list[tuple[str, dict]]:
    events = []
    for block in text.split("\n\n"):
        event = None
        data = None
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


class _Usage:
    """Minimal stand-in for the Anthropic SDK usage object (record_usage only
    reads these attributes via getattr)."""

    def __init__(self, input_tokens=0, output_tokens=0,
                 cache_creation_input_tokens=0, cache_read_input_tokens=0):
        self.input_tokens = input_tokens
        self.output_tokens = output_tokens
        self.cache_creation_input_tokens = cache_creation_input_tokens
        self.cache_read_input_tokens = cache_read_input_tokens


def test_estimate_cost_sonnet_plain_in_out():
    cost = budget.estimate_cost_usd(
        "claude-sonnet-4-20250514", input_tokens=1_000_000, output_tokens=1_000_000,
    )
    assert cost == pytest.approx(3.0 + 15.0)


def test_estimate_cost_haiku_plain_in_out():
    cost = budget.estimate_cost_usd(
        "claude-haiku-4-5", input_tokens=1_000_000, output_tokens=1_000_000,
    )
    assert cost == pytest.approx(1.0 + 5.0)


def test_estimate_cost_cache_write_is_1_25x_input():
    cost = budget.estimate_cost_usd(
        "claude-sonnet-4", cache_creation_input_tokens=1_000_000,
    )
    assert cost == pytest.approx(3.0 * 1.25)


def test_estimate_cost_cache_read_is_0_10x_input():
    cost = budget.estimate_cost_usd(
        "claude-sonnet-4", cache_read_input_tokens=1_000_000,
    )
    assert cost == pytest.approx(3.0 * 0.10)


def test_estimate_cost_all_token_buckets_combine():
    cost = budget.estimate_cost_usd(
        "claude-haiku-4-5",
        input_tokens=1_000_000,
        cache_creation_input_tokens=1_000_000,
        cache_read_input_tokens=1_000_000,
        output_tokens=1_000_000,
    )
    expected = 1.0 + (1.0 * 1.25) + (1.0 * 0.10) + 5.0
    assert cost == pytest.approx(expected)


def test_estimate_cost_unknown_model_falls_back_to_opus():
    cost = budget.estimate_cost_usd(
        "some-future-model-x9", input_tokens=1_000_000, output_tokens=1_000_000,
    )
    assert cost == pytest.approx(5.0 + 25.0)


def test_record_usage_accumulates_per_day():
    budget.record_usage("claude-haiku-4-5", _Usage(input_tokens=1_000_000))
    budget.record_usage("claude-haiku-4-5", _Usage(output_tokens=1_000_000))
    assert budget._LEDGER.cost_usd == pytest.approx(6.0)
    assert budget._LEDGER.calls == 2


def test_ai_within_budget_flips_false_when_exceeded(monkeypatch):
    monkeypatch.setenv("HEARTVAR_DAILY_USD_BUDGET", "3.0")
    assert budget.ai_within_budget() is True
    budget.record_usage("claude-haiku-4-5", _Usage(output_tokens=1_000_000))
    assert budget.ai_within_budget() is False


def test_ai_within_budget_true_when_budget_disabled(monkeypatch):
    monkeypatch.setenv("HEARTVAR_DAILY_USD_BUDGET", "0")
    budget.record_usage("claude-opus-4", _Usage(output_tokens=10_000_000))
    assert budget.ai_within_budget() is True


def test_record_usage_never_raises_on_bad_usage():
    cost = budget.record_usage("claude-haiku-4-5", object())
    assert cost == 0.0


def test_reserve_blocks_concurrent_burst(monkeypatch):
    monkeypatch.setenv("HEARTVAR_DAILY_USD_BUDGET", "0.10")
    assert budget.reserve(0.06) is True
    assert budget._LEDGER.reserved == pytest.approx(0.06)
    assert budget.ai_within_budget() is True
    assert budget.reserve(0.06) is False
    assert budget._LEDGER.reserved == pytest.approx(0.06)


def test_release_frees_the_reservation(monkeypatch):
    monkeypatch.setenv("HEARTVAR_DAILY_USD_BUDGET", "0.10")
    assert budget.reserve(0.06) is True
    budget.release(0.06)
    assert budget._LEDGER.reserved == pytest.approx(0.0)
    assert budget.reserve(0.06) is True


def test_reserve_always_admits_when_budget_disabled(monkeypatch):
    monkeypatch.setenv("HEARTVAR_DAILY_USD_BUDGET", "0")
    assert budget.reserve(9999.0) is True
    assert budget._LEDGER.reserved == pytest.approx(0.0)


def test_pricing_override_accepts_both_key_spellings(monkeypatch):
    monkeypatch.setenv(
        "HEARTVAR_MODEL_PRICING",
        json.dumps({"claude-sonnet-4-6": {"input": 9.0, "output": 40.0}}),
    )
    pricing = budget._load_pricing()
    assert pricing["claude-sonnet-4-6"] == {"in": 9.0, "out": 40.0}


def test_chat_returns_answer_on_success(monkeypatch):
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module.budget, "ai_within_budget", lambda: True)

    async def _fake_reply(context, question, history, mode="chat"):
        return "MYH7 R403Q is a well-established HCM variant."

    monkeypatch.setattr(app_module, "chat_reply", _fake_reply)
    client = TestClient(app_module.app)
    resp = client.post("/api/chat", json={
        "context": "MYH7 c.1208G>A",
        "question": "Is this pathogenic?",
        "history": [],
    })
    assert resp.status_code == 200, resp.text
    assert resp.json() == {"answer": "MYH7 R403Q is a well-established HCM variant."}


def test_chat_passes_history_through(monkeypatch):
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module.budget, "ai_within_budget", lambda: True)
    seen = {}

    async def _fake_reply(context, question, history, mode="chat"):
        seen["history"] = history
        return "ok"

    monkeypatch.setattr(app_module, "chat_reply", _fake_reply)
    client = TestClient(app_module.app)
    resp = client.post("/api/chat", json={
        "context": "ctx",
        "question": "follow-up?",
        "history": [{"role": "u", "text": "first"}, {"role": "a", "text": "reply"}],
    })
    assert resp.status_code == 200, resp.text
    assert seen["history"] == [
        {"role": "u", "text": "first"}, {"role": "a", "text": "reply"},
    ]


def test_chat_503_when_ai_unconfigured(monkeypatch):
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", False)

    async def _should_not_run(context, question, history):  # pragma: no cover
        raise AssertionError("chat_reply must not be called when AI unconfigured")

    monkeypatch.setattr(app_module, "chat_reply", _should_not_run)
    client = TestClient(app_module.app)
    resp = client.post("/api/chat", json={"context": "", "question": "hi"})
    assert resp.status_code == 503, resp.text


def test_chat_429_when_budget_exhausted(monkeypatch):
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module.budget, "reserve", lambda est: False)

    async def _should_not_run(context, question, history):  # pragma: no cover
        raise AssertionError("chat_reply must not be called when over budget")

    monkeypatch.setattr(app_module, "chat_reply", _should_not_run)
    client = TestClient(app_module.app)
    resp = client.post("/api/chat", json={"context": "", "question": "hi"})
    assert resp.status_code == 429, resp.text


def test_chat_422_on_empty_question(monkeypatch):
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module.budget, "ai_within_budget", lambda: True)
    client = TestClient(app_module.app)
    resp = client.post("/api/chat", json={"context": "ctx", "question": ""})
    assert resp.status_code == 422, resp.text


def _run_curate(body: dict) -> list[tuple[str, dict]]:
    client = TestClient(app_module.app)
    resp = client.post("/api/curate/stream", json=body)
    assert resp.status_code == 200, resp.text
    return _parse_sse(resp.text)


def test_curate_degrades_to_evidence_only_on_budget(monkeypatch):
    _install_fake_gather({})
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module.budget, "ai_within_budget", lambda: False)
    events = _run_curate({"gene": "MYH7", "hgvs_c": "c.1988G>A", "ai_mode": "server"})
    kinds = [e for e, _ in events]
    assert "stage2_complete" not in kinds, kinds
    payloads = [d for e, d in events if e == "db_only_complete"]
    assert len(payloads) == 1, kinds
    assert payloads[0]["ai_unavailable"] == "budget"
    assert payloads[0]["preliminary"] is True


def test_curate_degrades_to_evidence_only_when_unconfigured(monkeypatch):
    _install_fake_gather({})
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", False)
    events = _run_curate({"gene": "MYH7", "hgvs_c": "c.1988G>A", "ai_mode": "server"})
    kinds = [e for e, _ in events]
    assert "stage2_complete" not in kinds, kinds
    payloads = [d for e, d in events if e == "db_only_complete"]
    assert len(payloads) == 1, kinds
    assert payloads[0]["ai_unavailable"] == "unconfigured"


def test_curate_opt_out_has_null_ai_unavailable(monkeypatch):
    _install_fake_gather({})
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module.budget, "ai_within_budget", lambda: True)
    events = _run_curate({"gene": "MYH7", "hgvs_c": "c.1988G>A", "ai_mode": "none"})
    payloads = [d for e, d in events if e == "db_only_complete"]
    assert len(payloads) == 1
    assert payloads[0]["ai_unavailable"] is None


def test_load_key_prefers_api_txt(monkeypatch, tmp_path):
    (tmp_path / "API.txt").write_text("sk-ant-from-file-xyz\n", encoding="utf-8")
    monkeypatch.setattr(app_module, "PROJECT_ROOT", tmp_path)
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-ant-from-env-should-be-overridden")
    app_module._load_anthropic_key()
    assert app_module.os.environ["ANTHROPIC_API_KEY"] == "sk-ant-from-file-xyz"


def test_load_key_falls_back_to_env_when_no_api_txt(monkeypatch, tmp_path):
    assert not (tmp_path / "API.txt").exists()
    monkeypatch.setattr(app_module, "PROJECT_ROOT", tmp_path)
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-ant-env-only")
    app_module._load_anthropic_key()
    assert app_module.os.environ["ANTHROPIC_API_KEY"] == "sk-ant-env-only"


def _run_all():
    """Best-effort standalone runner. pytest is the supported path (the
    monkeypatch-based fixtures don't run here); this just smoke-checks the
    pure-budget tests so the file is importable + runnable."""
    pure = [
        v for k, v in sorted(globals().items())
        if k.startswith("test_") and "monkeypatch" not in v.__code__.co_varnames
    ]
    budget._LEDGER = budget._DailyLedger()
    for fn in pure:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(pure)} pure-budget tests passed "
          f"(run the full file with pytest for the endpoint tests)")


if __name__ == "__main__":
    _run_all()
