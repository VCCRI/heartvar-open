"""Tests for the chat assistant's "Generate a clinical report" mode.

The report is a different artefact from a chat answer, so it cannot share the
chat preamble: that one caps replies at ~120 words, bans headings and limits
lists to 3-4 bullets, which is the opposite of a 15-bullet structured report.
Rather than loosen those rules for every question (making all answers verbose),
``chat_reply`` takes a ``mode`` and swaps in a report preamble plus a larger
token ceiling.

What matters about a CLINICAL report is that nothing in it is invented, so most
of these tests pin the grounding rules rather than the prose:

  * the report preamble must keep the chat scoping (this variant only, answer
    exclusively from the context) and must instruct OMISSION of any line the
    loaded evidence cannot support. HeartVar has no penetrance or expressivity
    source, so the template's penetrance line would otherwise be filled with a
    PMID the model invented — a fabricated citation in a clinical document is
    the single worst failure available here;
  * identifiers (PMID / MIM / ClinVar accession) may only be echoed, never
    generated;
  * no VUS sub-tier ("VUS-3A") may be asserted — HeartVar computes the 5-tier
    ACMG call and nothing finer, so a sub-tier could only be guessed;
  * normal chat is byte-identical to before, so adding the mode cannot make
    ordinary answers longer or less scoped;
  * the budget RESERVATION scales with the mode's real output cap, or a report
    under-books the daily spend it is about to consume.
"""
from __future__ import annotations

import asyncio

import pytest
from starlette.testclient import TestClient

from backend import claude as claude_module
from backend import app as app_module
from backend.app import app
from backend.models import ChatRequest


def test_mode_defaults_to_chat():
    assert ChatRequest(question="why?").mode == "chat"


@pytest.mark.parametrize("mode", ["chat", "report"])
def test_valid_modes_accepted(mode):
    assert ChatRequest(question="q", mode=mode).mode == mode


@pytest.mark.parametrize("mode", ["REPORT", "essay", "", "chat ", "x" * 40])
def test_invalid_mode_rejected(mode):
    with pytest.raises(Exception):
        ChatRequest(question="q", mode=mode)


def test_endpoint_rejects_an_unknown_mode(monkeypatch):
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module.budget, "ai_within_budget", lambda: True)
    r = TestClient(app).post("/api/chat", json={
        "context": "c", "question": "q", "history": [], "mode": "freeform",
    })
    assert r.status_code == 422


class _FakeUsage:
    input_tokens = 10
    output_tokens = 10


class _FakeBlock:
    type = "text"

    def __init__(self, text):
        self.text = text


class _FakeMsg:
    def __init__(self, text):
        self.content = [_FakeBlock(text)]
        self.usage = _FakeUsage()


def _capture(monkeypatch):
    """Patch the Anthropic client and return the dict of kwargs it was called
    with, so the assembled system prompt and caps can be asserted."""
    seen = {}

    class _Messages:
        async def create(self, **kwargs):
            seen.update(kwargs)
            return _FakeMsg("stub")

    class _Client:
        messages = _Messages()

    monkeypatch.setattr(claude_module, "_get_client", lambda: _Client())
    monkeypatch.setattr(claude_module.budget, "record_usage", lambda *a, **k: None)
    return seen


def test_report_mode_uses_the_report_preamble(monkeypatch):
    seen = _capture(monkeypatch)
    asyncio.run(claude_module.chat_reply("CTX", "Generate a clinical report", mode="report"))
    system_text = seen["system"][0]["text"]
    assert claude_module._REPORT_SYSTEM_PREAMBLE in system_text
    assert claude_module._CHAT_SYSTEM_PREAMBLE not in system_text
    assert "CTX" in system_text


def test_chat_mode_is_unchanged(monkeypatch):
    """Adding the mode must not alter ordinary answers."""
    seen = _capture(monkeypatch)
    asyncio.run(claude_module.chat_reply("CTX", "why?"))
    assert seen["system"][0]["text"] == claude_module._CHAT_SYSTEM_PREAMBLE + "CTX"
    assert seen["max_tokens"] == claude_module.CHAT_MAX_TOKENS


def test_report_mode_raises_the_output_cap(monkeypatch):
    """A 15-bullet report does not fit the 400-token chat ceiling."""
    seen = _capture(monkeypatch)
    asyncio.run(claude_module.chat_reply("CTX", "report", mode="report"))
    assert seen["max_tokens"] == claude_module.REPORT_MAX_TOKENS
    assert claude_module.REPORT_MAX_TOKENS > claude_module.CHAT_MAX_TOKENS


def test_report_mode_allows_a_larger_context(monkeypatch):
    """The report needs facts a chat answer does not (transcript, mechanism,
    conditions, comparable variants), so its context cap is higher — but still
    a hard cap."""
    seen = _capture(monkeypatch)
    big = "x" * 40000
    asyncio.run(claude_module.chat_reply(big, "report", mode="report"))
    body = seen["system"][0]["text"]
    kept = len(body) - len(claude_module._REPORT_SYSTEM_PREAMBLE)
    assert kept == claude_module._REPORT_CONTEXT_MAX_CHARS
    assert claude_module._REPORT_CONTEXT_MAX_CHARS > claude_module._CHAT_CONTEXT_MAX_CHARS


def test_report_preamble_keeps_the_chat_scoping():
    """A report must be no less scoped than a chat answer."""
    p = claude_module._REPORT_SYSTEM_PREAMBLE.lower()
    assert "only" in p and "variant" in p
    for phrase in ("do not invent", "evidence context"):
        assert phrase in p, phrase


def test_report_preamble_requires_omission_of_unbacked_lines():
    """HeartVar has no penetrance/expressivity source. Without this rule the
    model fills those template lines with a PMID it made up."""
    p = claude_module._REPORT_SYSTEM_PREAMBLE.lower()
    assert "omit" in p
    assert "penetrance" in p


def test_report_preamble_forbids_generating_identifiers():
    p = claude_module._REPORT_SYSTEM_PREAMBLE.lower()
    assert "pmid" in p and "mim" in p
    assert "verbatim" in p


def test_report_preamble_forbids_a_vus_subtier():
    """HeartVar computes the 5-tier ACMG call and nothing finer, so 'VUS-3A'
    could only be a guess dressed as a classification."""
    p = claude_module._REPORT_SYSTEM_PREAMBLE
    assert "sub-tier" in p.lower()
    assert "3A" in p


def test_report_preamble_carries_the_exact_section_labels():
    """The curator pastes this into a lab system; the headings are the format."""
    p = claude_module._REPORT_SYSTEM_PREAMBLE
    assert "This variant is classified as" in p
    assert "Evidence in support of pathogenic classification:" in p
    assert "Evidence in support of benign classification:" in p
    assert "Additional information:" in p


def test_endpoint_forwards_the_mode(monkeypatch):
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module.budget, "ai_within_budget", lambda: True)
    seen = {}

    async def _fake_reply(context, question, history, mode="chat"):
        seen["mode"] = mode
        return "report text"

    monkeypatch.setattr(app_module, "chat_reply", _fake_reply)
    r = TestClient(app).post("/api/chat", json={
        "context": "c", "question": "Generate a clinical report",
        "history": [], "mode": "report",
    })
    assert r.status_code == 200, r.text
    assert seen["mode"] == "report"
    assert r.json() == {"answer": "report text"}


def test_report_reserves_the_larger_budget(monkeypatch):
    """The reservation exists so a burst cannot all pass the gate before any of
    them records spend. A report costs more, so it must book more."""
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module.budget, "ai_within_budget", lambda: True)
    booked = []

    def _estimate(model, in_tok, out_tok):
        booked.append(out_tok)
        return 0.0

    monkeypatch.setattr(app_module.budget, "estimate_call_cost", _estimate)
    monkeypatch.setattr(app_module.budget, "reserve", lambda c: True)
    monkeypatch.setattr(app_module.budget, "release", lambda c: None)

    async def _fake_reply(context, question, history, mode="chat"):
        return "ok"

    monkeypatch.setattr(app_module, "chat_reply", _fake_reply)
    client = TestClient(app)
    client.post("/api/chat", json={"context": "c", "question": "q", "history": []})
    client.post("/api/chat", json={
        "context": "c", "question": "q", "history": [], "mode": "report",
    })
    assert booked == [claude_module.CHAT_MAX_TOKENS, claude_module.REPORT_MAX_TOKENS]
