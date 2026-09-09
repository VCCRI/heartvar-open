"""Endpoint-wiring test for the no-key (AI-free) classification path.

`POST /api/curate/stream` with ai_mode="none" must run the deterministic
scoring pipeline over the hard-coded + supplementary criteria and emit a
`db_only_complete` SSE event carrying `classification`, `points_total`,
`preliminary: true`, and the full 28-code `criteria` list — WITHOUT making
any LLM call. The evidence gather (the only network path) is monkeypatched
to a canned state so the test is fully offline and deterministic.

No pytest dependency for the assertions themselves, but the FastAPI
TestClient needs the app importable. Runnable with pytest or directly
(``python -m backend.tests.test_no_ai_endpoint``).
"""
from __future__ import annotations

import json

import pytest
from starlette.testclient import TestClient

import backend.app as app_module


@pytest.fixture(autouse=True)
def _restore_gather():
    """_install_fake_gather reassigns the module-level _gather_evidence_sse;
    snapshot + restore it around each test so the fake never leaks to other
    test files (it previously persisted for the rest of the session)."""
    original = app_module._gather_evidence_sse
    yield
    app_module._gather_evidence_sse = original


def _install_fake_gather(evidence: dict):
    """Replace _gather_evidence_sse with an async generator that populates
    `state` from a canned evidence dict and emits no SSE events."""
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


def _run_no_ai(body: dict) -> dict:
    """POST the no-AI request and return the db_only_complete payload."""
    client = TestClient(app_module.app)
    resp = client.post("/api/curate/stream", json=body)
    assert resp.status_code == 200, resp.text
    events = _parse_sse(resp.text)
    kinds = [e for e, _ in events]
    assert "stage2_complete" not in kinds, kinds
    assert "ai_payload" not in kinds, kinds
    payloads = [d for e, d in events if e == "db_only_complete"]
    assert len(payloads) == 1, kinds
    return payloads[0]


def test_no_ai_endpoint_emits_preliminary_classification():
    _install_fake_gather({})
    payload = _run_no_ai({"gene": "MYH7", "hgvs_c": "c.1988G>A", "ai_mode": "none"})
    assert payload["preliminary"] is True
    assert payload["classification"] == "VUS"
    assert 0 <= payload["points_total"] <= 5
    assert len(payload["criteria"]) == 28


def test_no_ai_endpoint_tags_judgement_codes_not_assessed():
    """The 7 LLM codes, plus PM1/PP1/BS4 on this fixture.

    PM1, PP1 and BS4 became Python-authoritative on 2026-09-08 and left
    AI_EVALUATED_CRITERIA_CODES, but they are derived from curator input
    (PP1/BS4) or domain-P/LP evidence (PM1) and this request supplies neither, so
    they still come back not_assessed. That is the load-bearing part: the no-AI
    arm's criteria strings are IDENTICAL either side of the change, because these
    three always came from Python here. If one of them ever arrives as `not_met`
    instead, the no-AI arm has moved and the change was not the no-op it claims
    to be.
    """
    _install_fake_gather({})
    payload = _run_no_ai({"gene": "MYH7", "hgvs_c": "c.1988G>A", "ai_mode": "none"})
    not_assessed = {c["code"] for c in payload["criteria"] if c["status"] == "not_assessed"}
    assert not_assessed == (
        set(app_module.AI_EVALUATED_CRITERIA_CODES) | {"PM1", "PP1", "BS4"}
    )


def test_no_ai_endpoint_is_the_default_mode():
    _install_fake_gather({})
    payload = _run_no_ai({"gene": "MYH7", "hgvs_c": "c.1988G>A"})
    assert payload["preliminary"] is True
    assert payload["classification"] in {
        "Pathogenic", "Likely pathogenic", "VUS", "Likely benign", "Benign",
    }


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()


_MINIMAL_EVIDENCE = {"chdgene": {"ok": True, "listed": False}}


def _install_gather_emitting_scoring_ready(evidence: dict, waiting_on: list):
    """A gather that emits scoring_ready with a caller-chosen `waiting_on`."""
    async def _fake_gather(req, state):
        state["evidence"] = dict(evidence)
        state["variant_id"] = "1-1000-A-T"
        state["gene"] = req.gene or "MYH7"
        state["hgvs_c"] = req.hgvs_c
        state["timing"] = {}
        yield f'event: scoring_ready\ndata: {json.dumps({"waiting_on": waiting_on})}\n\n'
    app_module._gather_evidence_sse = _fake_gather


def _preliminary_events(waiting_on: list) -> list[str]:
    _install_gather_emitting_scoring_ready(_MINIMAL_EVIDENCE, waiting_on)
    client = TestClient(app_module.app)
    resp = client.post("/api/curate/stream",
                       json={"gene": "MYH7", "hgvs_c": "c.1208G>A",
                             "ai_mode": "none"})
    assert resp.status_code == 200, resp.text
    return [e for e, _ in _parse_sse(resp.text)]


def test_preliminary_fires_when_literature_is_still_outstanding():
    """The original case, kept so the change cannot regress it."""
    kinds = _preliminary_events(["pubtator3", "pmcoa"])
    assert "preliminary_classification" in kinds, kinds


def test_preliminary_still_fires_when_nothing_is_outstanding():
    """The case that silently stopped working. An empty `waiting_on` must still
    produce the early tier — on CHD7 this was the difference between a tier at
    11.8 s and no tier until 15.6 s."""
    kinds = _preliminary_events([])
    assert "preliminary_classification" in kinds, (
        "an empty waiting_on suppressed the early tier — this is the CHD7 "
        f"regression: {kinds}")


def test_scoring_ready_is_never_forwarded_to_the_client():
    """app.py consumes it as a trigger. Leaking it would put an internal
    diagnostic on a public stream."""
    for waiting_on in ([], ["pmcoa"]):
        kinds = _preliminary_events(waiting_on)
        assert "scoring_ready" not in kinds, kinds


def _events_for(evidence: dict, waiting_on: list) -> dict:
    """{event_name: payload} for one evidence-only run."""
    _install_gather_emitting_scoring_ready(evidence, waiting_on)
    client = TestClient(app_module.app)
    resp = client.post("/api/curate/stream",
                       json={"gene": "MYH7", "hgvs_c": "c.1208G>A",
                             "ai_mode": "none"})
    assert resp.status_code == 200, resp.text
    return {e: d for e, d in _parse_sse(resp.text)}


def test_the_preliminary_payload_carries_an_interpret_token(monkeypatch):
    """So the button is on screen with the tier, not minutes later."""
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module, "ai_auth_required", lambda: False)
    events = _events_for(_MINIMAL_EVIDENCE, ["pubtator3", "pmcoa"])
    prelim = events.get("preliminary_classification")
    assert prelim is not None, events.keys()
    assert prelim.get("interpret_token"), (
        "no interpret_token on the preliminary payload — the frontend renders "
        "the summary card without an 'Add AI interpretation' button")


def test_the_same_token_is_reused_by_the_final_event(monkeypatch):
    """Two tokens for one button would leave the browser holding the one with
    pre-literature criteria."""
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module, "ai_auth_required", lambda: False)
    events = _events_for(_MINIMAL_EVIDENCE, ["pubtator3", "pmcoa"])
    prelim = events["preliminary_classification"]
    final = events["db_only_complete"]
    assert final.get("interpret_token") == prelim.get("interpret_token"), (
        "the final event minted a NEW token; the browser is still holding the "
        "preliminary one")


def test_no_token_when_the_server_has_no_ai(monkeypatch):
    """The UI keeps its sign-in pill rather than offering a button that cannot
    work."""
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", False)
    events = _events_for(_MINIMAL_EVIDENCE, ["pubtator3"])
    prelim = events.get("preliminary_classification")
    assert prelim is not None
    assert not prelim.get("interpret_token")


def test_evidence_only_mode_leaves_the_button_enabled(monkeypatch):
    """`ai_unavailable` is the DISABLED reason. Opting out of AI is not a
    failure, so it must stay None or the button renders greyed out."""
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module, "ai_auth_required", lambda: False)
    events = _events_for(_MINIMAL_EVIDENCE, ["pmcoa"])
    assert events["preliminary_classification"].get("ai_unavailable") is None
