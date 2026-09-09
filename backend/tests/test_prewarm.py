"""Tests for C1 — activity-gated system-prompt cache pre-warming.

The per-curation cost is dominated by the COLD-cache WRITE of the ~18.9K-token
system prompt; warming keeps that block hot so real curations land as cheap
reads. This file covers the three contracts that matter:

  1. ``warm_system_cache`` sends ``system=_system_blocks()`` with ``max_tokens=1``,
     tallies its spend via ``budget.record_usage``, and SWALLOWS a client
     exception without raising (a warm failure must never disturb the server).
  2. The activity gate (``should_warm``): warms when activity is recent AND the
     budget is OK; SKIPS when idle (no recent activity) and SKIPS when
     ``ai_within_budget()`` is False — so an unused tool costs nothing.
  3. Default OFF: with ``HEARTVAR_CACHE_PREWARM`` unset, ``maybe_start_prewarm``
     schedules NO task and the lifespan startup path does not warm.

Fully offline — the Anthropic client is a fake, never the network; the gate is
driven by an injected ``now`` so no test sleeps on a real timer.

Runnable with pytest or directly
(``python -m backend.tests.test_prewarm``) — the __main__ runner executes the
async tests via asyncio.run and resets module state between them.
"""
from __future__ import annotations

import asyncio
from time import monotonic

import backend.claude as claude
from backend import budget, prewarm
from backend.claude import warm_system_cache


class _Usage:
    """Minimal stand-in for the SDK usage object (record_usage + the warm log
    only read these via getattr)."""

    def __init__(self, input_tokens=0, output_tokens=0,
                 cache_creation_input_tokens=0, cache_read_input_tokens=0):
        self.input_tokens = input_tokens
        self.output_tokens = output_tokens
        self.cache_creation_input_tokens = cache_creation_input_tokens
        self.cache_read_input_tokens = cache_read_input_tokens


class _Msg:
    def __init__(self, usage):
        self.usage = usage


class _FakeMessages:
    """Captures the kwargs of the single create() call the warmer makes."""

    def __init__(self, usage=None, raise_exc=None):
        self.calls: list[dict] = []
        self._usage = usage or _Usage(cache_creation_input_tokens=18900)
        self._raise = raise_exc

    async def create(self, **kwargs):
        self.calls.append(kwargs)
        if self._raise is not None:
            raise self._raise
        return _Msg(self._usage)


class _FakeClient:
    def __init__(self, messages):
        self.messages = messages


def _install_client(messages) -> None:
    claude._client = _FakeClient(messages)


def _reset_state() -> None:
    claude._client = None
    prewarm._last_activity = None
    budget._LEDGER = budget._DailyLedger()


def test_warm_sends_system_blocks_max_tokens_1_and_records_usage():
    """Warm sends the IDENTICAL production system block with max_tokens=1 and the
    pinned temperature, on the production model, and tallies its spend."""
    _reset_state()
    msgs = _FakeMessages(usage=_Usage(cache_creation_input_tokens=18900))
    _install_client(msgs)
    recorded: list[tuple[str, object]] = []
    orig_record = budget.record_usage
    budget.record_usage = lambda model, usage: recorded.append((model, usage))  # type: ignore[assignment]
    try:
        asyncio.run(warm_system_cache())
    finally:
        budget.record_usage = orig_record  # type: ignore[assignment]
        _reset_state()

    assert len(msgs.calls) == 1, "warm must make exactly one create() call"
    kw = msgs.calls[0]
    assert kw["max_tokens"] == 1, "warm must use the minimal 1-token cap"
    assert kw["model"] == claude.MODEL
    assert kw["temperature"] == claude.TEMPERATURE
    assert kw["system"] == claude._system_blocks(), "warm must send the prod system block"
    assert kw["messages"] == [{"role": "user", "content": "warm"}]
    assert len(recorded) == 1
    assert recorded[0][0] == claude.MODEL


def test_warm_swallows_client_exception_without_raising():
    """A failing client must be logged + swallowed — warming can never disturb
    the server (the next real curation just pays the cold write)."""
    _reset_state()
    msgs = _FakeMessages(raise_exc=RuntimeError("anthropic is down"))
    _install_client(msgs)
    try:
        asyncio.run(warm_system_cache())
    finally:
        _reset_state()
    assert len(msgs.calls) == 1, "the create() call was attempted before failing"


def test_warm_records_usage_into_real_ledger():
    """End-to-end: a warm that returns cache-creation tokens lands a (small) cost
    in the real budget ledger — proving the warm spend is genuinely tallied."""
    _reset_state()
    msgs = _FakeMessages(usage=_Usage(cache_creation_input_tokens=18900))
    _install_client(msgs)
    try:
        asyncio.run(warm_system_cache())
        assert budget._LEDGER.calls == 1
        assert budget._LEDGER.cost_usd > 0.0
    finally:
        _reset_state()


def test_gate_warms_when_recent_activity_and_budget_ok(monkeypatch):
    """Recent curation + budget OK -> the gate is open."""
    _reset_state()
    monkeypatch.setattr(budget, "ai_within_budget", lambda: True)
    base = monotonic()
    prewarm._last_activity = base
    assert prewarm.should_warm(now=base + 60.0) is True
    _reset_state()


def test_gate_skips_when_idle(monkeypatch):
    """No curation within the active window -> warming STOPS (idle costs nothing)."""
    _reset_state()
    monkeypatch.setattr(budget, "ai_within_budget", lambda: True)
    base = monotonic()
    prewarm._last_activity = base
    assert prewarm.should_warm(now=base + 1801.0) is False
    _reset_state()


def test_gate_skips_when_no_activity_at_all(monkeypatch):
    """No curation has happened since startup -> never warm (avoids warming an
    idle freshly-started server)."""
    _reset_state()
    monkeypatch.setattr(budget, "ai_within_budget", lambda: True)
    assert prewarm._last_activity is None
    assert prewarm.should_warm() is False
    _reset_state()


def test_gate_skips_when_over_budget(monkeypatch):
    """Recent activity but the daily budget is exhausted -> the gate is closed so
    warming can't push spend past the cap."""
    _reset_state()
    monkeypatch.setattr(budget, "ai_within_budget", lambda: False)
    base = monotonic()
    prewarm._last_activity = base
    assert prewarm.should_warm(now=base + 1.0) is False
    _reset_state()


def test_note_activity_arms_the_gate(monkeypatch):
    """note_activity() (called at curate start) sets the timestamp so the gate
    opens immediately afterwards."""
    _reset_state()
    monkeypatch.setattr(budget, "ai_within_budget", lambda: True)
    assert prewarm.should_warm() is False
    prewarm.note_activity()
    assert prewarm.should_warm() is True
    _reset_state()


def test_active_window_is_env_overridable(monkeypatch):
    """HEARTVAR_CACHE_PREWARM_ACTIVE_WINDOW widens/narrows the idle gate."""
    _reset_state()
    monkeypatch.setattr(budget, "ai_within_budget", lambda: True)
    monkeypatch.setenv("HEARTVAR_CACHE_PREWARM_ACTIVE_WINDOW", "10")
    base = monotonic()
    prewarm._last_activity = base
    assert prewarm.should_warm(now=base + 5.0) is True
    assert prewarm.should_warm(now=base + 11.0) is False
    monkeypatch.setenv("HEARTVAR_CACHE_PREWARM_ACTIVE_WINDOW", "not-a-number")
    assert prewarm._active_window() == 1800.0
    _reset_state()


def test_prewarm_disabled_by_default(monkeypatch):
    """With HEARTVAR_CACHE_PREWARM unset the feature is OFF."""
    monkeypatch.delenv("HEARTVAR_CACHE_PREWARM", raising=False)
    assert prewarm.prewarm_enabled() is False


def test_maybe_start_schedules_no_task_when_disabled(monkeypatch):
    """Default-off proof: maybe_start_prewarm returns None and creates NO task
    (so NO warm call ever happens) when the env is unset. We assert the negative
    by tracking asyncio.create_task — it must never be invoked, and any coroutine
    handed to it would be closed to avoid an un-awaited warning."""
    monkeypatch.delenv("HEARTVAR_CACHE_PREWARM", raising=False)
    created: list[object] = []

    def _spy(coro, *a, **k):  # pragma: no cover — must never run when disabled
        created.append(coro)
        coro.close()
        raise AssertionError("create_task must not be called when pre-warming is off")

    monkeypatch.setattr(asyncio, "create_task", _spy)

    async def _drive():
        return prewarm.maybe_start_prewarm()

    result = asyncio.run(_drive())
    assert result is None, "no task when disabled"
    assert created == [], "create_task must not be called when pre-warming is off"


def test_maybe_start_schedules_task_when_enabled(monkeypatch):
    """When HEARTVAR_CACHE_PREWARM is truthy the loop task IS created. We cancel
    it immediately so no real warming or sleeping happens in the test."""
    monkeypatch.setenv("HEARTVAR_CACHE_PREWARM", "1")
    assert prewarm.prewarm_enabled() is True

    async def _drive():
        task = prewarm.maybe_start_prewarm()
        assert task is not None
        task.cancel()
        try:
            await task
        except asyncio.CancelledError:
            pass
        return task

    task = asyncio.run(_drive())
    assert task.cancelled() or task.done()


def test_truthy_recognises_common_spellings(monkeypatch):
    for val in ("1", "true", "TRUE", "Yes", "on"):
        monkeypatch.setenv("HEARTVAR_CACHE_PREWARM", val)
        assert prewarm.prewarm_enabled() is True, val
    for val in ("0", "false", "no", "off", ""):
        monkeypatch.setenv("HEARTVAR_CACHE_PREWARM", val)
        assert prewarm.prewarm_enabled() is False, val


def test_loop_warms_then_cancels_cleanly(monkeypatch):
    """Drive ONE loop tick: a zero-interval loop with an open gate calls
    warm_system_cache exactly once before we cancel it. Proves the loop wires the
    gate to the warm fn and cancels cleanly (no orphaned task / exit warning)."""
    _reset_state()
    monkeypatch.setattr(prewarm, "_prewarm_interval", lambda: 0.0)
    monkeypatch.setattr(prewarm, "should_warm", lambda now=None: True)
    warmed = {"n": 0}

    async def _fake_warm():
        warmed["n"] += 1
        raise asyncio.CancelledError

    monkeypatch.setattr(prewarm, "warm_system_cache", _fake_warm)

    async def _drive():
        try:
            await prewarm._prewarm_loop()
        except asyncio.CancelledError:
            pass

    asyncio.run(_drive())
    assert warmed["n"] == 1
    _reset_state()


def test_loop_skips_warm_when_gate_closed(monkeypatch):
    """A closed gate (idle / over budget) means the loop ticks but does NOT warm."""
    _reset_state()
    ticks = {"n": 0}

    def _interval():
        return 0.0

    def _gate(now=None):
        ticks["n"] += 1
        if ticks["n"] >= 2:
            raise asyncio.CancelledError
        return False

    async def _should_not_warm():  # pragma: no cover
        raise AssertionError("warm must not run when the gate is closed")

    monkeypatch.setattr(prewarm, "_prewarm_interval", _interval)
    monkeypatch.setattr(prewarm, "should_warm", _gate)
    monkeypatch.setattr(prewarm, "warm_system_cache", _should_not_warm)

    async def _drive():
        try:
            await prewarm._prewarm_loop()
        except asyncio.CancelledError:
            pass

    asyncio.run(_drive())
    assert ticks["n"] >= 1
    _reset_state()


def _run_all():
    """Best-effort standalone runner for the monkeypatch-free async tests."""
    fns = [
        ("warm_sends", test_warm_sends_system_blocks_max_tokens_1_and_records_usage),
        ("warm_swallows", test_warm_swallows_client_exception_without_raising),
        ("warm_ledger", test_warm_records_usage_into_real_ledger),
    ]
    for name, fn in fns:
        fn()
        print(f"  ok  {name}")
    print(f"\n{len(fns)} non-monkeypatch tests passed "
          f"(run the full file with pytest for the gate/loop tests)")


if __name__ == "__main__":
    _run_all()
