"""Unit tests for backend.claude.stream_claude transient-failure retry.

A dropped connection / read-timeout mid-stream raises APITimeoutError, which
the Anthropic SDK does NOT retry once streaming has begun — silently dropping
the whole classification (observed: 2/3 runs in the 2026-06-04 eRepo
re-baseline each lost a variant this way). stream_claude restarts the stream
with backoff and yields a STREAM_RESTART sentinel so consumers discard the
failed attempt's partial chunks. Non-transient errors must propagate unchanged.

No pytest dependency — runnable with pytest or directly
(``python -m backend.tests.test_stream_claude_retry``).
"""
from __future__ import annotations

import asyncio
import os

import httpx
from anthropic import APITimeoutError

import backend.claude as claude
from backend.claude import STREAM_RESTART, stream_claude


def _timeout() -> APITimeoutError:
    return APITimeoutError(request=httpx.Request("POST", "https://api.anthropic.com"))


class _Usage:
    output_tokens = 10
    input_tokens = 20


class _FinalMsg:
    usage = _Usage()
    stop_reason = "end_turn"


class _FakeStream:
    """Async-context-manager stand-in for client.messages.stream(...)."""

    def __init__(self, chunks, fail_at=None, exc=None):
        self._chunks = chunks
        self._fail_at = fail_at
        self._exc = exc or _timeout()

    async def __aenter__(self):
        return self

    async def __aexit__(self, exc_type, exc, tb):
        return False

    @property
    def text_stream(self):
        return self._gen()

    async def _gen(self):
        for i, c in enumerate(self._chunks):
            if self._fail_at is not None and i == self._fail_at:
                raise self._exc
            yield c

    async def get_final_message(self):
        return _FinalMsg()


class _FakeClient:
    def __init__(self, streams):
        self._it = iter(streams)
        self.calls = 0

        class _Messages:
            def stream(_self, **kwargs):
                self.calls += 1
                return next(self._it)

        self.messages = _Messages()


def _with_client(streams):
    """Install a fake client for the duration of one collect()."""
    fake = _FakeClient(streams)
    claude._client = fake
    return fake


async def _collect():
    out = []
    async for c in stream_claude("prompt"):
        out.append(c)
    return out


def test_transient_then_success_restarts_and_recovers():
    """First attempt drops mid-stream; retry re-emits cleanly, with a
    STREAM_RESTART sentinel separating the discarded partial from the good run."""
    fake = _with_client([
        _FakeStream(["partial", "BOOM"], fail_at=1),
        _FakeStream(['{"tier"', ': "P"}'], fail_at=None),
    ])
    try:
        out = asyncio.run(_collect())
    finally:
        claude._client = None

    assert STREAM_RESTART in out, out
    idx = out.index(STREAM_RESTART)
    assert out[:idx] == ["partial"], out[:idx]
    assert "".join(out[idx + 1:]) == '{"tier": "P"}', out[idx + 1:]
    assert fake.calls == 2


def test_retries_exhausted_raises():
    """Every attempt fails transiently → APITimeoutError propagates after the
    max number of attempts, with a restart sentinel before each retry."""
    streams = [_FakeStream(["x", "BOOM"], fail_at=1) for _ in range(claude._CLAUDE_MAX_ATTEMPTS)]
    fake = _with_client(streams)
    try:
        raised = False
        try:
            asyncio.run(_collect())
        except APITimeoutError:
            raised = True
    finally:
        claude._client = None
    assert raised, "expected APITimeoutError after exhausting retries"
    assert fake.calls == claude._CLAUDE_MAX_ATTEMPTS


def test_non_retryable_propagates_without_retry():
    """A non-transient error is not retried — it surfaces on the first attempt."""
    fake = _with_client([_FakeStream(["x", "BOOM"], fail_at=1, exc=ValueError("nope"))])
    try:
        raised = False
        try:
            asyncio.run(_collect())
        except ValueError:
            raised = True
    finally:
        claude._client = None
    assert raised, "expected ValueError to propagate"
    assert fake.calls == 1, "non-retryable error must not trigger a retry"


def test_default_max_attempts_is_three():
    """Default cap is 3 attempts (= 2 retries) — covers single + double
    consecutive transient drops while bounding worst-case retry latency."""
    assert claude._CLAUDE_MAX_ATTEMPTS == 3


def test_env_override_max_attempts_respected():
    """HEARTVAR_CLAUDE_MAX_ATTEMPTS overrides the default via _int_env (the
    same parser the module-level constant uses), with a clean fallback on a
    non-integer value. Tested through the helper so we don't reload the module
    (reloading would rebind STREAM_RESTART and break the other tests)."""
    prev = os.environ.get("HEARTVAR_CLAUDE_MAX_ATTEMPTS")
    try:
        os.environ["HEARTVAR_CLAUDE_MAX_ATTEMPTS"] = "5"
        assert claude._int_env("HEARTVAR_CLAUDE_MAX_ATTEMPTS", 3) == 5
        os.environ["HEARTVAR_CLAUDE_MAX_ATTEMPTS"] = "not-a-number"
        assert claude._int_env("HEARTVAR_CLAUDE_MAX_ATTEMPTS", 3) == 3
        os.environ.pop("HEARTVAR_CLAUDE_MAX_ATTEMPTS", None)
        assert claude._int_env("HEARTVAR_CLAUDE_MAX_ATTEMPTS", 3) == 3
    finally:
        if prev is None:
            os.environ.pop("HEARTVAR_CLAUDE_MAX_ATTEMPTS", None)
        else:
            os.environ["HEARTVAR_CLAUDE_MAX_ATTEMPTS"] = prev


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
