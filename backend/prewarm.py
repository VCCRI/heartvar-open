"""Activity-gated system-prompt cache pre-warming (C1).

The per-curation cost is dominated by the COLD-cache WRITE of the ~18.9K-token
system prompt (≈$0.071/call at the 1.25× write multiplier). On a WARM cache the
SAME block is a cache READ (0.10×, ≈$0.006). Anthropic's ephemeral prompt cache
expires after its TTL (5-min default, or 1h via ``HEARTVAR_CACHE_TTL``); to keep
production curations landing as cheap reads we periodically re-send a minimal
request that re-writes the IDENTICAL system block before the TTL lapses.

DEFAULT OFF — with ``HEARTVAR_CACHE_PREWARM`` unset/falsey, NOTHING in this module
runs: ``maybe_start_prewarm`` does not create a task, no background loop ticks, no
extra Anthropic calls happen, and production behaviour is byte-identical to today.

Idle cost is bounded by an ACTIVITY GATE: the loop only warms when a real curation
happened within ``HEARTVAR_CACHE_PREWARM_ACTIVE_WINDOW`` seconds (default 1800 =
30 min). When the tool is idle longer than that window, warming STOPS and costs
nothing until the next real curation re-arms it. The loop also checks
``budget.ai_within_budget()`` before each warm, so warming can never push the
daily USD spend past the cap (and degrades to evidence-only-style silence exactly
like a real call would when over budget).

In-process state (single worker): the activity timestamp and the cache itself both
live in this process, so cache-warming REQUIRES a single uvicorn worker
(``--workers 1``). With multiple workers each has its own cache and its own
activity clock, so warming one worker does not warm the others.
"""

from __future__ import annotations

import asyncio
import logging
import os
from time import monotonic

from . import budget
from .claude import warm_system_cache

log = logging.getLogger("heartvar.prewarm")


def _truthy(name: str) -> bool:
    """A permissive truthy env check (``1``/``true``/``yes``/``on``, any case).
    Unset or blank or a falsey word -> False, so the feature is OFF by default."""
    return os.environ.get(name, "").strip().lower() in {"1", "true", "yes", "on"}


def _float_env(name: str, default: float) -> float:
    """Parse a float env var, falling back to ``default`` on unset/blank/invalid
    (mirrors backend.claude._int_env / the MAX_TOKENS pattern)."""
    try:
        return float(os.environ.get(name, "") or default)
    except ValueError:
        return default


def prewarm_enabled() -> bool:
    """Master switch. Default OFF — read live (not cached at import) so tests and
    deploys can flip it via the environment without reloading the module."""
    return _truthy("HEARTVAR_CACHE_PREWARM")


def _prewarm_interval() -> float:
    """Seconds between warm attempts. Default 240 (= 4 min), deliberately UNDER
    the 5-min ephemeral TTL so a warm always lands before the prior entry lapses.
    Pair with HEARTVAR_CACHE_TTL=1h to widen the warm-read window."""
    return _float_env("HEARTVAR_CACHE_PREWARM_INTERVAL", 240.0)


def _active_window() -> float:
    """Seconds of inactivity after which warming STOPS (bounds idle cost).
    Default 1800 (= 30 min)."""
    return _float_env("HEARTVAR_CACHE_PREWARM_ACTIVE_WINDOW", 1800.0)


_last_activity: float | None = None


def note_activity() -> None:
    """Record that a real curation just started. Called at the top of
    ``curate_stream``. Cheap (one timestamp write) and unconditional — it is safe
    to call even when pre-warming is disabled (it just updates a value nothing
    reads), so the curate path needs no feature-flag branch."""
    global _last_activity
    _last_activity = monotonic()


def _recent_activity(now: float | None = None) -> bool:
    """True if a real curation happened within the active window — the idle gate.
    ``now`` is injectable so tests can drive the gate deterministically without
    sleeping on real timers."""
    if _last_activity is None:
        return False
    now = monotonic() if now is None else now
    return (now - _last_activity) <= _active_window()


def should_warm(now: float | None = None) -> bool:
    """The pure gate predicate the loop consults each tick: warm only when a real
    curation is recent AND the daily AI budget still has room. Factored out so it
    is unit-testable without the loop's sleep. ``now`` is injectable for tests."""
    return _recent_activity(now) and budget.ai_within_budget()


async def _prewarm_loop() -> None:
    """Background loop: every interval, warm the cache iff the gate is open.

    The gate (``should_warm``) bounds idle cost — when the tool has been idle
    longer than the active window the loop wakes, sees no recent activity, and
    skips the warm (no Anthropic call), so an unused tool costs nothing. Cancelled
    cleanly on shutdown via ``asyncio.CancelledError`` (re-raised so the awaiter
    sees a clean cancellation, no orphaned task, no exit warning)."""
    interval = _prewarm_interval()
    log.info(
        "[prewarm] cache-warming ENABLED (interval=%.0fs, active_window=%.0fs) — "
        "requires a single worker (cache is in-process)",
        interval, _active_window(),
    )
    try:
        while True:
            await asyncio.sleep(interval)
            if should_warm():
                await warm_system_cache()
            else:
                log.debug("[prewarm] tick skipped (idle or over budget)")
    except asyncio.CancelledError:
        log.info("[prewarm] cache-warming loop cancelled (shutdown)")
        raise


def maybe_start_prewarm() -> asyncio.Task | None:
    """Start the background warming task IFF ``HEARTVAR_CACHE_PREWARM`` is truthy.

    Returns the created task (so the lifespan can cancel+await it on shutdown), or
    ``None`` when the feature is disabled — in which case NO task is scheduled and
    nothing changes. Call once from the FastAPI lifespan startup."""
    if not prewarm_enabled():
        return None
    return asyncio.create_task(_prewarm_loop())
