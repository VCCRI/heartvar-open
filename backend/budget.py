"""Server-side AI cost accounting + a daily spend guard.

HeartVar's AI now runs entirely on the OWNER's Anthropic key (no more
bring-your-own-key), so every AI call bills to one account. The owner caps
hard spend at the Anthropic console (e.g. a $100/month limit), but that
console cap is a brick wall: once hit, every request 400s mid-stream for the
rest of the month and the site's AI goes dark with no graceful degradation.

This module is the softer, in-app guard that sits in front of that wall:

  * Every Anthropic call reports its token usage here via ``record_usage``.
  * Costs are tallied per UTC day against ``HEARTVAR_DAILY_USD_BUDGET``.
  * ``ai_within_budget()`` lets request handlers fall back to the deterministic
    evidence-only result (with a notice) BEFORE the console cap is reached,
    instead of failing hard.

A daily budget of $X bounds the worst-case month to ~30·$X, so a daily figure
comfortably below console_cap/30 keeps normal months under the hard cap while
still allowing day-to-day bursts.

Default **$100.00/day** (raised from $20.00 on 2026-09-09 as usage grew beyond
one institute; $20.00 itself was raised from $3.00 on 2026-08-13). At current
pricing $100 is roughly 650 AI curations a day. The request-count limits were
raised to 1000 in the same change, so this is now the BINDING limit on AI use —
by design, because a dollar ceiling degrades gracefully to evidence-only where a
request ceiling returns 429.

⚠️ The worst-case month is now ~$600, which is ABOVE a $100 Anthropic console cap.
The console cap is the real backstop for this guard (see the caveat below), so it
must be raised to match, or the month will end with hard API failures rather than
the graceful evidence-only degradation this module exists to provide.

CAVEAT — the ledger is in-process (single worker, like the rate limiter and the
shared cache) and resets on restart. A restart mid-day clears the day's tally,
so a determined month-long abuser is ultimately bounded only by the console
cap, not by this guard. For a durable cross-restart budget, back this with a
shared store (Redis) the same way the limiter can be. Good enough for a niche
research tool; the console cap is the real backstop.
"""

from __future__ import annotations

import json
import logging
import os
from datetime import datetime, timezone

log = logging.getLogger("heartvar.budget")

_DEFAULT_PRICING: dict[str, dict[str, float]] = {
    "claude-opus-4": {"in": 5.00, "out": 25.00},
    "claude-sonnet-4": {"in": 3.00, "out": 15.00},
    "claude-haiku-4": {"in": 1.00, "out": 5.00},
}

_CACHE_WRITE_MULT = 1.25
_CACHE_READ_MULT = 0.10

_FALLBACK_PRICING = {"in": 5.00, "out": 25.00}


def _load_pricing() -> dict[str, dict[str, float]]:
    raw = os.environ.get("HEARTVAR_MODEL_PRICING", "").strip()
    if not raw:
        return dict(_DEFAULT_PRICING)
    try:
        parsed = json.loads(raw)
        out: dict[str, dict[str, float]] = {}
        for key, val in parsed.items():
            in_v = val["in"] if "in" in val else val["input"]
            out_v = val["out"] if "out" in val else val["output"]
            out[str(key)] = {"in": float(in_v), "out": float(out_v)}
        log.info("Loaded model pricing override (%d entries)", len(out))
        return out
    except (ValueError, KeyError, TypeError) as e:
        log.warning(
            "Ignoring malformed HEARTVAR_MODEL_PRICING (%s) — using defaults", e
        )
        return dict(_DEFAULT_PRICING)


_PRICING = _load_pricing()


def _prices_for(model: str) -> dict[str, float]:
    model = (model or "").strip()
    for prefix, prices in _PRICING.items():
        if prefix in model:
            return prices
    return _FALLBACK_PRICING


def estimate_cost_usd(
    model: str,
    input_tokens: int = 0,
    output_tokens: int = 0,
    cache_creation_input_tokens: int = 0,
    cache_read_input_tokens: int = 0,
) -> float:
    """USD cost of one call from its token usage. Cache-write tokens bill at
    1.25× input, cache-read tokens at 0.10× input; plain input/output bill at
    the model's full rate."""
    p = _prices_for(model)
    in_rate = p["in"] / 1_000_000.0
    out_rate = p["out"] / 1_000_000.0
    return (
        input_tokens * in_rate
        + cache_creation_input_tokens * in_rate * _CACHE_WRITE_MULT
        + cache_read_input_tokens * in_rate * _CACHE_READ_MULT
        + output_tokens * out_rate
    )


def _daily_budget_usd() -> float:
    """Daily spend ceiling in USD. 0 (or unset) disables the guard — AI is then
    bounded only by the Anthropic console cap."""
    try:
        return float(os.environ.get("HEARTVAR_DAILY_USD_BUDGET", "") or 100.0)
    except ValueError:
        return 100.0


class _DailyLedger:
    """In-process per-UTC-day running tally of AI spend. Resets automatically on
    the first call of a new UTC day.

    ``reserved`` holds the estimated cost of in-flight calls that have passed the
    budget gate but not yet recorded actual usage. It closes the check-then-spend
    race: on the single async worker, without it N concurrent requests would all
    read ``cost_usd`` before any of them recorded a cost and all proceed. With
    reservations, each admitted call first books an estimate against
    ``cost_usd + reserved``; the estimate is released when actual usage lands."""

    def __init__(self) -> None:
        self.day: str | None = None
        self.cost_usd: float = 0.0
        self.reserved: float = 0.0
        self.calls: int = 0
        self.input_tokens: int = 0
        self.output_tokens: int = 0

    def _roll(self, today: str) -> None:
        if self.day != today:
            if self.day is not None:
                log.info(
                    "AI budget day rollover %s -> %s (prior day: $%.4f over %d calls)",
                    self.day, today, self.cost_usd, self.calls,
                )
            self.day = today
            self.cost_usd = 0.0
            self.calls = 0
            self.input_tokens = 0
            self.output_tokens = 0

    def add(self, cost: float, in_tok: int, out_tok: int, today: str) -> None:
        self._roll(today)
        self.cost_usd += cost
        self.calls += 1
        self.input_tokens += in_tok
        self.output_tokens += out_tok

    def spent_today(self, today: str) -> float:
        self._roll(today)
        return self.cost_usd


_LEDGER = _DailyLedger()


def _today() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d")


def record_usage(model: str, usage: object) -> float:
    """Tally one Anthropic call against today's budget from its ``usage`` object
    (the SDK ``Message.usage`` / final-message usage). Tolerant of missing
    fields. Returns the call's estimated USD cost (also logged). Never raises —
    accounting must not break a response."""
    try:
        in_tok = int(getattr(usage, "input_tokens", 0) or 0)
        out_tok = int(getattr(usage, "output_tokens", 0) or 0)
        cache_w = int(getattr(usage, "cache_creation_input_tokens", 0) or 0)
        cache_r = int(getattr(usage, "cache_read_input_tokens", 0) or 0)
        cost = estimate_cost_usd(model, in_tok, out_tok, cache_w, cache_r)
        today = _today()
        _LEDGER.add(cost, in_tok + cache_w + cache_r, out_tok, today)
        budget = _daily_budget_usd()
        log.info(
            "[budget] %s call $%.4f | today $%.4f%s | in=%d cache_w=%d cache_r=%d out=%d",
            model, cost, _LEDGER.cost_usd,
            f"/${budget:.2f}" if budget > 0 else " (no cap)",
            in_tok, cache_w, cache_r, out_tok,
        )
        return cost
    except Exception as e:  # noqa: BLE001 — accounting must never break a call
        log.warning("budget.record_usage failed (ignored): %r", e)
        return 0.0


def ai_within_budget() -> bool:
    """True if today's AI spend (recorded + in-flight reservations) is under the
    daily cap, or the cap is disabled. A read-only check — prefer ``reserve``
    for a request that is about to spend, so concurrent in-flight calls can't all
    pass the gate before any of them records a cost."""
    budget = _daily_budget_usd()
    if budget <= 0:
        return True
    today = _today()
    return (_LEDGER.spent_today(today) + _LEDGER.reserved) < budget


def reserve(est_usd: float) -> bool:
    """Atomically (on the single-threaded event loop) admit a call that is about
    to spend ~``est_usd``: book the estimate against recorded + already-reserved
    spend and return True, or return False if it would exceed the daily cap.
    The caller MUST pair every True with a ``release(est_usd)`` (in a finally)
    once the call finishes. When the cap is disabled (budget <= 0) this always
    admits and books nothing. Returns True without booking when est is
    non-positive so a bad estimate can't wedge the gate."""
    budget = _daily_budget_usd()
    if budget <= 0 or est_usd <= 0:
        return True
    today = _today()
    if (_LEDGER.spent_today(today) + _LEDGER.reserved + est_usd) >= budget:
        return False
    _LEDGER.reserved += est_usd
    return True


def release(est_usd: float) -> None:
    """Release a prior ``reserve`` estimate (actual cost is tallied separately by
    ``record_usage``). Never lets the reservation pool go negative."""
    if est_usd and est_usd > 0:
        _LEDGER.reserved = max(0.0, _LEDGER.reserved - est_usd)


def estimate_call_cost(model: str, est_input_tokens: int, max_output_tokens: int) -> float:
    """Conservative upfront cost estimate for the reservation, priced at the
    model's full input+output rate (no cache discount — deliberately
    over-estimating so the reservation can't under-book a burst)."""
    return estimate_cost_usd(model, est_input_tokens, max_output_tokens)
