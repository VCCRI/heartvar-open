"""Shared retry + backoff wrapper for HeartVar's external-DB HTTP clients.

Several public-database clients (gnomAD, PanelApp, UniProt, GTEx, ProtVar,
Open Targets, MGI) historically made a SINGLE HTTP attempt with no retry, so
one transient blip — a read timeout, a dropped keep-alive connection, a one-off
429/5xx from the public service — permanently dropped that source for the whole
curation. Because ~9 such single-attempt calls fire in parallel on every
request, a different random one or two failed each run (gnomAD one run, PanelApp
the next). Routing them through ``request_with_retry`` gives them the same
resilience the NCBI / Ensembl / SpliceAI clients already had, which eliminates
the rotating "1-2 databases failed" symptom.
"""
import asyncio
import email.utils
import logging
import os
import random
import time

import httpx

from ._throttle import global_semaphore, throttle_host

log = logging.getLogger("heartvar.http_retry")

CONTACT_USER_AGENT = os.environ.get(
    "HEARTVAR_USER_AGENT",
    "HeartVar/1.0 (cardiac variant curation; +mailto:heartvar@victorchang.edu.au)",
)


def make_async_client(**kwargs) -> httpx.AsyncClient:
    """``httpx.AsyncClient`` with a default contact User-Agent header. Any
    User-Agent the caller passes in ``headers`` wins, so the NCBI / Ensembl /
    eRepo / PanelApp clients keep their own."""
    headers = {"User-Agent": CONTACT_USER_AGENT}
    headers.update(kwargs.pop("headers", None) or {})
    return httpx.AsyncClient(headers=headers, **kwargs)

_BACKOFF = (0.5, 1.5, 3.0)
_RETRY_AFTER_CAP = 10.0

_DEFAULT_RETRY_DEADLINE = 8.0

_DEFAULT_CONNECT_TIMEOUT = 3.0


def _connect_timeout() -> float:
    """Per-attempt CONNECT budget. Falls back on unset/blank/junk/non-positive."""
    raw = os.environ.get("HEARTVAR_HTTP_CONNECT_TIMEOUT", "").strip()
    if raw:
        try:
            value = float(raw)
            if value > 0:
                return value
        except ValueError:
            pass
        log.warning("HEARTVAR_HTTP_CONNECT_TIMEOUT=%r unusable — using %.1fs",
                    raw, _DEFAULT_CONNECT_TIMEOUT)
    return _DEFAULT_CONNECT_TIMEOUT


def _with_connect_cap(timeout):
    """Normalise a caller's ``timeout`` so CONNECT is bounded separately.

    A scalar becomes an ``httpx.Timeout`` that KEEPS the caller's budget for
    read/write/pool and caps connect. An explicit ``httpx.Timeout`` is returned
    untouched — the caller has already decided, and overriding that would be
    surprising. Done here rather than at ~20 call sites so no client can
    reintroduce the problem by passing a bare number.
    """
    if isinstance(timeout, httpx.Timeout):
        return timeout
    if not isinstance(timeout, (int, float)):
        return timeout
    connect = min(_connect_timeout(), float(timeout))
    return httpx.Timeout(float(timeout), connect=connect)


def _retry_deadline() -> float:
    """Total wall-clock budget for all attempts at one request.

    Falls back on unset/blank/junk/non-positive, following the codebase's
    convention that a typo in an env var must never change behaviour silently.
    """
    raw = os.environ.get("HEARTVAR_HTTP_RETRY_DEADLINE", "").strip()
    if raw:
        try:
            value = float(raw)
            if value > 0:
                return value
        except ValueError:
            pass
        log.warning("HEARTVAR_HTTP_RETRY_DEADLINE=%r unusable — using %.1fs",
                    raw, _DEFAULT_RETRY_DEADLINE)
    return _DEFAULT_RETRY_DEADLINE


def _retry_after_seconds(resp: httpx.Response) -> float | None:
    """Parse a ``Retry-After`` header (delta-seconds or HTTP-date) into seconds,
    capped at :data:`_RETRY_AFTER_CAP`. Returns ``None`` when the header is
    absent, non-positive, or unparseable — so hosts that don't send it (most of
    them) fall straight through to the static backoff schedule."""
    val = (resp.headers.get("retry-after") or "").strip()
    if not val:
        return None
    try:
        secs = float(val)
    except ValueError:
        try:
            dt = email.utils.parsedate_to_datetime(val)
        except (TypeError, ValueError):
            return None
        if dt is None:
            return None
        secs = dt.timestamp() - time.time()
    if secs <= 0:
        return None
    return min(secs, _RETRY_AFTER_CAP)


async def request_with_retry(client, method, url, *, timeout,
                             max_attempts=4, name="", logger=None,
                             deadline=None, **kwargs):
    """Issue an HTTP request with bounded retry + jittered backoff, governed by
    the process-wide per-host outbound throttle (see ._throttle).

    Retries on transport errors (``httpx.HTTPError``) and on HTTP 429 / >=500.
    A 200 — or any non-429 4xx — is treated as terminal and returned
    immediately (callers inspect ``status_code`` exactly as before). Returns
    the final ``httpx.Response``, or ``None`` if every attempt failed
    transiently (caller should then return its own ``ok: False`` record).

    Each attempt first waits for its turn under the host's rate limit, then
    issues the request inside the global concurrency cap, so the outbound rate
    to any public service stays bounded regardless of inbound traffic.

    ``deadline`` is the TOTAL wall-clock budget across all attempts (default
    :func:`_retry_deadline`). It exists because ``max_attempts`` bounded the
    count but nothing bounded the time — see the note beside
    ``_DEFAULT_RETRY_DEADLINE`` for the production measurement. Checked between
    attempts, so the worst case is the deadline plus one in-flight attempt.
    """
    lg = logger or log
    last_err = "no attempts"
    budget = _retry_deadline() if deadline is None else deadline
    timeout = _with_connect_cap(timeout)
    started = time.monotonic()
    for attempt in range(max_attempts):
        if attempt and (time.monotonic() - started) >= budget:
            lg.warning("%s gave up after %d attempt(s) — %.1fs retry deadline "
                       "reached (last: %s)", name or url, attempt,
                       budget, last_err)
            return None
        retry_after: float | None = None
        try:
            await throttle_host(url)
            async with global_semaphore():
                r = await client.request(method, url, timeout=timeout, **kwargs)
            if r.status_code != 429 and r.status_code < 500:
                return r
            last_err = f"HTTP {r.status_code}"
            retry_after = _retry_after_seconds(r)
        except httpx.HTTPError as e:
            last_err = f"{type(e).__name__}: {e}"
        lg.warning("%s attempt %d/%d failed: %s",
                   name or url, attempt + 1, max_attempts, last_err)
        if attempt < max_attempts - 1:
            base = _BACKOFF[min(attempt, len(_BACKOFF) - 1)] + random.uniform(0, 0.3)
            wait = max(base, retry_after) if retry_after else base
            remaining = budget - (time.monotonic() - started)
            if remaining <= 0:
                continue
            await asyncio.sleep(min(wait, remaining))
    lg.warning("%s exhausted %d attempts: %s", name or url, max_attempts, last_err)
    return None
