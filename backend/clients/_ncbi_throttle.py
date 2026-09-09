"""Process-wide rate limiter + retry helper for NCBI E-utilities calls.

NCBI publishes a 3 req/s cap unauthenticated and 10 req/s with an
``api_key``. HeartVar has multiple clients that target E-utilities
concurrently inside a single ``/api/curate`` request — currently
``medgen`` (esearch + elink + esummary) and ``pubmed`` (variant
esearch + efetch, and the separate gene-literature esearch + efetch).
Without a shared limiter the three first-wave requests fire at t=0,
the three second-wave requests fire at t≈0.34s, and NCBI returns
HTTP 429 ``{"error":"API rate limit exceeded", …}`` on the fourth
request inside the 1-second sliding window. The medgen client would
then surface that 429 as a misleading ``"gene symbol not found in
NCBI Gene"`` error or as a silent empty conditions list.

Separately, NCBI's elink occasionally returns HTTP 200 with a body
containing an unescaped C++ exception trace ("TXCLIENT … Read failed:
EOF") — valid-looking JSON in the prefix but with an unescaped newline
that breaks parsing. The retry helper treats parse failures and bodies
carrying a top-level ``ERROR`` key the same as 429: short backoff,
retry up to ``MAX_ATTEMPTS`` times.

Both the throttle and the retry helper are process-wide; concurrent
``/api/curate`` requests share the same throttle.

TWO BUCKETS. Only ``eutils.ncbi.nlm.nih.gov`` (used by ``pubmed`` and
``medgen``) honours the E-utilities 3→10 req/s ``api_key`` uplift. The
other NCBI services we hit — PubTator3 (``/research/pubtator3-api/``),
PMC OA (``/pmc/utils/oa/``) and BioC (``/research/bionlp/``) — are
separate hosts with their own, stricter, undocumented per-IP limits the
key does NOT raise. Driving them at the keyed 10 req/s rate is exactly
what earns a temporary block. They therefore get a SEPARATE, always-
conservative bucket (:func:`web_throttle` / :func:`web_get_json`) that
never speeds up when a key is present.
"""

from __future__ import annotations

import asyncio
import logging
import os
import time

import httpx

from ._http_retry import _with_connect_cap, _retry_after_seconds

log = logging.getLogger("heartvar.ncbi")

NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "").strip() or None

NCBI_TOOL = "heartvar"
NCBI_CONTACT_EMAIL = (
    os.environ.get("NCBI_EMAIL", "").strip() or "heartvar@victorchang.edu.au"
)
NCBI_USER_AGENT = f"{NCBI_TOOL}/1.0 (variant-curation; {NCBI_CONTACT_EMAIL})"

_MIN_INTERVAL = 0.105 if NCBI_API_KEY else 0.355
log.info(
    "NCBI E-utilities throttle: ~%.1f req/s (%s)",
    1.0 / _MIN_INTERVAL,
    "NCBI_API_KEY present"
    if NCBI_API_KEY
    else "no NCBI_API_KEY — set it to raise the cap 3->10 req/s",
)

_WEB_MIN_INTERVAL = 0.34

MAX_ATTEMPTS = 4
_BACKOFF = (0.4, 1.2, 3.0)

_lock = asyncio.Lock()
_last_call: float = 0.0

_web_lock = asyncio.Lock()
_web_last_call: float = 0.0


class NCBIError(Exception):
    """Raised when an E-utilities request fails after all retries, or
    fails for a non-transient reason (e.g. HTTP 4xx other than 429)."""


async def throttle() -> None:
    """Acquire the global NCBI rate-limit slot.

    Blocks until enough time has elapsed since the previous NCBI
    request that the next one can be issued without exceeding the
    documented per-IP cap. Called internally by :func:`eutils_get_json`
    / :func:`eutils_get_text`; callers should normally use those
    helpers rather than touching ``throttle`` directly.
    """
    global _last_call
    async with _lock:
        wait = (_last_call + _MIN_INTERVAL) - time.monotonic()
        if wait > 0:
            await asyncio.sleep(wait)
        _last_call = time.monotonic()


async def web_throttle() -> None:
    """Acquire the conservative rate-limit slot for NCBI's non-E-utilities
    web services (PubTator3, PMC OA, BioC).

    Held at ~3 req/s and independent of both the E-utilities throttle and
    ``NCBI_API_KEY`` — those services live on separate hosts that don't
    honour the key, so they must never inherit the keyed 10 req/s rate.
    Called internally by :func:`web_get_json`; the PMC client calls it
    directly because it issues raw XML/JSON GETs rather than going through
    the retry helper.
    """
    global _web_last_call
    async with _web_lock:
        wait = (_web_last_call + _WEB_MIN_INTERVAL) - time.monotonic()
        if wait > 0:
            await asyncio.sleep(wait)
        _web_last_call = time.monotonic()


def _retry_wait(attempt: int) -> float:
    """Return the backoff in seconds after ``attempt`` (0-indexed)."""
    if attempt < len(_BACKOFF):
        return _BACKOFF[attempt]
    return _BACKOFF[-1]


async def _get_json(
    client: httpx.AsyncClient,
    url: str,
    params: dict,
    timeout: float,
    throttle_fn,
) -> dict:
    """Shared throttled-GET-with-retry core used by both buckets.

    ``throttle_fn`` selects the rate-limit slot: :func:`throttle` for
    E-utilities (key-aware) or :func:`web_throttle` for the conservative
    non-E-utilities services. Retries up to :data:`MAX_ATTEMPTS` times on
    transient failures: HTTP 429, transport errors, JSON parse failures
    (NCBI returning an unescaped C++ exception trace inside a JSON
    envelope), and 200-status responses that carry a top-level
    ``"ERROR"`` field (also indicative of a backend hiccup).

    Raises :class:`NCBIError` on permanent failure (e.g. 4xx other
    than 429) or when all retries have been exhausted.
    """
    last_err: str = "no attempts made"
    for attempt in range(MAX_ATTEMPTS):
        await throttle_fn()
        try:
            r = await client.get(url, params=params, timeout=_with_connect_cap(timeout))
        except httpx.HTTPError as e:
            last_err = f"transport error: {e!r}"
            log.warning("NCBI %s attempt %d: %s", url.rsplit("/", 1)[-1], attempt + 1, last_err)
            await asyncio.sleep(_retry_wait(attempt))
            continue
        if r.status_code == 429:
            last_err = "HTTP 429 rate-limited"
            log.warning("NCBI %s attempt %d: 429 rate-limited", url.rsplit("/", 1)[-1], attempt + 1)
            await asyncio.sleep(max(_retry_wait(attempt), _retry_after_seconds(r) or 0.0))
            continue
        if r.status_code != 200:
            raise NCBIError(f"HTTP {r.status_code}")
        try:
            payload = r.json()
        except ValueError:
            last_err = f"malformed JSON ({r.text[:120]!r})"
            log.warning("NCBI %s attempt %d: %s", url.rsplit("/", 1)[-1], attempt + 1, last_err)
            await asyncio.sleep(_retry_wait(attempt))
            continue
        if isinstance(payload, dict) and payload.get("ERROR"):
            last_err = f"NCBI ERROR field: {str(payload['ERROR'])[:120]}"
            log.warning("NCBI %s attempt %d: %s", url.rsplit("/", 1)[-1], attempt + 1, last_err)
            await asyncio.sleep(_retry_wait(attempt))
            continue
        return payload
    raise NCBIError(f"all {MAX_ATTEMPTS} retries failed: {last_err}")


async def eutils_get_json(
    client: httpx.AsyncClient,
    url: str,
    params: dict,
    timeout: float,
) -> dict:
    """Throttled JSON GET against E-utilities (``eutils.ncbi.nlm.nih.gov``).

    Uses the key-aware :func:`throttle` (3 req/s, or 10 with
    ``NCBI_API_KEY``). Callers must inject ``api_key`` into ``params``
    themselves when set — see ``pubmed`` / ``medgen``.
    """
    return await _get_json(client, url, params, timeout, throttle)


async def web_get_json(
    client: httpx.AsyncClient,
    url: str,
    params: dict,
    timeout: float,
) -> dict:
    """Throttled JSON GET against NCBI's non-E-utilities services
    (PubTator3, BioC). Uses the conservative, key-independent
    :func:`web_throttle` — the API key does not raise these services'
    limits, so they must never be driven at the keyed eutils rate.
    """
    return await _get_json(client, url, params, timeout, web_throttle)


async def eutils_get_text(
    client: httpx.AsyncClient,
    url: str,
    params: dict,
    timeout: float,
) -> str:
    """Issue a throttled GET to ``url`` and return the raw response text.

    Used for ``retmode=xml`` E-utilities calls (e.g. PubMed efetch)
    where the caller does its own XML parsing. Retries on 429 and
    transport errors. Raises :class:`NCBIError` on permanent failure
    or after all retries.
    """
    last_err: str = "no attempts made"
    for attempt in range(MAX_ATTEMPTS):
        await throttle()
        try:
            r = await client.get(url, params=params, timeout=_with_connect_cap(timeout))
        except httpx.HTTPError as e:
            last_err = f"transport error: {e!r}"
            log.warning("NCBI %s attempt %d: %s", url.rsplit("/", 1)[-1], attempt + 1, last_err)
            await asyncio.sleep(_retry_wait(attempt))
            continue
        if r.status_code == 429:
            last_err = "HTTP 429 rate-limited"
            log.warning("NCBI %s attempt %d: 429 rate-limited", url.rsplit("/", 1)[-1], attempt + 1)
            await asyncio.sleep(max(_retry_wait(attempt), _retry_after_seconds(r) or 0.0))
            continue
        if r.status_code != 200:
            raise NCBIError(f"HTTP {r.status_code}")
        return r.text
    raise NCBIError(f"all {MAX_ATTEMPTS} retries failed: {last_err}")
