"""Process-wide async TTL + single-flight cache for expensive external lookups.

Two levers keep HeartVar inside Ensembl's / NCBI's per-IP rate limits under
multi-user load:

  * **TTL caching** — a resolved variant (VEP coordinates, PubMed/PubTator3
    hits, PMC full text) is stable for days-to-months, so a repeat curation of
    the same variant is served from memory instead of re-hitting the upstream.
  * **Single-flight** — if N concurrent ``/api/curate`` requests ask for the
    SAME key while the first call is still in flight, they all await that one
    upstream call rather than firing N identical requests at the service. This
    is what protects the shared NCBI token bucket when several curators land on
    the same variant at once.

SINGLE-WORKER ONLY. The cache lives in this process's memory, exactly like the
NCBI throttle (see ``_ncbi_throttle``). It holds only while the app runs as ONE
process: ``uvicorn ... --workers N``, gunicorn, or Azure scale-out across
instances each get their OWN cache (and own throttle), silently voiding both
guarantees. Cross-process sharing would need a shared backend (e.g. Redis).

── Horizontal-scaling migration (audit §7c; DESIGN ONLY — keep --workers=1 until
   a Redis instance is provisioned AND tested) ───────────────────────────────────
Three pieces of process-local state force single-worker today; to scale out,
move ALL three to Redis (no single one is sufficient — they protect different
limits):

  1. Inbound rate limiter (backend/app.py) — already wired as a safe win: set
     ``HEARTVAR_REDIS_URL`` and slowapi keys its windows in Redis (the Limiter
     is built with ``storage_uri`` when that env var is present). No interface
     change; defaults to in-process MemoryStorage.
  2. THIS TTL + single-flight cache — add a ``RedisTTLCache`` implementing the
     same ``get_or_set(key, factory, ttl)`` contract: ``GET`` key; on miss
     ``SET … PX=ttl*1000`` after the factory resolves. True cross-process
     single-flight needs a short distributed lock (``SET key:lock NX PX=…``) so
     only one worker calls the upstream; without it, accept best-effort
     per-worker coalescing (≤ N concurrent upstream calls for N workers — still
     far below unthrottled). Keep the existing namespaced-tuple keys, serialised
     to a stable string.
  3. Per-host outbound throttle (``_throttle.py``) + the NCBI 3/s bucket
     (``_ncbi_throttle.py``) — THE actual ban-prevention guarantee, and the one
     that MUST be shared before --workers>1: N workers with per-process throttles
     = N× the real outbound rate → an Ensembl/NCBI IP block. Replace the per-host
     min-interval + global semaphore with a Redis token-bucket / sliding-window
     (Lua ``INCR``+``PEXPIRE``) keyed per host, behind the same
     ``await throttle_host(url)`` signature so call sites don't change.

Select the backend via env (e.g. ``HEARTVAR_CACHE_BACKEND=memory|redis`` +
``HEARTVAR_REDIS_URL``). requirements: add ``redis`` (and ``limits[redis]`` for
slowapi). Until all three are Redis-backed, the Dockerfile MUST keep one worker.

No ``asyncio.Lock`` is used on purpose. The event loop is single-threaded, so
the synchronous bookkeeping between ``await`` points (store lookup, in-flight
registration) is atomic — a module-level ``asyncio.Lock`` would additionally
break tests that call ``asyncio.run`` more than once (the lock binds to the
first loop it touches).
"""
from __future__ import annotations

import asyncio
import logging
import time
from typing import Any, Awaitable, Callable, Hashable

log = logging.getLogger("heartvar.cache")

_DAY = 24 * 3600
TTL_VEP = 30 * _DAY
TTL_LITERATURE = 7 * _DAY
TTL_PMC = 90 * _DAY


class TTLCache:
    """In-process async cache with per-entry TTL and single-flight coalescing.

    ``get_or_set`` is the only entry point callers need. ``clear`` resets all
    state (tests use it via the autouse conftest fixture); ``stats`` exposes
    hit/miss/coalesced counters for observability.
    """

    def __init__(self, name: str, *, default_ttl: int = TTL_LITERATURE,
                 maxsize: int = 4096) -> None:
        self.name = name
        self._default_ttl = default_ttl
        self._maxsize = maxsize
        self._store: dict[Hashable, tuple[Any, float]] = {}
        self._inflight: dict[Hashable, asyncio.Future] = {}
        self._hits = 0
        self._misses = 0
        self._coalesced = 0

    async def get_or_set(
        self,
        key: Hashable,
        factory: Callable[[], Awaitable[Any]],
        *,
        ttl: int | None = None,
        should_cache: Callable[[Any], bool] | None = None,
    ) -> Any:
        """Return the cached value for ``key`` or compute it via ``factory``.

        ``factory`` is an argument-less coroutine factory run only on a miss.
        ``should_cache`` decides whether a freshly computed value is stored
        (default: always); a value that fails the predicate is still RETURNED
        to this caller and to any coalesced waiters, but not retained — so a
        transient failure isn't pinned for the whole TTL. Exceptions propagate
        to every waiter and are never cached.
        """
        now = time.monotonic()

        entry = self._store.get(key)
        if entry is not None and entry[1] > now:
            self._hits += 1
            return entry[0]
        if entry is not None:
            self._store.pop(key, None)

        inflight = self._inflight.get(key)
        if inflight is not None:
            self._coalesced += 1
            return await inflight

        self._misses += 1
        fut: asyncio.Future = asyncio.get_running_loop().create_future()
        self._inflight[key] = fut
        try:
            value = await factory()
        except BaseException as exc:
            self._inflight.pop(key, None)
            if not fut.done():
                fut.set_exception(exc)
                fut.exception()
            raise

        self._inflight.pop(key, None)
        if should_cache is None or should_cache(value):
            eff_ttl = self._default_ttl if ttl is None else ttl
            self._store[key] = (value, now + eff_ttl)
            self._evict()
        if not fut.done():
            fut.set_result(value)
        return value

    def _evict(self) -> None:
        """Bound the store by dropping the soonest-to-expire entries."""
        overflow = len(self._store) - self._maxsize
        if overflow <= 0:
            return
        for k in sorted(self._store, key=lambda k: self._store[k][1])[:overflow]:
            self._store.pop(k, None)

    def clear(self) -> None:
        self._store.clear()
        self._inflight.clear()
        self._hits = self._misses = self._coalesced = 0

    def stats(self) -> dict[str, int]:
        return {
            "entries": len(self._store),
            "inflight": len(self._inflight),
            "hits": self._hits,
            "misses": self._misses,
            "coalesced": self._coalesced,
        }


EXTERNAL_CACHE = TTLCache("external", default_ttl=TTL_LITERATURE)


def _ok(result: Any) -> bool:
    """should_cache predicate for the ``{"ok": bool, ...}`` client shape —
    only successful lookups are retained."""
    return isinstance(result, dict) and result.get("ok") is True
