"""Process-wide OUTBOUND rate governance for HeartVar's live public-database
clients.

Every /api/curate request fans out to ~13 external services in parallel, and
concurrent curations multiply that. Without a ceiling, a burst of users could
exceed a public service's rate limit and get HeartVar's *server* IP throttled
or blocked — breaking the tool for everyone. This module bounds the outbound
rate PER HOST (independent of how many inbound users there are) and caps total
simultaneous outbound connections. That is the real ban-prevention guarantee;
the inbound per-IP limiter in app.py only limits how often a user can ask.

It complements (does not replace) the service-specific throttles already in
place, which keep their own correct limits and do NOT route through here:
  - _ncbi_throttle.py — PubMed / MedGen / PubTator3 / PMC obey NCBI's 3 req/s;
  - erepo_client.py   — ClinGen eRepo at 1 req/s.

Wired in via request_with_retry (covers gnomAD, Ensembl, UniProt, GTEx,
PanelApp, ProtVar, Open Targets, MGI) and called directly by the
SpliceAI client.

PROCESS-LOCAL — this throttle MUST become Redis-backed before running more than
one worker (N workers = N× the real outbound rate → an upstream IP block). See
the "Horizontal-scaling migration" design block in ``_cache.py`` (audit §7c).
"""
from __future__ import annotations

import asyncio
import os
import time
from urllib.parse import urlparse

_HOST_MIN_INTERVAL = {
    "rest.ensembl.org": 0.10,
    "grch37.rest.ensembl.org": 0.10,
    "gnomad.broadinstitute.org": 0.20,
    "spliceailookup.broadinstitute.org": 0.20,
    "spliceai-38-xwkwwwxdwq-uc.a.run.app": 0.20,
    "spliceai-37-xwkwwwxdwq-uc.a.run.app": 0.20,
}
_DEFAULT_MIN_INTERVAL = 0.20

try:
    _MAX_CONCURRENT = max(1, int(os.environ.get("HEARTVAR_MAX_OUTBOUND_CONCURRENCY", "30")))
except ValueError:
    _MAX_CONCURRENT = 30

_host_last_call: dict[str, float] = {}
_host_locks: dict[str, asyncio.Lock] = {}
_global_sem: asyncio.Semaphore | None = None


def global_semaphore() -> asyncio.Semaphore:
    """Lazily create the global concurrency semaphore (on the running loop)."""
    global _global_sem
    if _global_sem is None:
        _global_sem = asyncio.Semaphore(_MAX_CONCURRENT)
    return _global_sem


async def throttle_host(url: str) -> None:
    """Block until it is this host's turn under its per-host min-interval.

    Serialises calls to a single host so the process-wide outbound rate to it
    never exceeds 1/interval, no matter how many inbound users there are. A
    single curation makes ~1 call per host, so it waits ~0; the throttle only
    queues when concurrent curations target the same host."""
    host = urlparse(url).netloc
    if not host:
        return
    interval = _HOST_MIN_INTERVAL.get(host, _DEFAULT_MIN_INTERVAL)
    lock = _host_locks.setdefault(host, asyncio.Lock())
    async with lock:
        wait = (_host_last_call.get(host, 0.0) + interval) - time.monotonic()
        if wait > 0:
            await asyncio.sleep(wait)
        _host_last_call[host] = time.monotonic()
