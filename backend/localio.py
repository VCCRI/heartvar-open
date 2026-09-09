"""A dedicated executor for LOCAL DATA ACCESS, and reusable read-only SQLite.

WHY, AND IT WAS MEASURED RATHER THAN REASONED. The loop-probe reading of
2026-08-30 (RYR2, cold container) reproduced a 20.2 s stall in which 15 sources
completed inside a 160 ms window — while the LARGEST single event-loop block in
the whole run was 0.533 s. So the loop was never blocked. The work was queued
somewhere else, and that somewhere is the DEFAULT thread pool.

Roughly fourteen clients dispatch ``asyncio.to_thread`` for SQLite / tabix /
bigWig / file reads against the SMB mount, and each one opened a FRESH
``sqlite3.connect`` per call. ``asyncio.to_thread`` uses the loop's default
executor, sized ``min(32, os.cpu_count() + 4)``. Two consequences, both bad:

  1. The DB work serialises against itself once the pool is full.
  2. **asyncio routes every ``getaddrinfo`` through that same default executor.**
     So the moment SQLite-over-SMB fills it, DNS for every outbound HTTP client
     queues behind it — which is exactly why ``pubmed``, ``pubtator3`` and
     ``gene_literature`` flushed together with the DB sources despite being pure
     network calls that never touch the mount.

This module fixes both halves:

  * ``run_local`` is a drop-in for ``asyncio.to_thread`` that dispatches to a
    SEPARATE, larger pool. Saturating local I/O can then never starve DNS.
  * ``connect_ro`` caches one read-only connection PER THREAD PER DATABASE.
    Because the pool's threads are long-lived, a call stops paying
    ``sqlite3.connect`` — open, header read, schema read, i.e. several SMB round
    trips — on every single query.

SIZING. The work is I/O-bound on a network filesystem, not CPU-bound, so the
right worker count tracks concurrent sources rather than cores. The default is
32: comfortably above the ~14 concurrent DB sources plus their follow-up
queries, and each idle thread costs only its stack.

⚠ READ-ONLY, ALWAYS, and that is a correctness requirement rather than caution:
SQLite cannot take write locks on Azure Files. Every caller here is a reader —
the databases are built by the monthly job, never by the webapp.

⚠ STALENESS. A cached connection holds its file handle open, so if the monthly
job replaces a database underneath a running container, that container keeps
serving the OLD file until it restarts. That is the safe direction (a torn read
would be worse) and ``monthly-restart-webapp.yml`` restarts the app after the
build, so the window closes on its own.
"""
from __future__ import annotations

import asyncio
import functools
import logging
import os
import sqlite3
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, TypeVar

log = logging.getLogger("heartvar.localio")

T = TypeVar("T")

THREAD_NAME_PREFIX = "heartvar-localio"

_DEFAULT_MAX_WORKERS = 32

_pool: ThreadPoolExecutor | None = None
_pool_lock = threading.Lock()


def log_safe(value: object) -> str:
    """Neutralise user-controlled text before it enters a log record.

    Strips CR/LF/tab so a crafted gene / HGVS / chat question cannot forge or
    split log lines (CWE-117 log injection). Applied to request-derived values
    at every log call that interpolates them.

    Lives here rather than in app.py because the CLIENTS need it too — CodeQL
    alert #60 was a `gene` reaching a log line in clients/clinvar.py, which
    cannot import app.py. app.py's `_log_safe` is now an alias for this, so
    there is one implementation rather than two free to drift.
    """
    return (
        str(value)
        .replace("\r\n", " ")
        .replace("\n", " ")
        .replace("\r", " ")
        .replace("\t", " ")
    )


def max_workers() -> int:
    """Worker count for the local-I/O pool.

    Deliberately NOT derived from ``cpu_count``: these threads sit in SMB reads,
    so the useful count tracks concurrent sources, not cores. Env-tunable, with
    the codebase's usual fall-back-on-junk behaviour.
    """
    raw = os.environ.get("HEARTVAR_LOCALIO_WORKERS", "").strip()
    if raw:
        try:
            value = int(raw)
            if value > 0:
                return value
            log.warning("HEARTVAR_LOCALIO_WORKERS=%r must be > 0 — using %d",
                        raw, _DEFAULT_MAX_WORKERS)
        except ValueError:
            log.warning("HEARTVAR_LOCALIO_WORKERS=%r is not an int — using %d",
                        raw, _DEFAULT_MAX_WORKERS)
    return _DEFAULT_MAX_WORKERS


def executor() -> ThreadPoolExecutor:
    """The shared local-I/O pool, created on first use.

    Lazy and double-checked rather than built at import: importing this module
    must not spawn threads in a test process or a CLI script that never does any
    local I/O.
    """
    global _pool
    if _pool is None:
        with _pool_lock:
            if _pool is None:
                workers = max_workers()
                _pool = ThreadPoolExecutor(
                    max_workers=workers,
                    thread_name_prefix=THREAD_NAME_PREFIX,
                )
                log.info("[localio] dedicated local-I/O pool: %d workers "
                         "(the default executor is left free for getaddrinfo)",
                         workers)
    return _pool


async def run_local(fn: Callable[..., T], /, *args: Any, **kwargs: Any) -> T:
    """Run ``fn`` on the local-I/O pool. Drop-in for ``asyncio.to_thread``.

    Same signature and same semantics — result returned, exception propagated —
    so migrating a call site is a rename. The ONLY difference is which pool it
    lands on, and that difference is the entire point: work sent here cannot
    delay the DNS lookups asyncio performs on the default executor.

    NOTE: unlike ``asyncio.to_thread`` this does not copy the caller's
    ``contextvars`` context. No caller relies on it — these are file and SQLite
    reads that take their arguments explicitly.
    """
    loop = asyncio.get_running_loop()
    return await loop.run_in_executor(
        executor(), functools.partial(fn, *args, **kwargs)
    )


_local = threading.local()


def connect_ro(db_path: str | os.PathLike[str]) -> sqlite3.Connection | None:
    """A cached read-only connection to ``db_path``, or None if unusable.

    Returns None rather than raising for a missing file, matching how the
    clients already treat an absent database (serve nothing, let the caller fall
    back) instead of turning it into a 500.
    """
    path = os.fspath(db_path)
    cache: dict[str, tuple[sqlite3.Connection, tuple[int, int] | None]] = getattr(
        _local, "conns", None)
    if cache is None:
        cache = _local.conns = {}

    entry = cache.get(path)
    if entry is not None:
        con, cached_id = entry
        if _file_identity(path) != cached_id:
            log.info("[localio] %s was replaced on disk — reopening", path)
            cache.pop(path, None)
            try:
                con.close()
            except sqlite3.Error:
                pass
        else:
            try:
                con.execute("SELECT 1").fetchone()
                return con
            except sqlite3.Error:
                cache.pop(path, None)
                try:
                    con.close()
                except sqlite3.Error:
                    pass

    if not os.path.isfile(path):
        return None
    try:
        con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        _tune_for_slow_storage(con, path)
    except sqlite3.Error as exc:
        log.warning("[localio] cannot open %r read-only (%r)", path, exc)
        return None
    cache[path] = (con, _file_identity(path))
    return con


def _file_identity(path: str) -> tuple[int, int] | None:
    """(device, inode) for ``path``, or None if it cannot be stat'd.

    This is what detects an atomically-published replacement: the data build
    renames a new file over the old one, so the path is unchanged and only the
    inode moves. Returning None on a failed stat means "unknown", which
    compares unequal to any recorded identity and so forces a reopen — the safe
    direction, since the alternative is serving stale data forever.
    """
    try:
        st = os.stat(path)
    except OSError:
        return None
    return (st.st_dev, st.st_ino)


_PAGE_CACHE_KIB = int(os.environ.get("HEARTVAR_SQLITE_CACHE_KIB", "8192"))


def _tune_for_slow_storage(con: sqlite3.Connection, path: str) -> None:
    """Make a read-only connection cheap to re-use over a network mount.

    Read-only and idempotent: every PRAGMA here is a per-connection hint, none
    writes to the database (which matters — the mount is mounted read-only and
    we have no write access to it).
    """
    try:
        con.execute(f"PRAGMA cache_size=-{_PAGE_CACHE_KIB}")
        con.execute("PRAGMA mmap_size=268435456")
        con.execute("PRAGMA read_uncommitted=1")
    except sqlite3.Error as exc:
        log.debug("[localio] could not tune %r (%r)", path, exc)


def close_all() -> None:
    """Close this thread's cached connections. For tests and shutdown."""
    cache: dict[str, tuple[sqlite3.Connection, tuple[int, int] | None]] = getattr(
        _local, "conns", None) or {}
    for con, _identity in cache.values():
        try:
            con.close()
        except sqlite3.Error:
            pass
    cache.clear()
