"""Tests for the dedicated local-I/O executor and the read-only SQLite cache.

WHY THIS MODULE EXISTS — measured, not assumed. The loop-probe reading of
2026-08-30 (RYR2, cold container) reproduced a 20.2 s stall with 15 sources
completing inside a 160 ms window, while the LARGEST single loop block in the
whole run was 0.533 s. The loop was never blocked; the work was queued off it.

~14 clients dispatch ``asyncio.to_thread`` for SQLite over the SMB mount, each
opening a FRESH ``sqlite3.connect`` per call, into the DEFAULT pool of
``min(32, cpu_count + 4)``. asyncio also routes every ``getaddrinfo`` through
that SAME default pool — so once the DB work fills it, DNS for the outbound HTTP
clients queues behind it, which is why ``pubmed``/``pubtator3``/
``gene_literature`` flushed together with the DB sources despite being pure
network.

Two fixes, two contracts:
  1. Local data access gets its OWN executor, so saturating it can never starve
     DNS or anything else on the default pool. The third test here is the actual
     regression test for the bug.
  2. Read-only connections are cached PER THREAD and reused, so a call stops
     paying `sqlite3.connect` (open + header + schema read = several SMB round
     trips) every single time.
"""
from __future__ import annotations

import asyncio
import os
import sqlite3
import threading
import time

import pytest

from backend import localio


def test_run_local_returns_the_result():
    assert asyncio.run(localio.run_local(lambda a, b: a + b, 2, 3)) == 5


def test_run_local_propagates_exceptions():
    def boom():
        raise ValueError("propagated")

    with pytest.raises(ValueError, match="propagated"):
        asyncio.run(localio.run_local(boom))


def test_run_local_passes_keyword_arguments():
    assert asyncio.run(localio.run_local(lambda *, x: x * 2, x=21)) == 42


def test_run_local_does_not_run_on_the_loop_thread():
    loop_thread = threading.get_ident()
    worker = asyncio.run(localio.run_local(threading.get_ident))
    assert worker != loop_thread


def test_run_local_uses_the_dedicated_pool_not_the_default():
    """Thread name is the observable that tells the two pools apart."""
    name = asyncio.run(localio.run_local(lambda: threading.current_thread().name))
    assert localio.THREAD_NAME_PREFIX in name, (
        f"expected a {localio.THREAD_NAME_PREFIX!r} thread, got {name!r} — "
        "local I/O is still landing on the default executor"
    )


def test_saturating_local_io_does_not_starve_the_default_executor():
    """Fill the local-I/O pool; default-pool work must still complete promptly.

    This is the bug, in miniature. Before the split, these shared one pool: the
    SQLite-over-SMB calls filled it and asyncio's getaddrinfo — and therefore
    every outbound HTTP source — waited behind them. ``asyncio.to_thread`` stands
    in for getaddrinfo here because both dispatch to the DEFAULT executor.
    """
    async def main():
        release = threading.Event()
        started = threading.Barrier(2)

        def hog():
            with contextlib_suppress():
                started.wait(timeout=5)
            release.wait(timeout=10)
            return "hog"

        hogs = [asyncio.create_task(localio.run_local(hog))
                for _ in range(localio.max_workers() * 2)]
        await asyncio.sleep(0.2)

        t0 = time.perf_counter()
        got = await asyncio.wait_for(asyncio.to_thread(lambda: "dns"), timeout=5)
        elapsed = time.perf_counter() - t0

        release.set()
        await asyncio.gather(*hogs, return_exceptions=True)
        return got, elapsed

    got, elapsed = asyncio.run(main())
    assert got == "dns"
    assert elapsed < 2.0, (
        f"default-executor work took {elapsed:.2f}s while local I/O was "
        "saturated — the pools are not actually separated"
    )


import contextlib


def contextlib_suppress():
    return contextlib.suppress(Exception)


@pytest.fixture
def db(tmp_path):
    path = tmp_path / "t.db"
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE t (k TEXT PRIMARY KEY, v TEXT)")
    con.execute("INSERT INTO t VALUES ('a', '1')")
    con.commit()
    con.close()
    return str(path)


def test_connect_ro_reads(db):
    con = localio.connect_ro(db)
    assert con.execute("SELECT v FROM t WHERE k='a'").fetchone()[0] == "1"


def test_connect_ro_reuses_the_same_connection_on_one_thread(db):
    """The whole point: stop paying connect() per call over SMB."""
    assert localio.connect_ro(db) is localio.connect_ro(db)


def test_connect_ro_is_per_thread(db):
    """sqlite3 connections are not safe to share across threads by default."""
    seen = {}

    def grab(tag):
        seen[tag] = localio.connect_ro(db)

    t1 = threading.Thread(target=grab, args=("a",))
    t2 = threading.Thread(target=grab, args=("b",))
    t1.start(); t1.join(); t2.start(); t2.join()
    assert seen["a"] is not seen["b"], "one connection was shared across threads"


def test_connect_ro_is_read_only(db):
    con = localio.connect_ro(db)
    with pytest.raises(sqlite3.OperationalError):
        con.execute("INSERT INTO t VALUES ('b', '2')")


def test_connect_ro_recovers_from_a_closed_connection(db):
    """A cached handle that has gone bad must be replaced, not re-served."""
    first = localio.connect_ro(db)
    first.close()
    second = localio.connect_ro(db)
    assert second is not first
    assert second.execute("SELECT v FROM t WHERE k='a'").fetchone()[0] == "1"


def test_connect_ro_returns_none_for_a_missing_file(tmp_path):
    assert localio.connect_ro(str(tmp_path / "nope.db")) is None


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))


def _one_row_db(path, value):
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE t (v TEXT)")
    con.execute("INSERT INTO t VALUES (?)", (value,))
    con.commit()
    con.close()


def test_an_atomically_replaced_database_is_reopened(tmp_path):
    live = tmp_path / "live.db"
    _one_row_db(live, "old")
    localio.close_all()
    try:
        first = localio.connect_ro(live)
        assert first.execute("SELECT v FROM t").fetchone()[0] == "old"

        staged = tmp_path / "live.db.stage.1234"
        _one_row_db(staged, "new")
        os.replace(staged, live)

        second = localio.connect_ro(live)
        assert second.execute("SELECT v FROM t").fetchone()[0] == "new", (
            "still reading the pre-rename inode — a data build would be "
            "invisible until the app restarted")
    finally:
        localio.close_all()


def test_an_unchanged_database_keeps_the_same_connection(tmp_path):
    """The staleness check must not throw away the page cache on every call —
    that would undo the whole point of connect_ro."""
    live = tmp_path / "stable.db"
    _one_row_db(live, "x")
    localio.close_all()
    try:
        a = localio.connect_ro(live)
        b = localio.connect_ro(live)
        assert a is b, "reopened an unchanged file — page cache thrown away"
    finally:
        localio.close_all()
