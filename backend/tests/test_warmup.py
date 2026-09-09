"""Tests for start-up warming of the local data mount.

WHY. Measured 2026-08-30 on a cold production container: the RYR2 curation took
42.6 s, of which ~20 s was the DB gather and ~15.9 s was `same_site`. On a WARM
container the same curation's `vep` source drops from 18.58 s to 273 ms. So a
large part of the cost is FIRST ACCESS to files on the SMB mount, and every
container restart — i.e. every deploy — hands that bill to whichever curator
happens to click first.

Warming moves it to start-up, where nobody is waiting. The lifespan already does
exactly this for the PanelApp HPO descendant cache ("so the first real curation
doesn't pay the ~9 s cold-start cost"); this generalises it to the mount.

It also LOGS a per-file timing, which makes the next deployment log say which
file is slow — readable in the admin page, without arming the loop probe.

Contracts:
  1. Every SQLite database under the data dir gets touched.
  2. A missing, corrupt or unreadable file NEVER raises — warming is
     best-effort and must not stop the server from starting.
  3. It runs on the local-I/O pool, not the event loop and not the default pool.
  4. It can be turned off.
"""
from __future__ import annotations

import asyncio
import sqlite3
import threading

import pytest

import time

from backend import localio, warmup


def _make_db(path, table="t"):
    con = sqlite3.connect(path)
    con.execute(f"CREATE TABLE {table} (k TEXT PRIMARY KEY)")
    con.execute(f"INSERT INTO {table} VALUES ('x')")
    con.commit()
    con.close()


def test_warms_every_database_it_finds(tmp_path):
    _make_db(tmp_path / "a.db")
    _make_db(tmp_path / "b.db")
    report = asyncio.run(warmup.warm_local_data(tmp_path))
    assert set(report["warmed"]) == {"a.db", "b.db"}
    assert all(isinstance(v, float) for v in report["seconds"].values())


def test_a_corrupt_database_does_not_raise(tmp_path):
    _make_db(tmp_path / "good.db")
    (tmp_path / "bad.db").write_bytes(b"this is not a database at all")
    report = asyncio.run(warmup.warm_local_data(tmp_path))
    assert "good.db" in report["warmed"]
    assert "bad.db" in report["failed"], "a corrupt file must be reported, not raised"


def test_a_missing_directory_does_not_raise(tmp_path):
    report = asyncio.run(warmup.warm_local_data(tmp_path / "nope"))
    assert report["warmed"] == []


def test_warming_runs_on_the_local_io_pool(tmp_path):
    """It must not occupy the default pool — that is what starves DNS.

    Spies on _warm_one_db rather than on localio.connect_ro. Warming no longer
    uses connect_ro (it opens and CLOSES, so it cannot leave a handle that
    denies the data build's rename — see warmup._warm_one_db), and the pool this
    test is actually about is the one run_local dispatches to.
    """
    _make_db(tmp_path / "a.db")
    seen: list[str] = []
    real = warmup._warm_one_db

    def spy(path):
        seen.append(threading.current_thread().name)
        return real(path)

    warmup._warm_one_db = spy
    try:
        asyncio.run(warmup.warm_local_data(tmp_path))
    finally:
        warmup._warm_one_db = real
    assert seen, "_warm_one_db was never called"
    assert all(localio.THREAD_NAME_PREFIX in n for n in seen), (
        f"warming ran on the wrong pool: {seen}"
    )


def test_warming_leaves_no_open_connection(tmp_path):
    """The bug that failed a production data build. scripts/_dbbuild.publish
    renames the new database over the old one, and on Azure Files SMB an open
    handle on the target DENIES that rename:
        PermissionError: [Errno 13] Permission denied:
          '/app/data/clinvar.db.stage.33' -> '/app/data/clinvar.db'
    connect_ro caches per (thread, database) and close_all() has no production
    caller, so warming used to leave a handle open on every database for the
    life of the process."""
    _make_db(tmp_path / "a.db")
    localio.close_all()
    asyncio.run(warmup.warm_local_data(tmp_path))
    cache = getattr(localio._local, "conns", None) or {}
    assert not any("a.db" in path for path in cache), (
        "warming left a cached connection open — a data build's atomic rename "
        f"onto this file would be denied: {sorted(cache)}")


def test_can_be_disabled(monkeypatch):
    monkeypatch.setenv("HEARTVAR_WARM_LOCAL_DATA", "0")
    assert warmup.warming_enabled() is False


def test_enabled_by_default(monkeypatch):
    monkeypatch.delenv("HEARTVAR_WARM_LOCAL_DATA", raising=False)
    assert warmup.warming_enabled() is True, (
        "warming is the fix for the cold-start cost — it must be ON unless "
        "explicitly disabled"
    )


def test_lifespan_schedules_the_warmup():
    """An unused warmer warms nothing."""
    import ast
    from pathlib import Path
    tree = ast.parse((Path(__file__).resolve().parents[1] / "app.py")
                     .read_text(encoding="utf-8"))
    src = next(ast.unparse(n) for n in ast.walk(tree)
               if isinstance(n, ast.AsyncFunctionDef) and n.name == "lifespan")
    assert "warm_local_data" in src, (
        "backend/app.py's lifespan must schedule warm_local_data()"
    )


def test_warming_is_concurrency_bounded(tmp_path):
    """⚠ MEASURED REGRESSION, 2026-08-30. Unbounded, warming read every DB and
    the 847 MB FASTA off the mount at once — so a curation arriving DURING
    start-up contended with it and took 31.7 s, against 5.4 s once warming was
    done. Warming exists to protect the first curation, so it must not be the
    thing that slows it down.

    Slower to finish, far less destructive while it runs.
    """
    for n in "abcdefgh":
        _make_db(tmp_path / f"{n}.db")

    lock = threading.Lock()
    state = {"now": 0, "peak": 0}
    real = warmup._warm_one_db

    def counting(path):
        with lock:
            state["now"] += 1
            state["peak"] = max(state["peak"], state["now"])
        try:
            time.sleep(0.05)
            return real(path)
        finally:
            with lock:
                state["now"] -= 1

    warmup._warm_one_db = counting
    try:
        report = asyncio.run(warmup.warm_local_data(tmp_path))
    finally:
        warmup._warm_one_db = real

    assert len(report["warmed"]) == 8, "every file must still be warmed"
    assert state["peak"] <= warmup.warm_concurrency(), (
        f"peak concurrency {state['peak']} exceeded the "
        f"{warmup.warm_concurrency()} slot limit — warming is unthrottled"
    )


def test_warm_concurrency_defaults_low_and_is_tunable(monkeypatch):
    monkeypatch.delenv("HEARTVAR_WARM_CONCURRENCY", raising=False)
    assert warmup.warm_concurrency() == 2
    monkeypatch.setenv("HEARTVAR_WARM_CONCURRENCY", "5")
    assert warmup.warm_concurrency() == 5


def test_warm_concurrency_falls_back_on_junk(monkeypatch):
    for raw in ("nonsense", "0", "-3", ""):
        monkeypatch.setenv("HEARTVAR_WARM_CONCURRENCY", raw)
        assert warmup.warm_concurrency() == 2, f"{raw!r} should fall back"


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
