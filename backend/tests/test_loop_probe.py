"""Tests for the event-loop probe — the instrument for the ~19 s production stall.

WHY THIS EXISTS. Browser-side profiling on 2026-08-30 caught fifteen sources
completing in the SAME MILLISECOND at ~19 s. That is the signature of a BLOCKED
EVENT LOOP, not of fifteen slow lookups: `pubmed`/`pubtator3` are pure network
and landed in the same 5 ms window as the local-SQLite sources and the VEP
subprocess, while `erepo` (also network) finished cleanly at 204 ms. A network
completion can only be dragged out like that by a loop that is not running
callbacks.

Profiling then ELIMINATED the written-up hypothesis (~10 synchronous
json.load/read_text calls): measured inside the real image, the whole cold path
is ~250 ms of CPU (gencc's 11.3 MB read_text+parse is 23.5 ms; hgvs+cdot+pysam
imports total 122 ms; `hgvs.parser.Parser()` is 1 ms). So the block is I/O
LATENCY on the SMB mount, on a path that runs on the loop — and two candidates
remain (`hgvs_resolver.resolve` called bare inside `async fetch_hgvs`, and
gencc's `read_text` inside `async _ensure_loaded`). Deciding between them needs
the loop itself to name the callback.

asyncio already does exactly that: with debug mode ON it logs
``Executing <Handle ...> took N seconds`` for any callback exceeding
``slow_callback_duration``. Stdlib, no new dependency, no new log plumbing.

The four contracts:
  1. Default OFF — with ``HEARTVAR_LOOP_PROBE`` unset the loop is NOT touched,
     so production behaviour is byte-identical (debug mode costs real overhead
     on every callback, so it must never be on by accident).
  2. Enabled — sets loop debug AND ``slow_callback_duration``. Both are
     required: asyncio gates the slow-callback check on ``get_debug()``, so the
     threshold alone silently does nothing.
  3. The threshold is env-tunable and falls back safely on junk input.
  4. The ``asyncio`` logger can actually emit WARNING, so the lines reach the
     CONTAINER LOG STREAM. This is the trap the offline-VEP banner already fell
     into — ``offline_status()`` writes to a file on the mount, so IT reading
     the log stream never saw it. An instrument nobody can read is not one.

Fully offline; no sleeps, no network. Runnable with pytest or directly
(``python -m backend.tests.test_loop_probe``).
"""
from __future__ import annotations

import asyncio
import logging

import pytest

from backend import loop_probe


@pytest.fixture(autouse=True)
def _clean_logging_state():
    """Strip the probe's asyncio forwarder and restore the logger level.

    Arming the probe mutates PROCESS-WIDE logging state that outlives the test.
    Without this, a test that arms the probe leaves the forwarder installed and
    the next "probe is OFF" assertion passes or fails for the wrong reason — so
    the fixture exists to keep those assertions honest, not to soften them.
    """
    asyncio_log = logging.getLogger("asyncio")
    saved_level = asyncio_log.level
    yield
    for h in list(asyncio_log.handlers):
        if isinstance(h, loop_probe._ForwardToHeartvar):
            asyncio_log.removeHandler(h)
    asyncio_log.setLevel(saved_level)


def _in_loop(fn):
    """Run ``fn(loop)`` inside a fresh event loop and return its result.

    A fresh loop per test matters: ``set_debug`` is loop state, so a leaked
    debug-mode loop would make a later "default OFF" assertion pass for the
    wrong reason.
    """
    async def _main():
        return fn(asyncio.get_running_loop())
    return asyncio.run(_main())


def test_probe_disabled_by_default_leaves_the_loop_untouched(monkeypatch):
    monkeypatch.delenv("HEARTVAR_LOOP_PROBE", raising=False)

    def check(loop):
        before = loop.slow_callback_duration
        enabled = loop_probe.maybe_enable_loop_probe()
        return enabled, loop.get_debug(), loop.slow_callback_duration, before

    enabled, debug, after, before = _in_loop(check)
    assert enabled is False
    assert debug is False, "debug mode must not be enabled by default"
    assert after == before, "slow_callback_duration must not be touched"


def test_probe_off_for_falsey_values(monkeypatch):
    for raw in ("", "0", "false", "no", "off", "  "):
        monkeypatch.setenv("HEARTVAR_LOOP_PROBE", raw)
        assert loop_probe.probe_enabled() is False, f"{raw!r} should be OFF"


def test_probe_enables_debug_and_threshold_together(monkeypatch):
    monkeypatch.setenv("HEARTVAR_LOOP_PROBE", "1")
    monkeypatch.delenv("HEARTVAR_LOOP_PROBE_THRESHOLD", raising=False)

    def check(loop):
        enabled = loop_probe.maybe_enable_loop_probe()
        return enabled, loop.get_debug(), loop.slow_callback_duration

    enabled, debug, threshold = _in_loop(check)
    assert enabled is True
    assert debug is True, "slow-callback logging is gated on get_debug()"
    assert threshold == 0.25


def test_probe_on_for_truthy_values(monkeypatch):
    for raw in ("1", "true", "TRUE", "yes", "on", " On "):
        monkeypatch.setenv("HEARTVAR_LOOP_PROBE", raw)
        assert loop_probe.probe_enabled() is True, f"{raw!r} should be ON"


def test_threshold_is_env_tunable(monkeypatch):
    monkeypatch.setenv("HEARTVAR_LOOP_PROBE", "1")
    monkeypatch.setenv("HEARTVAR_LOOP_PROBE_THRESHOLD", "1.5")
    assert _in_loop(lambda loop: (loop_probe.maybe_enable_loop_probe(),
                                  loop.slow_callback_duration))[1] == 1.5


def test_threshold_falls_back_on_junk(monkeypatch):
    monkeypatch.setenv("HEARTVAR_LOOP_PROBE", "1")
    monkeypatch.setenv("HEARTVAR_LOOP_PROBE_THRESHOLD", "not-a-number")
    assert _in_loop(lambda loop: (loop_probe.maybe_enable_loop_probe(),
                                  loop.slow_callback_duration))[1] == 0.25


def test_asyncio_logger_can_emit_warning_when_probe_is_on(monkeypatch):
    """The whole point is a line in the CONTAINER LOG STREAM."""
    monkeypatch.setenv("HEARTVAR_LOOP_PROBE", "1")
    logging.getLogger("asyncio").setLevel(logging.CRITICAL)

    _in_loop(lambda loop: loop_probe.maybe_enable_loop_probe())

    assert logging.getLogger("asyncio").isEnabledFor(logging.WARNING), (
        "asyncio logs the slow-callback line at WARNING; if that logger is "
        "above WARNING the instrument is silent"
    )


def test_asyncio_logger_untouched_when_probe_is_off(monkeypatch):
    monkeypatch.delenv("HEARTVAR_LOOP_PROBE", raising=False)
    logging.getLogger("asyncio").setLevel(logging.CRITICAL)

    _in_loop(lambda loop: loop_probe.maybe_enable_loop_probe())

    assert not logging.getLogger("asyncio").isEnabledFor(logging.WARNING), (
        "with the probe OFF nothing may reconfigure logging"
    )


def _lifespan_source() -> str:
    """``app.py``'s ``lifespan`` source, read from disk.

    Read as TEXT rather than imported: ``backend.app`` pulls fastapi, anthropic,
    pysam and the rest of the runtime, and this contract needs none of it. The
    repo already checks call-graph facts this way (see
    ``test_gencc_local_only.test_the_builder_still_owns_the_refresh``).
    """
    import ast
    from pathlib import Path
    app_py = Path(__file__).resolve().parents[1] / "app.py"
    tree = ast.parse(app_py.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.AsyncFunctionDef) and node.name == "lifespan":
            return ast.unparse(node)
    raise AssertionError("no `async def lifespan` found in backend/app.py")


def test_lifespan_arms_the_probe():
    """An instrument that is never called cannot report anything.

    The probe has to be armed from the lifespan specifically, because
    ``maybe_enable_loop_probe`` needs a RUNNING loop — arming it at import time
    would hit the no-running-loop branch and silently do nothing.
    """
    import ast
    tree = ast.parse(_lifespan_source())
    called = {
        n.func.id for n in ast.walk(tree)
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
    } | {
        n.func.attr for n in ast.walk(tree)
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
    }
    assert "maybe_enable_loop_probe" in called, (
        "backend/app.py's lifespan must call maybe_enable_loop_probe() so "
        "HEARTVAR_LOOP_PROBE has somewhere to take effect"
    )


def test_app_imports_the_probe():
    """The lifespan call has to resolve to this module, not to a stray name."""
    from pathlib import Path
    app_py = (Path(__file__).resolve().parents[1] / "app.py").read_text(
        encoding="utf-8")
    assert "from .loop_probe import maybe_enable_loop_probe" in app_py


def test_asyncio_records_reach_the_heartvar_logger_tree(monkeypatch):
    """The slow-callback line must land in the DEPLOYMENT LOG FILE.

    ⚠ THE TRAP THIS EXISTS FOR. ``_record_deploy_event`` attaches its
    ``FileHandler`` to ``logging.getLogger("heartvar")`` (app.py), and
    /api/admin/log-content serves ONLY that file — categories `db_builder` and
    `deployments`, from HEARTVAR_LOGS_DIR. It does NOT serve container stdout.
    We have no Azure access (deploy doc item A0 is still an open ask to IT), so
    the admin page is the ONLY readable surface.

    asyncio logs the verdict on the ``asyncio`` logger, which is NOT under
    ``heartvar`` and therefore never reaches that file. Arming the probe without
    this forwarding produces a diagnostic whose output nobody can read — the
    exact failure mode of the offline-VEP banner.
    """
    monkeypatch.setenv("HEARTVAR_LOOP_PROBE", "1")
    seen: list[str] = []

    class _Capture(logging.Handler):
        def emit(self, record):
            seen.append(record.getMessage())

    handler = _Capture()
    heartvar_log = logging.getLogger("heartvar")
    heartvar_log.addHandler(handler)
    try:
        _in_loop(lambda loop: loop_probe.maybe_enable_loop_probe())
        logging.getLogger("asyncio").warning(
            "Executing <Task ... coro=<fetch_vep() defined at "
            "backend/clients/ensembl_vep.py:1> ...> took 18.900 seconds")
    finally:
        heartvar_log.removeHandler(handler)

    assert any("took 18.900 seconds" in m for m in seen), (
        "asyncio's slow-callback verdict never reached the `heartvar` logger, so "
        "it will not be in the deployment log file the admin page serves"
    )


def test_asyncio_forwarding_not_installed_when_probe_is_off(monkeypatch):
    monkeypatch.delenv("HEARTVAR_LOOP_PROBE", raising=False)
    seen: list[str] = []

    class _Capture(logging.Handler):
        def emit(self, record):
            seen.append(record.getMessage())

    handler = _Capture()
    heartvar_log = logging.getLogger("heartvar")
    heartvar_log.addHandler(handler)
    try:
        _in_loop(lambda loop: loop_probe.maybe_enable_loop_probe())
        logging.getLogger("asyncio").warning("should not be forwarded")
    finally:
        heartvar_log.removeHandler(handler)

    assert not any("should not be forwarded" in m for m in seen), (
        "with the probe OFF nothing may re-route asyncio's records"
    )


if __name__ == "__main__":  # pragma: no cover
    import pytest as _pytest
    raise SystemExit(_pytest.main([__file__, "-q"]))
