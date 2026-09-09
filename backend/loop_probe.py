"""Event-loop probe — the instrument for the ~19 s production stall.

WHAT IT IS FOR. Browser-side profiling of production (2026-08-30) caught fifteen
sources completing in the SAME MILLISECOND at ~19 s. That is a BLOCKED EVENT
LOOP, not fifteen slow lookups: ``pubmed``/``pubtator3`` are pure network and
landed in the same 5 ms window as the local-SQLite sources and the VEP
subprocess, while ``erepo`` (also network) finished cleanly at 204 ms. A network
completion can only be dragged out that way by a loop that is not running
callbacks.

WHAT WAS ALREADY RULED OUT, so nobody re-runs it. The written-up hypothesis was
"~10 synchronous json.load/read_text calls over SMB". Measured inside the real
image, that whole cold path is ~250 ms of CPU:

    gencc 11.3 MB read_text + json.loads      23.5 ms
    hgnc_alias 1.5 MB                          7.4 ms
    hpo_labels / gene_id_map.json.gz      2.6 / 9.6 ms
    hgvs + cdot + pysam imports              122 ms  (IPython/psycopg2 never load)
    hgvs.parser.Parser()                       1 ms
    cdot make_provider + AssemblyMapper      104 ms
    one c.->g. resolve                         1.8 ms

So the block is I/O LATENCY on the SMB mount, on a path that runs on the loop.
Two candidates remain, both in the 19 s flush group and both blocking:

  1. ``clients/vep_offline.py:808`` — ``hgvs_resolver.resolve(...)`` called BARE
     inside ``async def fetch_hgvs``. The only candidate that is NEW:
     ``HEARTVAR_HGVS_RESOLVER=1`` was flipped on 2026-08-29/30.
  2. ``clients/gencc.py:136`` — 11.3 MB ``read_text()`` inside
     ``async def _ensure_loaded``. Reading from the mount predates the stall.

HOW IT WORKS. asyncio already names the blocker: with debug mode ON it logs
``Executing <Handle ... created at FILE:LINE> took N seconds`` for any callback
over ``slow_callback_duration``, and in debug mode the handle carries its
creation site — which is the file and line we are trying to identify. Stdlib, no
new dependency, nothing to parse.

DEFAULT OFF, and that matters more than usual: debug mode wraps EVERY callback
with timing and capture, so it is real per-callback overhead, not a free flag. It
is a diagnostic to be switched on for a run and switched off again, not a
setting to leave on.

⚠ IT LOGS TO THE CONTAINER LOG STREAM, deliberately. ``offline_status()`` writes
its banner to a file on the mount inside a ``try/except OSError``, so IT reading
the container log stream never saw whether offline VEP was on. An instrument
nobody can read is not an instrument, so this one only ever touches stdlib
logging and asserts the ``asyncio`` logger can actually emit WARNING.

USAGE. Set ``HEARTVAR_LOOP_PROBE=1``, restart, run ONE curation, read the log
for ``took`` — the largest line names the blocking callback. Then unset it.
``HEARTVAR_LOOP_PROBE_THRESHOLD`` (seconds, default 0.25) sets the floor.
"""
from __future__ import annotations

import asyncio
import logging
import os

log = logging.getLogger("heartvar.loop_probe")

_DEFAULT_THRESHOLD = 0.25


def _truthy(name: str) -> bool:
    """A permissive truthy env check — same semantics as ``prewarm._truthy``.
    Unset, blank or a falsey word -> False, so the probe is OFF by default."""
    return os.environ.get(name, "").strip().lower() in {"1", "true", "yes", "on"}


def probe_enabled() -> bool:
    """Master switch. Default OFF. Read live (not cached at import) so a deploy
    can flip it via the environment without a code change."""
    return _truthy("HEARTVAR_LOOP_PROBE")


def _threshold() -> float:
    """Seconds a callback may hold the loop before asyncio logs it.

    Falls back to the default on unset/blank/junk, following the codebase's
    ``_float_env`` convention — a typo in an env var must not stop the server.
    A non-positive value would log every callback, so it is also rejected.
    """
    raw = os.environ.get("HEARTVAR_LOOP_PROBE_THRESHOLD", "").strip()
    if not raw:
        return _DEFAULT_THRESHOLD
    try:
        value = float(raw)
    except ValueError:
        log.warning("HEARTVAR_LOOP_PROBE_THRESHOLD=%r is not a number — using "
                    "%.2fs", raw, _DEFAULT_THRESHOLD)
        return _DEFAULT_THRESHOLD
    if value <= 0:
        log.warning("HEARTVAR_LOOP_PROBE_THRESHOLD=%r must be > 0 — using "
                    "%.2fs", raw, _DEFAULT_THRESHOLD)
        return _DEFAULT_THRESHOLD
    return value


class _ForwardToHeartvar(logging.Handler):
    """Re-dispatch an ``asyncio`` record through the ``heartvar`` logger tree.

    ⚠ WITHOUT THIS THE VERDICT IS UNREADABLE, which is the whole failure mode
    this module's docstring warns about — and it very nearly shipped.
    ``app._record_deploy_event`` attaches its ``FileHandler`` to
    ``logging.getLogger("heartvar")``, and ``/api/admin/log-content`` serves ONLY
    those files (categories ``db_builder`` and ``deployments`` under
    ``HEARTVAR_LOGS_DIR``). It does NOT serve container stdout. With no Azure
    access — deploy-doc item A0 is still an open ask to IT — the admin page is
    the ONLY readable surface, so a line that exists only on stdout is a line
    nobody can read.

    asyncio's logger sits OUTSIDE the ``heartvar`` tree, so its records never
    reach that handler. Forwarding is done with a handler rather than by
    re-parenting the logger so it is order-independent: the deployment
    FileHandler is added AFTER the probe is armed, and propagation is resolved
    at emit time, so the record lands in whatever handlers exist by then.

    No recursion risk: ``heartvar``'s handlers write to a file and to stdout,
    neither of which logs back to ``asyncio``.
    """

    def emit(self, record: logging.LogRecord) -> None:
        try:
            logging.getLogger("heartvar.loop_probe.asyncio").handle(record)
        except Exception:  # pragma: no cover - logging must never raise
            self.handleError(record)


def _forward_asyncio_into_heartvar(asyncio_log: logging.Logger) -> None:
    """Install the forwarder once, idempotently (the lifespan can run twice in
    tests, and a duplicate handler would double every line)."""
    if any(isinstance(h, _ForwardToHeartvar) for h in asyncio_log.handlers):
        return
    asyncio_log.addHandler(_ForwardToHeartvar())


def maybe_enable_loop_probe() -> bool:
    """Turn on asyncio's slow-callback reporting IFF ``HEARTVAR_LOOP_PROBE``.

    Returns True when the probe was armed, False when it was not — in which case
    NOTHING is touched: not the loop's debug flag, not ``slow_callback_duration``,
    not the logging config. Call from the FastAPI lifespan, where a loop is
    already running.

    BOTH settings are required. asyncio gates the slow-callback check on
    ``loop.get_debug()``, so setting ``slow_callback_duration`` alone silently
    does nothing at all — the exact shape of a diagnostic that looks armed and
    reports nothing.
    """
    if not probe_enabled():
        return False
    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        log.warning("[loop-probe] HEARTVAR_LOOP_PROBE is set but there is no "
                    "running event loop — probe NOT armed")
        return False

    threshold = _threshold()
    loop.set_debug(True)
    loop.slow_callback_duration = threshold

    asyncio_log = logging.getLogger("asyncio")
    if not asyncio_log.isEnabledFor(logging.WARNING):
        asyncio_log.setLevel(logging.WARNING)
    _forward_asyncio_into_heartvar(asyncio_log)

    log.warning(
        "[loop-probe] ARMED — asyncio debug mode ON, slow_callback_duration="
        "%.2fs. Every callback is now timed, which is REAL OVERHEAD: run one "
        "curation, grep the log for 'took', then unset HEARTVAR_LOOP_PROBE.",
        threshold,
    )
    return True
