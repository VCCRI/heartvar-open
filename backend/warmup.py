"""Warm the local data mount at start-up, so no curator pays first access.

MEASURED, 2026-08-30, on a cold production container. The RYR2 curation took
42.6 s end to end: ~20 s of DB gather and ~15.9 s of ``same_site``. Run again on
the now-WARM container, the same curation's ``vep`` source fell from 18.58 s to
273 ms. So a large share of the cost is FIRST ACCESS to files on the SMB mount —
and because every deploy restarts the container, that bill lands on whichever
curator happens to click first.

The lifespan already does precisely this for one cache ("Fire the JAX HPO
descendant cache warm-up in the background so the first real curation doesn't pay
the ~9 s cold-start cost"). This generalises it to the mount.

WHY IT ALSO LOGS TIMINGS. Every warmed file gets a per-file duration in the
deployment log, which /api/admin/log-content serves. So the next deploy states
plainly which file on the mount is slow, readable without Azure access and
without arming the loop probe. The diagnostic is a side effect of the fix.

⚠ BEST EFFORT, ALWAYS. Warming runs in the background, never blocks readiness,
and swallows every per-file failure. A corrupt database or a mount that is not
attached yet must not stop the server from starting — it only means the first
curation pays what it used to.
"""
from __future__ import annotations

import asyncio
import logging
import os
import sqlite3
import time
from pathlib import Path

from . import localio
from .clients._paths import DATA_DIR

log = logging.getLogger("heartvar.warmup")

_EXTRA_ENV_PATHS = (
    "HEARTVAR_CDOT_DB",
    "GENCC_SNAPSHOT_PATH",
    "HGNC_ALIAS_MAP_PATH",
    "HEARTVAR_EREPO_TSV",
)


def warming_enabled() -> bool:
    """ON by default — it is the fix, not an experiment. ``0``/``false``/``no``/
    ``off`` disables it."""
    return os.environ.get("HEARTVAR_WARM_LOCAL_DATA", "").strip().lower() \
        not in {"0", "false", "no", "off"}


_DEFAULT_WARM_CONCURRENCY = 2


def warm_concurrency() -> int:
    """How many files warming may read AT ONCE. Deliberately small.

    ⚠ MEASURED, and it is the reason this knob exists. Unbounded, warming pulled
    every database and the 847 MB FASTA off the mount simultaneously — so a
    curation that arrived DURING start-up contended with it and took 31.7 s,
    against 5.4 s once warming had finished. Warming exists to protect the first
    curation; saturating the mount to do it defeats the purpose.

    2 is slow enough to leave the mount real headroom and still finishes start-up
    warming in well under the time a first curation used to cost. Raise it only
    if the deployment log shows warming still unfinished when traffic arrives.
    """
    raw = os.environ.get("HEARTVAR_WARM_CONCURRENCY", "").strip()
    if raw:
        try:
            value = int(raw)
            if value > 0:
                return value
            log.warning("HEARTVAR_WARM_CONCURRENCY=%r must be > 0 — using %d",
                        raw, _DEFAULT_WARM_CONCURRENCY)
        except ValueError:
            log.warning("HEARTVAR_WARM_CONCURRENCY=%r is not an int — using %d",
                        raw, _DEFAULT_WARM_CONCURRENCY)
    return _DEFAULT_WARM_CONCURRENCY


def _data_dir() -> Path:
    """The mount. ``/app/data`` in production, which is where every default
    client path resolves to as well."""
    return Path(os.environ.get("HEARTVAR_DATA_DIR") or DATA_DIR)


def _warm_one_db(path: Path) -> None:
    """Touch a SQLite file so its header, schema and root pages are cached.

    ``SELECT 1`` alone would not read the schema, and the schema read is the part
    that costs several SMB round trips, so the table list is queried instead.
    """
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        con.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    finally:
        con.close()


def _warm_one_file(path: Path) -> None:
    """Read the head of a plain file, enough to pull it into the page cache."""
    with open(path, "rb") as fh:
        fh.read(1024 * 1024)


async def warm_local_data(data_dir: str | os.PathLike[str] | None = None) -> dict:
    """Touch every local database (and the big known artefacts), concurrently.

    Returns ``{"warmed": [...], "failed": [...], "seconds": {name: secs}}``.
    Never raises: each item is guarded, and the caller is a background task.
    """
    root = Path(data_dir) if data_dir is not None else _data_dir()
    warmed: list[str] = []
    failed: list[str] = []
    seconds: dict[str, float] = {}

    targets: list[tuple[str, Path, bool]] = []
    try:
        for p in sorted(root.glob("*.db")):
            targets.append((p.name, p, True))
    except OSError as exc:
        log.warning("[warmup] cannot list %s (%r)", root, exc)

    if data_dir is None:
        for env in _EXTRA_ENV_PATHS:
            raw = os.environ.get(env, "").strip()
            if not raw:
                continue
            p = Path(raw)
            targets.append((p.name, p, p.suffix == ".db"))

    gate = asyncio.Semaphore(warm_concurrency())

    async def one(label: str, path: Path, is_db: bool) -> None:
        async with gate:
            await _warm_guarded(label, path, is_db)

    async def _warm_guarded(label: str, path: Path, is_db: bool) -> None:
        t0 = time.perf_counter()
        try:
            await localio.run_local(_warm_one_db if is_db else _warm_one_file, path)
        except Exception as exc:  # noqa: BLE001 — best effort, by design
            failed.append(label)
            seconds[label] = time.perf_counter() - t0
            log.warning("[warmup] %s FAILED after %.2fs (%r)",
                        label, seconds[label], exc)
            return
        seconds[label] = time.perf_counter() - t0
        warmed.append(label)

    if targets:
        await asyncio.gather(*(one(*t) for t in targets))
        slow = sorted(seconds.items(), key=lambda kv: -kv[1])
        log.info("[warmup] warmed %d/%d local files (concurrency %d) — "
                 "total file-time %.2fs, slowest: %s",
                 len(warmed), len(targets), warm_concurrency(),
                 sum(seconds.values()),
                 ", ".join(f"{n}={s:.2f}s" for n, s in slow[:8]) or "none")
    else:
        log.info("[warmup] no local data files found under %s", root)

    return {"warmed": warmed, "failed": failed, "seconds": seconds}


async def warm_resolver() -> None:
    """Build the HGVS resolver's providers and open the FASTA up front.

    THE SINGLE BIGGEST COLD COST MEASURED: ~18.5 s, and it lands on ``fetch_vep``
    because ``hgvs_resolver.resolve()`` is called from ``async fetch_hgvs``. The
    cdot provider build itself is fast (0.36 s in production), so the remainder is
    ``_fasta()`` opening the 847 MB bgzipped genome and its .fai/.gzi on the
    mount. Both are ``lru_cache``d, so doing it once here spends it where nobody
    is waiting.
    """
    from .clients import hgvs_resolver
    if not hgvs_resolver.resolver_enabled():
        return
    t0 = time.perf_counter()
    try:
        await localio.run_local(hgvs_resolver._providers)
        await localio.run_local(hgvs_resolver._fasta)
    except Exception as exc:  # noqa: BLE001 — best effort
        log.warning("[warmup] resolver warm FAILED after %.2fs (%r)",
                    time.perf_counter() - t0, exc)
        return
    log.info("[warmup] resolver + FASTA ready in %.2fs (this is the cost a "
             "first curation used to pay)", time.perf_counter() - t0)
