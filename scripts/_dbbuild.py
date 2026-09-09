"""_dbbuild.py — build SQLite caches safely on any deploy filesystem.

The deploy target mounts ``data/`` from an Azure Files (SMB) share. SQLite's
locking relies on POSIX ``fcntl`` byte-range locks, which SMB/CIFS/NFS network
mounts do not reliably honour — so building a database *directly on the mount*
can raise ``sqlite3.OperationalError: database is locked`` (or a "disk I/O
error") even with a single writer, and definitely when two builder replicas
overlap.

The fix every builder shares:

  1. Build the DB on **local disk** (an ephemeral scratch dir — never the
     mount), where SQLite locking works normally.
  2. When it's finished, ``publish()`` it onto the final path with a plain
     file copy + a same-directory atomic rename. The share only ever sees a
     sequential write and a metadata rename, both of which it handles fine —
     no SQLite locking ever touches the network mount.

Usage in a builder::

    import _dbbuild
    ...
    tmp_path = _dbbuild.staging_db_path(db_path)   # local scratch, not the mount
    conn = _dbbuild.connect(tmp_path)              # adds busy_timeout
    conn.execute("PRAGMA journal_mode = OFF")      # existing bulk pragmas stay
    ... build ...
    conn.close()
    _dbbuild.publish(tmp_path, db_path)            # copy to mount + atomic rename

The scratch directory defaults to the OS temp dir (``$TMPDIR`` / ``/tmp``,
which is local container disk). Override with ``HEARTVAR_DBBUILD_TMPDIR`` if
that filesystem is too small for a large DB — but it MUST point at local disk,
never the ``data/`` mount, or the SMB locking problem comes straight back.
"""

from __future__ import annotations

import atexit
import os
import shutil
import time
import sqlite3
import tempfile
from pathlib import Path

_BUSY_TIMEOUT_MS = 60_000

_SCRATCH_DIRS: list[str] = []


def staging_db_path(final_path: os.PathLike | str) -> Path:
    """Return a path on LOCAL disk to build a SQLite DB before publishing it.

    The returned path lives in a fresh scratch directory (cleaned up at process
    exit), NOT next to ``final_path`` — which may be a network/SMB mount where
    SQLite locking is unreliable. Keeps the DB's basename so log lines and the
    published file match.
    """
    final_path = Path(final_path)
    parent = os.environ.get("HEARTVAR_DBBUILD_TMPDIR") or None
    if parent:
        Path(parent).mkdir(parents=True, exist_ok=True)
    scratch = tempfile.mkdtemp(prefix="heartvar-dbbuild-", dir=parent)
    _SCRATCH_DIRS.append(scratch)
    return Path(scratch) / final_path.name


def connect(db_path: os.PathLike | str, *, busy_timeout_ms: int = _BUSY_TIMEOUT_MS) -> sqlite3.Connection:
    """``sqlite3.connect`` with a busy_timeout so a transient lock waits rather
    than instantly raising ``OperationalError('database is locked')``.

    Callers keep setting their own bulk-load pragmas (journal_mode/synchronous/
    temp_store) afterwards exactly as before — this only adds the timeout.
    """
    conn = sqlite3.connect(db_path, timeout=busy_timeout_ms / 1000.0)
    conn.execute(f"PRAGMA busy_timeout = {int(busy_timeout_ms)}")
    return conn


def _replace_with_retry(staged: Path, final_path: Path,
                        attempts: int = 10, delay: float = 3.0) -> None:
    """``os.replace`` the staged file into place, retrying on a denied rename.

    ⚠ AN OPEN READER CAN DENY THIS, and on Azure Files SMB it does. Observed
    2026-08-31 on a ClinVar build:
        PermissionError: [Errno 13] Permission denied:
          '/app/data/clinvar.db.stage.33' -> '/app/data/clinvar.db'
    Renaming over an existing file needs delete access to the target, and SMB
    refuses while another process holds it open. The app is exactly such a
    process: it serves curations out of these databases continuously.

    The permanent case was warming holding a connection for the life of the
    process (fixed in backend/warmup.py), but the transient case is unfixable
    and legitimate — a curation may simply be mid-query when the build lands.
    So this retries rather than failing a build that did all its work correctly:
    30 s of attempts covers any single query, and a build is a monthly event
    where 30 s costs nothing.

    On final failure the error names the cause, because "Permission denied"
    reads like a filesystem ACL problem and is not one.
    """
    last: OSError | None = None
    for attempt in range(attempts):
        try:
            os.replace(staged, final_path)
            if attempt:
                print(f"    published on attempt {attempt + 1} "
                      f"(the target was held open)")
            return
        except PermissionError as exc:
            last = exc
            if attempt < attempts - 1:
                print(f"    rename denied (attempt {attempt + 1}/{attempts}) — "
                      f"a reader holds {final_path.name} open; retrying in "
                      f"{delay:.0f}s")
                time.sleep(delay)
    raise RuntimeError(
        f"could not publish {final_path.name}: the rename was denied for "
        f"{attempts * delay:.0f}s. On Azure Files SMB this means a process "
        f"holds the target open, NOT a filesystem permission problem. The "
        f"database was built correctly and is at {staged}. Last error: {last}"
    )

def publish(local_path: os.PathLike | str, final_path: os.PathLike | str) -> None:
    """Move a finished DB from local disk onto ``final_path`` atomically.

    ``os.replace`` cannot cross filesystems (local scratch -> network mount
    raises ``EXDEV``), so copy to a pid-unique sibling *on the destination
    filesystem* and then same-directory rename it into place. That rename is
    atomic and is a plain metadata op the SMB share handles fine; the pid
    suffix keeps two overlapping builds from clobbering each other's staging
    file. Readers never observe a partially-copied DB.

    A prior run hard-killed mid-copy (SIGKILL / OOM — the ``finally`` never
    runs) can leave a ``<db>.stage.<pid>`` orphan on the mount, so we reap any
    stale staging siblings of *this* DB before copying. That's safe because
    build_all.sh serialises builders under a single-writer lock, so no live
    peer owns one — restoring the self-healing the old fixed-name temp had.
    """
    local_path = Path(local_path)
    final_path = Path(final_path)
    final_path.parent.mkdir(parents=True, exist_ok=True)
    for orphan in final_path.parent.glob(f"{final_path.name}.stage.*"):
        try:
            orphan.unlink()
        except OSError:
            pass
    staged = final_path.with_name(f"{final_path.name}.stage.{os.getpid()}")
    try:
        shutil.copyfile(local_path, staged)
        _replace_with_retry(staged, final_path)
    finally:
        if staged.exists():
            try:
                staged.unlink()
            except OSError:
                pass


@atexit.register
def _cleanup_scratch() -> None:
    for d in _SCRATCH_DIRS:
        shutil.rmtree(d, ignore_errors=True)
