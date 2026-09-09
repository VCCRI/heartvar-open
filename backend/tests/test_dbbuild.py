"""Tests for scripts/_dbbuild.py — the shared "build SQLite on local disk,
publish to the (possibly SMB) data mount" helper every DB builder now uses.

The deploy target mounts data/ from Azure Files (SMB), where SQLite's fcntl
locking is unreliable and building a DB in place raised "database is locked".
These lock in the two properties that fix depends on:
  1. staging_db_path() never returns a path inside the destination directory
     (so no SQLite locking ever touches the mount), and
  2. publish() moves a finished DB across directories/filesystems atomically,
     leaving no partial file behind.
"""
from __future__ import annotations

import errno
import importlib.util
import os
import sqlite3
from pathlib import Path

import pytest

_MOD = Path(__file__).resolve().parent.parent.parent / "scripts" / "_dbbuild.py"
_spec = importlib.util.spec_from_file_location("dbbuild_under_test", _MOD)
db = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(db)


def _build_toy_db(path: Path, rows) -> None:
    """Do exactly what a builder does: create a table and insert rows."""
    conn = db.connect(path)
    conn.execute("PRAGMA journal_mode = OFF")
    conn.execute("CREATE TABLE t (k TEXT, v INTEGER)")
    conn.executemany("INSERT INTO t (k, v) VALUES (?, ?)", rows)
    conn.commit()
    conn.close()


def test_staging_path_is_not_on_the_destination_mount(tmp_path):
    mount = tmp_path / "app" / "data"
    mount.mkdir(parents=True)
    final = mount / "uniprot.db"

    staged = db.staging_db_path(final)

    assert mount not in staged.parents
    assert staged.name == "uniprot.db"


def test_staging_path_respects_tmpdir_override(tmp_path, monkeypatch):
    scratch_root = tmp_path / "local-scratch"
    monkeypatch.setenv("HEARTVAR_DBBUILD_TMPDIR", str(scratch_root))

    staged = db.staging_db_path(tmp_path / "data" / "gtex.db")

    assert scratch_root in staged.parents
    assert scratch_root.is_dir()


def test_connect_sets_busy_timeout(tmp_path):
    conn = db.connect(tmp_path / "x.db", busy_timeout_ms=12_345)
    try:
        (got,) = conn.execute("PRAGMA busy_timeout").fetchone()
        assert got == 12_345
    finally:
        conn.close()


def test_publish_uses_copy_then_same_dir_rename_not_cross_fs_replace(tmp_path, monkeypatch):
    real_replace = db.os.replace

    def guarded_replace(src, dst, *a, **k):
        if os.path.dirname(str(src)) != os.path.dirname(str(dst)):
            raise OSError(errno.EXDEV, "simulated cross-device rename")
        return real_replace(src, dst, *a, **k)

    monkeypatch.setattr(db.os, "replace", guarded_replace)

    staged = db.staging_db_path(tmp_path / "data" / "clinvar.db")
    final = tmp_path / "data" / "clinvar.db"
    assert staged.parent != final.parent
    _build_toy_db(staged, [("a", 1), ("b", 2)])

    db.publish(staged, final)

    assert final.exists()
    conn = sqlite3.connect(final)
    try:
        assert conn.execute("SELECT COUNT(*) FROM t").fetchone()[0] == 2
    finally:
        conn.close()
    assert not list(final.parent.glob("*.stage.*"))


def test_publish_reaps_orphan_stage_files_from_a_prior_crash(tmp_path):
    final = tmp_path / "data" / "clinvar.db"
    final.parent.mkdir(parents=True)
    (final.parent / "clinvar.db.stage.99999").write_bytes(b"orphan from a SIGKILL")
    other = final.parent / "uniprot.db.stage.88888"
    other.write_bytes(b"someone else's staging")

    staged = db.staging_db_path(final)
    _build_toy_db(staged, [("x", 1)])
    db.publish(staged, final)

    assert not (final.parent / "clinvar.db.stage.99999").exists()
    assert other.exists()
    assert not list(final.parent.glob("clinvar.db.stage.*"))


def test_publish_overwrites_existing_destination(tmp_path):
    final = tmp_path / "data" / "medgen.db"
    final.parent.mkdir(parents=True)
    final.write_bytes(b"stale-not-a-db")

    staged = db.staging_db_path(final)
    _build_toy_db(staged, [("only", 42)])
    db.publish(staged, final)

    conn = sqlite3.connect(final)
    try:
        assert conn.execute("SELECT v FROM t WHERE k='only'").fetchone()[0] == 42
    finally:
        conn.close()


def test_publish_leaves_no_partial_file_when_copy_fails(tmp_path):
    final = tmp_path / "data" / "mgi.db"
    final.parent.mkdir(parents=True)
    missing = tmp_path / "does-not-exist.db"

    with pytest.raises(OSError):
        db.publish(missing, final)

    assert not final.exists()
    assert not list(final.parent.glob("*.stage.*"))
