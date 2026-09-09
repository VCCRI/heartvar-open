"""Shared pytest fixtures for the backend test suite.

The only fixture here resets the process-wide external-lookup cache before
every test. Several client tests drive the SAME key through different mocked
upstream responses (e.g. test_pmcoa runs the same PMID as available, then
not-open-access, then not-in-PMC); without a reset the second test would be
served the first test's cached result. Clearing before each test keeps them
independent regardless of collection order.

Autouse so individual test files need no per-test boilerplate. Standalone
``python -m backend.tests.test_*`` runs don't load conftest, so the few files
that reuse a key across their own ``__main__`` runner clear the cache inline.

This module also redirects HEARTVAR_LOGS_DIR at import time — see the comment
below; that one is not a fixture because its target is read at import.
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest

os.environ.setdefault(
    "HEARTVAR_LOGS_DIR", tempfile.mkdtemp(prefix="heartvar-test-logs-")
)

from backend.clients._cache import EXTERNAL_CACHE  # noqa: E402


@pytest.fixture(autouse=True)
def _reset_external_cache():
    EXTERNAL_CACHE.clear()
    yield
    EXTERNAL_CACHE.clear()


import pytest as _pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def require_file(relative_path: str, why: str = ""):
    """Skip unless `relative_path` exists under the repo root. Returns the Path."""
    target = REPO_ROOT / relative_path
    if not target.exists():
        _pytest.skip(f"{relative_path} is not present in this checkout"
                     + (f" — {why}" if why else ""))
    return target


def require_db(name: str):
    """Skip unless a built SQLite cache is present. Build with
    ``python3 scripts/build_<name>_db.py``."""
    target = REPO_ROOT / "data" / f"{name}.db"
    if not target.exists():
        _pytest.skip(f"data/{name}.db not built — run scripts/build_{name}_db.py")
    return target


def require_unstripped_prose(text: str, marker: str):
    """Skip when the prose a test asserts on has been stripped.

    The public snapshot removes comments and truncates docstrings, so a test
    that checks the WORDING of a comment is meaningless there. It still runs in
    full against the source tree, which is where a stale comment gets written.
    """
    if marker not in (text or ""):
        _pytest.skip("comment prose is stripped in this tree; "
                     "this assertion only applies to the source repository")
