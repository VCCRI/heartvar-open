"""Unit tests for the shared local-DB path helper (codebase audit P3.2 +
deployment-risk audit Part 1A): one source of truth for the repo-root data/
directory, plus an env override for deploy flexibility.

No pytest dependency — runnable directly
(``python -m backend.tests.test_client_paths``).
"""
from __future__ import annotations

import os

from backend.clients import _paths
from backend.clients import clinvar, biogrid, fetal_heart


def test_data_path_under_repo_root_data_dir():
    p = _paths.data_path("clinvar.db")
    assert p == _paths.PROJECT_ROOT / "data" / "clinvar.db"
    assert p.parent.name == "data"


def test_db_path_default_falls_back_to_data_dir():
    os.environ.pop("HEARTVAR_TEST_DB_PATH", None)
    assert _paths.db_path("x.db", "HEARTVAR_TEST_DB_PATH") == _paths.DATA_DIR / "x.db"


def test_db_path_honours_env_override():
    os.environ["HEARTVAR_TEST_DB_PATH"] = "/mnt/elsewhere/clinvar.db"
    try:
        assert str(_paths.db_path("x.db", "HEARTVAR_TEST_DB_PATH")) == "/mnt/elsewhere/clinvar.db"
    finally:
        os.environ.pop("HEARTVAR_TEST_DB_PATH", None)


def test_clients_resolve_to_repo_root_data_dir():
    assert clinvar.DB_PATH == _paths.DATA_DIR / "clinvar.db"
    assert biogrid.DB_PATH == _paths.DATA_DIR / "biogrid.db"
    assert fetal_heart.DB_PATH == _paths.DATA_DIR / "fetal_heart.db"
    assert clinvar.DB_PATH.parent == biogrid.DB_PATH.parent == fetal_heart.DB_PATH.parent


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
