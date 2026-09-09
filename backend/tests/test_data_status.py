"""Tests for GET /api/data-status — the About page's "Last update" date.

The About page used to carry a hand-typed "Data is updated monthly. Last update:
21/07/2026." That claim was false: the monthly refresh had been skipping every
already-present cache since first build (see monthly_db_update.txt and
backend/tests/test_build_all_cadence.py), so the single freshness signal users
could see was the one nobody was maintaining.

The date is now derived from ``data/build_stamp.json``, written by
scripts/build_all.sh. The properties that make it trustworthy, and which these
tests pin down:

  * a run that rebuilt nothing must NOT advance the date (otherwise a crash-retry
    makes a stale mirror look fresh — the original bug, re-expressed);
  * a missing or corrupt stamp yields no date at all rather than a guess;
  * the endpoint is public (the About page is public) and leaks no build topology.
"""
from __future__ import annotations

import json

import pytest
from starlette.testclient import TestClient

from backend.app import PROJECT_ROOT, app


@pytest.fixture
def client():
    with TestClient(app) as c:
        yield c


def _write_stamp(tmp_path, monkeypatch, payload) -> None:
    p = tmp_path / "build_stamp.json"
    p.write_text(payload if isinstance(payload, str) else json.dumps(payload))
    monkeypatch.setenv("HEARTVAR_BUILD_STAMP_PATH", str(p))


FULL_STAMP = {
    "last_run_utc": "2026-09-01T00:14:02+00:00",
    "last_refresh_utc": "2026-09-01T00:14:02+00:00",
    "mode": "--monthly (refresh monthly-cadence sources)",
    "refreshed": ["clinvar", "uniprot"],
    "built": ["clinvar", "uniprot", "opentargets"],
    "kept": ["gnomad_freq", "gtex"],
    "failed": [],
}


def test_reports_the_last_refresh(client, tmp_path, monkeypatch):
    _write_stamp(tmp_path, monkeypatch, FULL_STAMP)
    body = client.get("/api/data-status").json()
    assert body["last_refresh_utc"] == "2026-09-01T00:14:02+00:00"
    assert body["refreshed_count"] == 2
    assert body["built_count"] == 3
    assert body["failed_count"] == 0


def test_a_noop_run_does_not_advance_the_date(client, tmp_path, monkeypatch):
    """A crash-retry that rebuilt nothing keeps the older refresh date.

    This is the whole point: --resume skips everything present, and if that
    advanced the displayed date, the page would report a fresh mirror on the exact
    failure that was hiding staleness in the first place.
    """
    _write_stamp(tmp_path, monkeypatch, {
        "last_run_utc": "2026-09-15T03:00:00+00:00",
        "last_refresh_utc": "2026-09-01T00:14:02+00:00",
        "refreshed": [], "built": [], "kept": ["clinvar"], "failed": [],
    })
    body = client.get("/api/data-status").json()
    assert body["last_refresh_utc"] == "2026-09-01T00:14:02+00:00"
    assert body["last_run_utc"] == "2026-09-15T03:00:00+00:00"
    assert body["built_count"] == 0


def test_missing_stamp_falls_back_to_file_mtime(client, tmp_path, monkeypatch):
    """No stamp, but data files present → report when they were last written.

    This is the state of a mount built before build_all.sh wrote stamps, and of a
    build still in progress (the stamp lands at the end). A real date from the
    files beats a blank, and it upgrades itself to the exact stamp on the next
    completed build.
    """
    import os
    import time

    monkeypatch.setenv("HEARTVAR_BUILD_STAMP_PATH", str(tmp_path / "absent.json"))
    monkeypatch.setenv("HEARTVAR_DATA_DIR", str(tmp_path))
    db = tmp_path / "clinvar.db"
    db.write_text("x")
    when = time.mktime((2026, 7, 21, 4, 30, 0, 0, 0, -1))
    os.utime(db, (when, when))

    body = client.get("/api/data-status").json()
    assert body["source"] == "file_mtime"
    assert body["last_refresh_utc"].startswith("2026-07-2")


def test_stamp_wins_over_mtime(client, tmp_path, monkeypatch):
    """The stamp is authoritative: it knows whether anything was actually
    rebuilt, which an mtime cannot."""
    _write_stamp(tmp_path, monkeypatch, FULL_STAMP)
    monkeypatch.setenv("HEARTVAR_DATA_DIR", str(tmp_path))
    (tmp_path / "clinvar.db").write_text("x")
    body = client.get("/api/data-status").json()
    assert body["source"] == "build_stamp"
    assert body["last_refresh_utc"] == "2026-09-01T00:14:02+00:00"


def test_versioned_caches_do_not_count_as_a_refresh(client, tmp_path, monkeypatch):
    """gnomAD/GTEx/AlphaFold are provisioned once and never refreshed, so their
    mtime is a first-provision date. Counting it would report the mirror as
    refreshed when it never was."""
    monkeypatch.setenv("HEARTVAR_BUILD_STAMP_PATH", str(tmp_path / "absent.json"))
    monkeypatch.setenv("HEARTVAR_DATA_DIR", str(tmp_path))
    for name in ("gnomad_freq.db", "gtex.db", "gnomad_constraint.db"):
        (tmp_path / name).write_text("x")
    assert client.get("/api/data-status").json() == {}


def test_no_stamp_and_no_data_yields_no_date(client, tmp_path, monkeypatch):
    """Empty mount → {} → the frontend shows "Data is updated monthly." undated.
    Never a fabricated date."""
    monkeypatch.setenv("HEARTVAR_BUILD_STAMP_PATH", str(tmp_path / "absent.json"))
    monkeypatch.setenv("HEARTVAR_DATA_DIR", str(tmp_path))
    assert client.get("/api/data-status").json() == {}


def test_corrupt_stamp_yields_no_date(client, tmp_path, monkeypatch):
    """A half-written file must not 500 a public page."""
    _write_stamp(tmp_path, monkeypatch, '{"last_refresh_utc": "2026-09-01T0')
    assert client.get("/api/data-status").json() == {}


def test_non_object_stamp_yields_no_date(client, tmp_path, monkeypatch):
    _write_stamp(tmp_path, monkeypatch, "[1, 2, 3]")
    assert client.get("/api/data-status").json() == {}


def test_endpoint_is_public_and_leaks_no_source_lists(client, tmp_path, monkeypatch):
    """The About page is public, so this must need no session. It reports counts
    only — which sources exist and which failed is operator detail."""
    _write_stamp(tmp_path, monkeypatch, FULL_STAMP)
    resp = client.get("/api/data-status")
    assert resp.status_code == 200
    body = resp.json()
    assert "refreshed" not in body and "built" not in body and "kept" not in body
    assert set(body) == {
        "last_refresh_utc", "last_run_utc", "source",
        "refreshed_count", "built_count", "failed_count",
    }


def test_about_page_ships_no_hardcoded_date():
    """Regression guard: the date must come from the build, not the markup.

    If someone re-hardcodes it, the page can silently disagree with reality again
    — which is exactly what happened with 21/07/2026.
    """
    html = (PROJECT_ROOT / "index.html").read_text(encoding="utf-8")
    assert "src-update-note" in html, "the About freshness note has moved or gone"
    note_line = next(
        line for line in html.splitlines() if "src-update-note" in line
    )
    assert "Last update:" not in note_line, (
        "the About note hardcodes a date again — it must be injected from "
        "/api/data-status by static/heartvar.datastamp.js"
    )


def test_vep_is_a_reported_source():
    from backend.app import _DB_SOURCE_ARTIFACTS, _DB_SOURCE_CADENCE
    assert _DB_SOURCE_ARTIFACTS["vep"] == "data/vep/.heartvar_vep_manifest.json"
    assert _DB_SOURCE_CADENCE["vep"] == "versioned"


def test_vep_artifact_matches_build_all_output_for():
    """This registry mirrors output_for() in build_all.sh. Drift means the portal
    reports on a file the builder never writes — which is how a source silently
    reads as "not built" while the app is happily using it."""
    import re
    from pathlib import Path as _Path
    from backend.app import _DB_SOURCE_ARTIFACTS

    script = (_Path(__file__).resolve().parents[2] / "scripts" / "build_all.sh").read_text()
    match = re.search(r'^\s*vep\)\s*echo "\$DATA_DIR/([^"]+)";;', script, re.MULTILINE)
    assert match, "output_for() has no vep arm"
    assert _DB_SOURCE_ARTIFACTS["vep"] == f"data/{match.group(1)}"


def test_vep_is_not_in_the_env_override_map():
    """HEARTVAR_VEP_DATA is a DIRECTORY; that map expects the file the reader
    opens. Listing it would make the portal report a bogus path — the exact bug
    the map's own comment warns about."""
    from backend.app import _DB_SOURCE_ENV_OVERRIDES
    assert "vep" not in _DB_SOURCE_ENV_OVERRIDES


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
