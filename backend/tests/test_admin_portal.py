"""Tests for the admin portal — /api/admin/log-index, log-content, db-status.

The portal shows build and deployment logs and per-source data status. Three
things make it worth pinning:

  * **It is the only authenticated-and-privileged surface in the app.** Every
    other endpoint is either public (evidence, the rule-based classification) or
    gated on merely being signed in (the AI paths). These three answer only to
    an allow-listed address, and the allow-list is an env var — so "unset" must
    mean closed, not open. That is the failure the auth work already learned
    once: ``require_auth`` existed but was attached to nothing.
  * **It reads files chosen by the client.** ``category`` and ``name`` come off
    the query string. Neither is ever joined onto a path: the category selects a
    constant directory name from a literal map, and the filename is matched
    against that directory's listing, so only a file the server itself listed
    can be opened.
  * **The logs contain curation queries.** The deployment handler captures every
    heartvar INFO record, and the app logs gene + HGVS. So a leak here is a leak
    of what was curated, not just of server internals.

No network and no filesystem writes outside tmp_path: _LOG_BASE_DIR is
redirected per test.
"""
from __future__ import annotations

import pathlib

import pytest
from starlette.testclient import TestClient

import backend.app as app_module
from backend import ratelimit

ADMIN = {
    "oid": "11111111-1111-1111-1111-111111111111",
    "preferred_username": "Admin.Person@victorchang.edu.au",
    "name": "Admin Person",
}
CURATOR = {
    "oid": "22222222-2222-2222-2222-222222222222",
    "preferred_username": "curator@hospital.org.au",
    "name": "Ordinary Curator",
}

ADMIN_ENDPOINTS = [
    "/api/admin/log-index",
    "/api/admin/log-content?category=deployments&name=x.log",
    "/api/admin/db-status",
]


@pytest.fixture
def client():
    return TestClient(app_module.app)


@pytest.fixture(autouse=True)
def _no_inbound_rate_limit(monkeypatch):
    monkeypatch.setattr(ratelimit.limiter, "enabled", False)
    yield


@pytest.fixture
def logs_dir(tmp_path, monkeypatch):
    """Point the log store at tmp_path and return it."""
    base = tmp_path / "logs"
    (base / "db_builder").mkdir(parents=True)
    (base / "deployments").mkdir(parents=True)
    monkeypatch.setattr(app_module, "_LOG_BASE_DIR", base)
    return base


def signed_in_as(monkeypatch, account: dict | None):
    monkeypatch.setattr(app_module, "session_account", lambda request: account)


def admin_list(monkeypatch, value: str):
    monkeypatch.setenv("HEARTVAR_ADMIN_EMAILS", value)


@pytest.mark.parametrize("url", ADMIN_ENDPOINTS)
def test_anonymous_is_refused(client, monkeypatch, logs_dir, url):
    signed_in_as(monkeypatch, None)
    admin_list(monkeypatch, "admin.person@victorchang.edu.au")
    assert client.get(url).status_code == 403


@pytest.mark.parametrize("url", ADMIN_ENDPOINTS)
def test_an_ordinary_signed_in_curator_is_refused(client, monkeypatch, logs_dir, url):
    """Being signed in buys the AI paths, not the portal."""
    signed_in_as(monkeypatch, dict(CURATOR))
    admin_list(monkeypatch, "admin.person@victorchang.edu.au")
    assert client.get(url).status_code == 403


@pytest.mark.parametrize("url", ADMIN_ENDPOINTS)
def test_unset_allow_list_closes_the_portal(client, monkeypatch, logs_dir, url):
    """An unconfigured allow-list must deny everyone — including an account that
    would be an admin on another deployment. Fail-closed, like ai_auth_required."""
    signed_in_as(monkeypatch, dict(ADMIN))
    monkeypatch.delenv("HEARTVAR_ADMIN_EMAILS", raising=False)
    assert client.get(url).status_code == 403


def test_an_empty_allow_list_closes_the_portal(client, monkeypatch, logs_dir):
    signed_in_as(monkeypatch, dict(ADMIN))
    admin_list(monkeypatch, "   ,  , ")
    assert client.get("/api/admin/log-index").status_code == 403


def test_the_allow_list_is_case_insensitive(client, monkeypatch, logs_dir):
    """Token claims vary in case between providers; the address is the same."""
    signed_in_as(monkeypatch, dict(ADMIN))
    admin_list(monkeypatch, "ADMIN.PERSON@VICTORCHANG.EDU.AU")
    assert client.get("/api/admin/log-index").status_code == 200


def test_the_email_claim_is_honoured_too(client, monkeypatch, logs_dir):
    """Google tokens carry `email` where Entra carries `preferred_username`."""
    signed_in_as(monkeypatch, {"oid": "x", "email": "admin.person@victorchang.edu.au"})
    admin_list(monkeypatch, "admin.person@victorchang.edu.au")
    assert client.get("/api/admin/log-index").status_code == 200


def test_a_substring_of_an_admin_address_is_not_an_admin(client, monkeypatch, logs_dir):
    signed_in_as(monkeypatch, {"oid": "x", "email": "person@victorchang.edu.au"})
    admin_list(monkeypatch, "admin.person@victorchang.edu.au")
    assert client.get("/api/admin/log-index").status_code == 403


def test_auth_status_reports_admin(client, monkeypatch):
    """The nav renders the Admin link off this flag, so it must track the
    server's own answer rather than being assumed by the browser."""
    signed_in_as(monkeypatch, dict(ADMIN))
    admin_list(monkeypatch, "admin.person@victorchang.edu.au")
    assert client.get("/api/auth/status").json()["is_admin"] is True
    signed_in_as(monkeypatch, dict(CURATOR))
    assert client.get("/api/auth/status").json()["is_admin"] is False


@pytest.fixture
def as_admin(monkeypatch):
    signed_in_as(monkeypatch, dict(ADMIN))
    admin_list(monkeypatch, "admin.person@victorchang.edu.au")
    yield


def write_log(directory: pathlib.Path, name: str, body: str = "hello", mtime: int | None = None):
    p = directory / name
    p.write_text(body, encoding="utf-8")
    if mtime is not None:
        import os
        os.utime(p, (mtime, mtime))
    return p


def test_index_lists_newest_first_and_caps_at_the_retention_count(
    client, logs_dir, as_admin,
):
    d = logs_dir / "deployments"
    for i in range(app_module._LOG_MAX_FILES + 3):
        write_log(d, f"run-{i}.log", mtime=1_700_000_000 + i * 60)
    body = client.get("/api/admin/log-index").json()
    names = [e["name"] for e in body["deployments"]]
    assert len(names) == app_module._LOG_MAX_FILES
    assert names[0] == f"run-{app_module._LOG_MAX_FILES + 2}.log"
    assert "run-0.log" not in names
    assert body["db_builder"] == []


def test_index_covers_only_the_allowed_categories(client, logs_dir, as_admin):
    (logs_dir / "secrets").mkdir()
    write_log(logs_dir / "secrets", "keys.log")
    assert set(client.get("/api/admin/log-index").json()) == {"db_builder", "deployments"}


def test_index_is_empty_not_an_error_when_no_logs_exist(client, logs_dir, as_admin):
    body = client.get("/api/admin/log-index").json()
    assert body == {"db_builder": [], "deployments": []}


def test_content_returns_the_file(client, logs_dir, as_admin):
    write_log(logs_dir / "db_builder", "build.log", "line one\nline two\n")
    body = client.get(
        "/api/admin/log-content", params={"category": "db_builder", "name": "build.log"},
    ).json()
    assert body["content"] == "line one\nline two\n"
    assert body["truncated"] is False


@pytest.mark.parametrize("name", [
    "../../../etc/passwd",
    "../deployments/run.log",
    "..%2f..%2fetc%2fpasswd",
    "build.log/../../secret.log",
    "/etc/passwd",
    "build.txt",
    "build.log.bak",
    ".log",
    "",
])
def test_content_rejects_traversal_and_odd_names(client, logs_dir, as_admin, name):
    write_log(logs_dir / "db_builder", "build.log")
    r = client.get(
        "/api/admin/log-content", params={"category": "db_builder", "name": name},
    )
    assert r.status_code in (400, 404, 422), f"{name!r} returned {r.status_code}"


def test_content_rejects_an_unknown_category(client, logs_dir, as_admin):
    write_log(logs_dir / "db_builder", "build.log")
    r = client.get(
        "/api/admin/log-content", params={"category": "secrets", "name": "build.log"},
    )
    assert r.status_code == 400


def test_content_cannot_escape_via_a_symlink(client, logs_dir, as_admin, tmp_path):
    """The filename regex alone would pass a symlink; the resolve()/relative_to()
    check is what stops it pointing outside the category directory."""
    outside = tmp_path / "outside.log"
    outside.write_text("SECRET", encoding="utf-8")
    link = logs_dir / "db_builder" / "sneaky.log"
    try:
        link.symlink_to(outside)
    except (OSError, NotImplementedError):
        pytest.skip("symlinks unavailable on this platform")
    r = client.get(
        "/api/admin/log-content", params={"category": "db_builder", "name": "sneaky.log"},
    )
    assert r.status_code in (400, 404)
    assert "SECRET" not in r.text


def test_content_404s_for_a_missing_file(client, logs_dir, as_admin):
    r = client.get(
        "/api/admin/log-content", params={"category": "deployments", "name": "nope.log"},
    )
    assert r.status_code == 404


def test_content_tails_a_large_file(client, logs_dir, as_admin):
    """A deployment log grows for the life of the process — the response must be
    bounded, and the END is the part an admin wants."""
    cap = app_module._LOG_MAX_READ_BYTES
    body = "".join(f"line {i:07d}\n" for i in range(cap // 8))
    assert len(body) > cap
    write_log(logs_dir / "deployments", "big.log", body)
    data = client.get(
        "/api/admin/log-content", params={"category": "deployments", "name": "big.log"},
    ).json()
    assert data["truncated"] is True
    assert len(data["content"]) <= cap + 200
    assert data["content"].startswith("... (truncated")
    assert body.rstrip("\n").rsplit("\n", 1)[-1] in data["content"]


def test_content_does_not_leak_the_filesystem_path_on_error(
    client, logs_dir, as_admin, monkeypatch,
):
    """CWE-209: the OS error text carries the absolute path of the data mount."""
    write_log(logs_dir / "db_builder", "build.log")

    def _boom(*a, **k):
        raise OSError("EACCES: /mnt/heartvar-data/logs/db_builder/build.log")

    monkeypatch.setattr(pathlib.Path, "open", _boom)
    r = client.get(
        "/api/admin/log-content", params={"category": "db_builder", "name": "build.log"},
    )
    assert r.status_code == 500
    assert "/mnt/" not in r.text and "EACCES" not in r.text


def test_db_status_reports_presence_and_cadence(client, logs_dir, as_admin, tmp_path):
    body = client.get("/api/admin/db-status").json()
    names = {s["name"] for s in body["sources"]}
    assert "clinvar" in names and "alphafold" in names
    by_name = {s["name"]: s for s in body["sources"]}
    assert by_name["clinvar"]["cadence"] == "monthly"
    assert by_name["alphafold"]["cadence"] == "static"
    assert by_name["gnomad_freq"]["cadence"] == "versioned"
    for s in body["sources"]:
        assert set(s) >= {
            "name", "last_updated", "present", "artifact_path", "cadence",
            "failed_in_last_run", "error_snippet",
        }


def test_db_status_survives_a_corrupt_build_stamp(
    client, logs_dir, as_admin, tmp_path, monkeypatch,
):
    stamp = tmp_path / "build_stamp.json"
    stamp.write_text("{not json", encoding="utf-8")
    monkeypatch.setenv("HEARTVAR_BUILD_STAMP_PATH", str(stamp))
    body = client.get("/api/admin/db-status").json()
    assert body["last_run"] is None
    assert body["sources"]


@pytest.fixture
def fake_roots(tmp_path, monkeypatch):
    """Redirect PROJECT_ROOT at tmp_path so presence can be tested by creating
    and withholding artifacts, without touching the real mirror."""
    (tmp_path / "data").mkdir()
    (tmp_path / "backend" / "data").mkdir(parents=True)
    monkeypatch.setattr(app_module, "PROJECT_ROOT", tmp_path)
    for var in (
        "HGNC_ALIAS_MAP_PATH", "HEARTVAR_EREPO_TSV", "GENCC_SNAPSHOT_PATH",
        "CLINGEN_GV_PATH", "CLINVAR_DB_PATH", "HEARTVAR_BUILD_STAMP_PATH",
    ):
        monkeypatch.delenv(var, raising=False)
    return tmp_path


def source_row(body: dict, name: str) -> dict:
    return {s["name"]: s for s in body["sources"]}[name]


def test_db_status_honours_the_env_override_the_reader_uses(
    client, logs_dir, as_admin, fake_roots, monkeypatch,
):
    """HGNC_ALIAS_MAP_PATH is what hgnc_alias.py opens, so it is what decides
    whether the map is really there."""
    mount_copy = fake_roots / "data" / "hgnc_alias_map.json"
    mount_copy.write_text("{}", encoding="utf-8")
    monkeypatch.setenv("HGNC_ALIAS_MAP_PATH", str(mount_copy))

    row = source_row(client.get("/api/admin/db-status").json(), "hgnc_alias")
    assert row["present"] is True
    assert row["last_updated"] is not None


def test_db_status_finds_a_mount_copy_with_no_env_override(
    client, logs_dir, as_admin, fake_roots,
):
    """build_all.sh copies these onto the mount unconditionally, so the copy is a
    legitimate answer even on a deploy that forgot the override — which is the
    failure that hid gencc."""
    (fake_roots / "data" / "gencc_submissions.json").write_text("{}", encoding="utf-8")

    assert source_row(client.get("/api/admin/db-status").json(), "gencc")["present"] is True


def test_db_status_reports_the_path_it_actually_found(
    client, logs_dir, as_admin, fake_roots,
):
    """An admin chasing a missing artifact needs the path that was checked, not
    the one the registry happens to declare."""
    (fake_roots / "data" / "erepo_all.tsv").write_text("x", encoding="utf-8")

    row = source_row(client.get("/api/admin/db-status").json(), "erepo")
    assert row["artifact_path"] == "data/erepo_all.tsv"


def test_db_status_still_reports_a_genuinely_absent_artifact(
    client, logs_dir, as_admin, fake_roots,
):
    """The fix widens where the portal looks; it must not make presence
    unconditional."""
    body = client.get("/api/admin/db-status").json()
    for name in ("gencc", "hgnc_alias", "erepo"):
        assert source_row(body, name)["present"] is False
        assert source_row(body, name)["last_updated"] is None


def test_db_status_finds_a_mounted_db_at_its_declared_path(
    client, logs_dir, as_admin, fake_roots,
):
    """Regression guard: the large DBs already resolved correctly, because
    PROJECT_ROOT/data IS the mount."""
    (fake_roots / "data" / "clinvar.db").write_bytes(b"x")

    assert source_row(client.get("/api/admin/db-status").json(), "clinvar")["present"] is True


def test_db_status_finds_an_in_image_artifact_that_ships_in_the_repo(
    client, logs_dir, as_admin, fake_roots,
):
    """hpo_labels is tracked in git and read from a fixed backend/data path, so
    the in-image location must stay a valid answer."""
    (fake_roots / "backend" / "data" / "hpo_labels.json").write_text("{}", encoding="utf-8")

    assert source_row(client.get("/api/admin/db-status").json(), "hpo_labels")["present"] is True


def test_db_status_does_not_fall_back_when_an_override_is_set(
    client, logs_dir, as_admin, fake_roots, monkeypatch,
):
    """The readers resolve ``override or default`` — never both. So an override
    naming a file that is not there means the app cannot read it, whatever the
    in-image copy holds. Reporting "present" off a copy the reader will never
    open is the false positive that mirrors the bug this fixes."""
    (fake_roots / "backend" / "data" / "clingen_gene_validity.json").write_text(
        "{}", encoding="utf-8"
    )
    monkeypatch.setenv("CLINGEN_GV_PATH", str(fake_roots / "data" / "not_there.json"))

    row = source_row(client.get("/api/admin/db-status").json(), "clingen_gv")
    assert row["present"] is False


def test_db_status_ignores_heartvar_structure_dir_for_the_alphafold_manifest(
    client, logs_dir, as_admin, fake_roots, monkeypatch,
):
    """HEARTVAR_STRUCTURE_DIR moves the structures app.py *serves*; the manifest
    reader (clients/alphafold.py) resolves data/alphafold/manifest.json and does
    not consult it. So the portal must not either — honouring a variable the
    reader ignores is the same mismatch this fixes, pointing the other way."""
    manifest = fake_roots / "data" / "alphafold" / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("{}", encoding="utf-8")
    monkeypatch.setenv("HEARTVAR_STRUCTURE_DIR", str(fake_roots / "elsewhere"))

    row = source_row(client.get("/api/admin/db-status").json(), "alphafold")
    assert row["present"] is True
    assert row["artifact_path"] == "data/alphafold/manifest.json"


def test_db_status_env_override_map_only_names_known_sources():
    """A typo'd or stale key would silently never be consulted."""
    assert set(app_module._DB_SOURCE_ENV_OVERRIDES) <= set(app_module._DB_SOURCE_ARTIFACTS)


def test_the_deploy_workflow_points_the_gitignored_copies_at_the_mount():
    """hgnc_alias_map.json, erepo_all.tsv and gencc_submissions.json are
    gitignored, so the image has no copy: without the override the reader falls
    back to an in-image path that cannot exist and silently goes to the network.
    GENCC_SNAPSHOT_PATH was missing for exactly this reason. clingen_gv is
    excluded — it is tracked in git, so its in-image copy is real."""
    from .conftest import require_file
    wf = require_file(".github/workflows/restart-webapp.yml").read_text(
        encoding="utf-8"
    )
    for var in ("HGNC_ALIAS_MAP_PATH", "HEARTVAR_EREPO_TSV", "GENCC_SNAPSHOT_PATH"):
        assert f'{var}="/app/data/' in wf, f"{var} is not set to a mount path by the deploy"


def test_db_status_artifact_map_matches_the_cadence_map():
    """Both mirror scripts/build_all.sh. A source in one and not the other
    silently reports the wrong cadence (the .get() default) or is invisible."""
    assert set(app_module._DB_SOURCE_ARTIFACTS) == set(app_module._DB_SOURCE_CADENCE)


_REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent


def test_no_shell_script_carries_a_utf8_bom():
    """A BOM before `#!` stops the kernel recognising the shebang: running the
    script gives "#!/usr/bin/env bash: No such file or directory". It arrives
    silently from Windows editors, and .gitattributes `eol=lf` does not strip it
    — it only normalises line endings. entrypoint_builder.sh and build_all.sh are
    the builder container's entrypoint and the monthly data build, so a BOM there
    breaks the data refresh rather than anything visible in the app."""
    offenders = [
        str(p.relative_to(_REPO_ROOT))
        for p in sorted((_REPO_ROOT / "scripts").rglob("*.sh"))
        if p.read_bytes().startswith(b"\xef\xbb\xbf")
    ]
    assert not offenders, f"UTF-8 BOM found in: {offenders}"


def test_no_shell_script_has_crlf_line_endings():
    """CRLF breaks bash inside the Linux containers the same way — `\\r` becomes
    part of the last token on every line."""
    offenders = [
        str(p.relative_to(_REPO_ROOT))
        for p in sorted((_REPO_ROOT / "scripts").rglob("*.sh"))
        if b"\r\n" in p.read_bytes()
    ]
    assert not offenders, f"CRLF line endings found in: {offenders}"


def test_content_cannot_reach_a_file_outside_the_category(client, logs_dir, as_admin):
    """The filename is matched against a directory listing, so a name that exists
    in ANOTHER allowed category is still not reachable from this one."""
    write_log(logs_dir / "deployments", "elsewhere.log", "OTHER CATEGORY")
    r = client.get(
        "/api/admin/log-content",
        params={"category": "db_builder", "name": "elsewhere.log"},
    )
    assert r.status_code == 404
    assert "OTHER CATEGORY" not in r.text


def test_the_suite_never_writes_deployment_logs_into_the_repo_data_dir():
    """THE 2026-08-28 CONFUSION, pinned.

    app._LOG_BASE_DIR is resolved at IMPORT time from HEARTVAR_LOGS_DIR and
    defaults to <repo>/data/logs — which in a dev checkout is the real ~24 GB
    mirror. So importing the app in a test wrote a live
    data/logs/deployments/<ts>.log, and _record_deploy_event then attaches a
    FileHandler, so the rest of the run's records landed in it too: pytest tmp
    paths, `ip=testclient`, mocked failures like RuntimeError('anthropic is
    down'), log-truncation fixture content.

    The admin portal serves that same directory as "deployments", so five such
    files were read as Azure deployment logs while checking whether a
    redeployment had landed. A test artefact that impersonates production
    evidence is worse than no log at all.

    conftest.py redirects HEARTVAR_LOGS_DIR at import time (it must be import
    time — the target is a module-level constant). This asserts the redirect is
    actually in effect, so deleting it fails here instead of silently polluting
    the mirror again.
    """
    log_base = pathlib.Path(app_module._LOG_BASE_DIR).resolve()
    repo_data = (_REPO_ROOT / "data").resolve()

    assert not log_base.is_relative_to(repo_data), (
        f"the suite writes log files into the real data directory ({log_base}). "
        "conftest.py must set HEARTVAR_LOGS_DIR before backend.app is imported"
    )
