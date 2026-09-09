"""Guards on scripts/entrypoint_builder.sh — when the container exits, and when
it must not.

THE BUG THESE PIN: the image runs in two places that want opposite lifetimes. As
a webapp sidecar Azure requires the container to keep running, so the entrypoint
ended in `exec sleep infinity`. But every real data build runs it as a Container
App Job, where the container is expected to EXIT — so the execution sat in
`Running` forever and the 2026-08-25 06:43 build was reported as
`timed out waiting for data-build-job-exec1` five hours after it had
in fact finished at 06:44:51. No amount of extra headroom on the poll deadline
could have fixed that; the wait could never succeed.

The entrypoint is copied into a temp tree with a STUB build_all.sh so nothing
here can start a real build. SCRIPT_DIR is derived from BASH_SOURCE, so the copy
resolves its data dir and its sibling script inside the temp tree. `sleep` is
stubbed too — the sidecar path ends in `exec sleep infinity`, and BSD sleep on a
macOS dev machine rejects `infinity` outright, so a test that waited for the
process to hang would pass on the deploy target and fail on the laptop. Asserting
which branch was TAKEN is the portable question anyway.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
ENTRYPOINT = ROOT / "scripts" / "entrypoint_builder.sh"

def _tree(tmp_path: Path, *, build_all_exit: int = 0,
          vep_exit: int | None = None) -> Path:
    """A temp scripts/ + data/ tree with the real entrypoint and a stub builder.

    The stub records every invocation so a test can assert WHICH pass ran, and
    exits with the code the test asked for.
    """
    scripts = tmp_path / "scripts"
    scripts.mkdir()
    shutil.copy(ENTRYPOINT, scripts / "entrypoint_builder.sh")
    (tmp_path / "data").mkdir()

    vep_line = (
        f'if [[ "$*" == *"--only vep"* ]]; then exit {vep_exit}; fi\n'
        if vep_exit is not None else ""
    )
    stub = scripts / "build_all.sh"
    stub.write_text(
        "#!/usr/bin/env bash\n"
        f'echo "STUB build_all.sh $*" >> "{tmp_path}/calls.log"\n'
        + vep_line
        + f"exit {build_all_exit}\n"
    )
    stub.chmod(0o755)

    bindir = tmp_path / "bin"
    bindir.mkdir()
    sleeper = bindir / "sleep"
    sleeper.write_text(
        "#!/usr/bin/env bash\n"
        f'echo "STUB sleep $*" >> "{tmp_path}/calls.log"\nexit 0\n'
    )
    sleeper.chmod(0o755)
    return scripts


def _run(scripts: Path, env: dict | None = None, timeout: int = 30):
    bindir = scripts.parent / "bin"
    full = {"PATH": f"{bindir}:{os.environ['PATH']}"}
    full.update(env or {})
    return subprocess.run(
        ["bash", str(scripts / "entrypoint_builder.sh")],
        capture_output=True, text=True, timeout=timeout, env=full,
    )


def _calls(tmp_path: Path) -> str:
    log = tmp_path / "calls.log"
    return log.read_text() if log.exists() else ""


def _stamp_this_month(tmp_path: Path) -> None:
    """A build stamp recording a refresh in the current calendar month, so the
    monthly pass is not due."""
    now = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    (tmp_path / "data" / "build_stamp.json").write_text(
        '{\n  "last_run_utc": "%s",\n  "last_refresh_utc": "%s"\n}\n' % (now, now)
    )


@pytest.mark.parametrize("job_env", [
    {"CONTAINER_APP_JOB_NAME": "data-build-job"},
    {"CONTAINER_APP_JOB_EXECUTION_NAME": "data-build-job-exec1"},
    {"EXIT_WHEN_DONE": "1"},
])
def test_exits_instead_of_sleeping_when_run_as_a_job(tmp_path, job_env):
    """A Container App Job execution stays `Running` until its container exits.
    Sleeping forever is what made a finished build look like a timeout."""
    scripts = _tree(tmp_path)
    _stamp_this_month(tmp_path)
    res = _run(scripts, job_env)
    assert res.returncode == 0, res.stdout + res.stderr
    assert "Running as a Container App Job" in res.stdout, res.stdout
    assert "Entering idle sleep" not in res.stdout, res.stdout
    assert "STUB sleep" not in _calls(tmp_path), _calls(tmp_path)


def test_still_sleeps_when_it_is_a_sidecar(tmp_path):
    """The sidecar deployment genuinely needs this — Azure kills a sidecar that
    exits. Fixing the job case must not break it."""
    scripts = _tree(tmp_path)
    _stamp_this_month(tmp_path)
    res = _run(scripts, {})
    assert "Entering idle sleep" in res.stdout, res.stdout
    assert "STUB sleep infinity" in _calls(tmp_path), _calls(tmp_path)


def test_stay_alive_wins_over_job_detection(tmp_path):
    """An escape hatch for the day the platform stops injecting those vars."""
    scripts = _tree(tmp_path)
    _stamp_this_month(tmp_path)
    res = _run(scripts, {"CONTAINER_APP_JOB_NAME": "j", "STAY_ALIVE": "1"})
    assert "STAY_ALIVE set" in res.stdout, res.stdout
    assert "STUB sleep infinity" in _calls(tmp_path), _calls(tmp_path)


def test_a_failed_data_build_exits_nonzero(tmp_path):
    scripts = _tree(tmp_path, build_all_exit=1)
    res = _run(scripts, {"EXIT_WHEN_DONE": "1"})
    assert res.returncode == 1, res.stdout + res.stderr


def test_a_failed_vep_completion_pass_exits_nonzero(tmp_path):
    """The exact state on the mount: cache verified, manifest absent, so the
    finisher runs. If its failure does not reach the container's exit status the
    execution reports success while offline VEP is still broken."""
    scripts = _tree(tmp_path, build_all_exit=0, vep_exit=1)
    _stamp_this_month(tmp_path)
    state = tmp_path / "data" / "vep" / ".heartvar_vep_state"
    state.mkdir(parents=True)
    (state / "cache.json").write_text('{"component": "cache"}\n')
    res = _run(scripts, {"EXIT_WHEN_DONE": "1"})
    assert "--only vep" in _calls(tmp_path), _calls(tmp_path)
    assert res.returncode == 1, res.stdout + res.stderr


def test_the_sidecar_does_not_attempt_the_vep_finisher(tmp_path):
    """The sidecar boots on every webapp restart, and every push to main rebuilds
    the image and restarts the webapp. Attempting an install there turned each
    deploy into a failed vep step and a failure email — two of 2026-08-26's came
    from pushes rather than from anyone asking for a build."""
    scripts = _tree(tmp_path, vep_exit=1)
    _stamp_this_month(tmp_path)
    state = tmp_path / "data" / "vep" / ".heartvar_vep_state"
    state.mkdir(parents=True)
    (state / "cache.json").write_text('{"component": "cache"}\n')
    res = _run(scripts, {})
    assert "--only vep" not in _calls(tmp_path), _calls(tmp_path)
    assert "this is the webapp" in res.stdout, res.stdout
    assert "Entering idle sleep" in res.stdout, res.stdout


def test_the_job_still_finishes_a_partial_install(tmp_path):
    """The gate must not disable the one path that can complete it."""
    scripts = _tree(tmp_path, vep_exit=0)
    _stamp_this_month(tmp_path)
    state = tmp_path / "data" / "vep" / ".heartvar_vep_state"
    state.mkdir(parents=True)
    (state / "cache.json").write_text('{"component": "cache"}\n')
    res = _run(scripts, {"CONTAINER_APP_JOB_NAME": "data-build-job"})
    assert "--only vep" in _calls(tmp_path), _calls(tmp_path)
    assert res.returncode == 0, res.stdout + res.stderr


def test_a_refresh_already_done_this_month_is_skipped(tmp_path):
    """A container start is not a refresh, and a container start is also every
    webapp restart — four of them in one afternoon on 2026-08-25."""
    scripts = _tree(tmp_path)
    _stamp_this_month(tmp_path)
    res = _run(scripts, {"EXIT_WHEN_DONE": "1"})
    assert "SKIPPING the data refresh" in res.stdout, res.stdout
    assert "--monthly" not in _calls(tmp_path), _calls(tmp_path)


def test_a_first_install_is_never_started_automatically(tmp_path):
    """No cache state means the ~23 GB download has never run. Finishing a
    partial install is safe; starting one behind someone's back is not."""
    scripts = _tree(tmp_path)
    _stamp_this_month(tmp_path)
    res = _run(scripts, {"EXIT_WHEN_DONE": "1"})
    assert "--only vep" not in _calls(tmp_path), _calls(tmp_path)
    assert "offline VEP is NOT installed" in res.stdout, res.stdout


@pytest.mark.parametrize("request_env", [
    {"CHECK_ONLY": "1"},
    {"ONLY": "vep"},
    {"CHECK_ONLY": "1", "ONLY": "vep"},
])
def test_an_explicit_request_overrides_the_once_a_month_guard(tmp_path, request_env):
    """CHECK_ONLY or ONLY means somebody asked for something specific. The
    calendar has no opinion about that."""
    scripts = _tree(tmp_path)
    _stamp_this_month(tmp_path)
    res = _run(scripts, {"EXIT_WHEN_DONE": "1", **request_env})
    assert "the once-per-calendar-month guard does not apply" in res.stdout, res.stdout
    assert "SKIPPING the data refresh" not in res.stdout, res.stdout
    assert "--monthly" in _calls(tmp_path), _calls(tmp_path)


def test_a_skipped_build_is_never_reported_as_a_successful_one(tmp_path):
    """"completed successfully" must mean it RAN and succeeded. Printing it after
    the skip branch is what made the no-op dispatch read as a green build."""
    scripts = _tree(tmp_path)
    _stamp_this_month(tmp_path)
    res = _run(scripts, {"EXIT_WHEN_DONE": "1"})
    assert "SKIPPING the data refresh" in res.stdout, res.stdout
    assert "completed successfully" not in res.stdout, (
        "the entrypoint claims success for a build it never ran:\n" + res.stdout
    )
    assert "was NOT RUN" in res.stdout, res.stdout


def test_check_only_never_triggers_the_vep_completion_pass(tmp_path):
    """CHECK_ONLY means install nothing. The completion pass is an INSTALL, and
    build_all.sh has already probed vep on the same pass, so running it here
    would both disobey the flag and probe twice."""
    scripts = _tree(tmp_path)
    _stamp_this_month(tmp_path)
    state = tmp_path / "data" / "vep" / ".heartvar_vep_state"
    state.mkdir(parents=True)
    (state / "cache.json").write_text('{"component": "cache"}')
    res = _run(scripts, {"EXIT_WHEN_DONE": "1", "CHECK_ONLY": "1", "ONLY": "vep"})
    assert "not attempting the offline-VEP" in res.stdout, res.stdout
    assert "--only vep" not in _calls(tmp_path), _calls(tmp_path)


def test_a_plain_container_start_is_still_skipped(tmp_path):
    """The regression direction. Bypassing the guard for explicit requests must
    not bypass it for the case it exists for — a webapp restart, which happens on
    every push to main."""
    scripts = _tree(tmp_path)
    _stamp_this_month(tmp_path)
    res = _run(scripts, {"EXIT_WHEN_DONE": "1"})
    assert "SKIPPING the data refresh" in res.stdout, res.stdout
    assert "--monthly" not in _calls(tmp_path), _calls(tmp_path)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))


def test_a_sidecar_start_ignores_a_stray_only(tmp_path):
    """A webapp container start is the scheduled refresh and must cover the whole
    monthly set, whatever is left in the app settings."""
    scripts = _tree(tmp_path)
    res = _run(scripts, {"ONLY": "uniprot"})
    assert "IGNORING ONLY='uniprot'" in res.stdout, res.stdout
    calls = _calls(tmp_path)
    assert "--monthly" in calls, calls
    assert "uniprot" not in calls, (
        "the stray ONLY reached build_all.sh anyway: " + calls
    )


def test_a_sidecar_start_ignores_a_stray_skip(tmp_path):
    """Same for SKIP — a scheduled refresh that quietly omits ClinVar is the
    failure this guard exists to stop."""
    scripts = _tree(tmp_path)
    res = _run(scripts, {"SKIP": "clinvar"})
    assert "IGNORING" in res.stdout and "clinvar" in res.stdout, res.stdout
    assert "clinvar" not in _calls(tmp_path), _calls(tmp_path)


def test_a_job_dispatch_keeps_its_only(tmp_path):
    """The 2026-08-28 contract is untouched: a job execution asked for something
    specific and gets it."""
    scripts = _tree(tmp_path)
    res = _run(scripts, {"EXIT_WHEN_DONE": "1", "ONLY": "vep"})
    assert "IGNORING ONLY" not in res.stdout, res.stdout
