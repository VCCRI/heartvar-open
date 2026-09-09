"""Guards on scripts/build_all.sh — build ORDER and refresh CADENCE.

Both of the failures these cover actually happened, in the 2026-08-01 monthly
refresh (see monthly_db_update.txt):

  1. **Order.** ``build_opentargets_db.py`` reads
     ``backend/data/hgnc_complete_set.txt``, which only
     ``build_hgnc_alias_db.py`` downloads — and it is gitignored AND
     dockerignored, so it does not exist in a fresh container. Open Targets was
     sequenced *before* hgnc_alias, so it died on FileNotFoundError on every
     single build. Nothing caught it because each builder is allowed to fail
     independently (the app degrades to a live fallback).

  2. **Cadence.** The scheduled job ran ``--resume``, which skips every artifact
     already on the mount. Twelve sources that should refresh monthly — ClinVar,
     UniProt and the rest — were reported "PRESENT" and never rebuilt, so the
     mirror silently aged.

These are shell-level bugs in a script no Python test would otherwise touch, and
both are invisible in a passing build log, so they are asserted here rather than
left to review. The cadence table is also the contract behind
the cadence table (kept in the deployment documentation, out of the
repo since 2026-08-27); if the two drift, this file should be the thing that
fails.

The bash functions are evaluated in isolation (extracted from the script, never
sourced) so nothing here can start a real data build.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

BUILD_ALL = Path(__file__).resolve().parents[2] / "scripts" / "build_all.sh"
SCRIPT = BUILD_ALL.read_text()


def _run_labels() -> list[str]:
    """Labels in the order build_all.sh runs them, Python and shell steps alike.

    `run_sh` must be matched too. Without it a shell step escapes
    test_every_run_label_has_a_cadence_and_an_output entirely — exactly the kind
    of quiet gap in the cadence table this file exists to prevent.
    """
    return re.findall(r"^run(?:_sh)?\s+(\w+)", SCRIPT, re.MULTILINE)


def _extract_function(name: str) -> str:
    match = re.search(rf"^{name}\(\)\s*\{{.*?^\}}", SCRIPT, re.MULTILINE | re.DOTALL)
    assert match, f"{name}() not found in build_all.sh"
    return match.group(0)


def _call_bash_function(name: str, arg: str) -> str:
    """Evaluate one extracted function in a throwaway bash and call it.

    Extracting rather than sourcing build_all.sh is deliberate: sourcing would
    execute the script and kick off real multi-hour downloads.
    """
    body = _extract_function(name)
    script = f'DATA_DIR=/d; BACKEND_DATA=/b\n{body}\n{name} "{arg}"'
    out = subprocess.run(
        ["bash", "-c", script], capture_output=True, text=True, timeout=30,
    )
    assert out.returncode == 0, out.stderr
    return out.stdout.strip()


HGNC_DEPENDENTS = ("opentargets", "gene_id_map")


@pytest.mark.parametrize("label", HGNC_DEPENDENTS)
def test_hgnc_dependents_run_after_hgnc_alias(label):
    labels = _run_labels()
    assert label in labels, f"{label} is not wired into build_all.sh"
    assert labels.index(label) > labels.index("hgnc_alias"), (
        f"{label} reads backend/data/hgnc_complete_set.txt, which only "
        "build_hgnc_alias_db.py downloads (it is gitignored AND dockerignored, so "
        "it is absent in a fresh container). Running it before hgnc_alias is the "
        "2026-08-01 FileNotFoundError. Move it after."
    )


def test_uniprot_runs_before_alphafold():
    """AlphaFold structure selection reads data/uniprot.db."""
    labels = _run_labels()
    assert labels.index("uniprot") < labels.index("alphafold")


EXPECTED_CADENCE = {
    "uniprot": "monthly",
    "clinvar": "monthly",
    "medgen": "monthly",
    "mgi": "monthly",
    "panelapp": "monthly",
    "hpo_labels": "monthly",
    "hgnc_alias": "monthly",
    "erepo": "monthly",
    "gencc": "monthly",
    "clingen_gv": "monthly",
    "opentargets": "monthly",
    "biogrid": "monthly",
    "gnomad_constraint": "versioned",
    "gnomad_freq": "versioned",
    "gtex": "versioned",
    "alphafold": "static",
    "spliceai": "static",
    "fetal_heart": "static",
    "gene_id_map": "image",
    "vep": "versioned",
}


@pytest.mark.parametrize("label,expected", sorted(EXPECTED_CADENCE.items()))
def test_cadence_classification(label, expected):
    assert _call_bash_function("cadence_for", label) == expected


def test_clinvar_and_uniprot_are_monthly():
    """The headline regression: these two were skipped on 2026-08-01."""
    for label in ("clinvar", "uniprot"):
        assert _call_bash_function("cadence_for", label) == "monthly"


def test_unknown_label_defaults_to_monthly():
    """A new builder must default to being REFRESHED, not to going stale —
    forgetting a cadence entry should be the safe mistake."""
    assert _call_bash_function("cadence_for", "some_new_source_2027") == "monthly"


def test_gnomad_freq_is_never_monthly():
    """8+ GB and multi-hour for a frozen numbered release. Rebuilding this every
    month is how a monthly job becomes something IT switches off."""
    assert _call_bash_function("cadence_for", "gnomad_freq") != "monthly"


def test_every_run_label_has_a_cadence_and_an_output():
    """Every source must appear explicitly in BOTH tables.

    A missing cadence entry silently falls through to `monthly`, and a missing
    output entry makes --monthly unable to tell "already built" from "absent" —
    so it rebuilds unconditionally. Neither is dangerous, but both mean a source
    was added without anyone deciding how often it should refresh.
    """
    cadence_arms = _extract_function("cadence_for")
    output_arms = _extract_function("output_for")
    for label in _run_labels():
        assert re.search(rf"[\s|(]{re.escape(label)}[|)]", cadence_arms), (
            f"{label} has no explicit cadence_for entry — decide whether it is "
            "monthly / versioned / static / image and add it."
        )
        assert re.search(rf"^\s*{re.escape(label)}\)", output_arms, re.MULTILINE), (
            f"{label} has no output_for entry, so --monthly cannot tell whether "
            "it is already built."
        )


def test_monthly_and_resume_are_distinct_flags():
    """--resume must keep its crash-retry meaning (skip everything present); the
    scheduled refresh is --monthly. Collapsing them recreates the 2026-08-01 no-op."""
    assert "--monthly) MONTHLY=1" in SCRIPT
    assert "--resume) RESUME=1" in SCRIPT


def test_entrypoint_uses_monthly_not_resume():
    """The builder container's start IS the scheduled refresh."""
    entrypoint = (BUILD_ALL.parent / "entrypoint_builder.sh").read_text()
    assert "--monthly" in entrypoint
    assert not re.search(r"build_all\.sh\s+--resume", entrypoint), (
        "entrypoint_builder.sh must not hardcode --resume — that is the "
        "2026-08-01 bug (every monthly source reported PRESENT and was skipped)."
    )


def test_vep_is_opt_in_only():
    """The ~23 GB Ensembl cache must not be startable by a routine build.

    `versioned` cadence is not sufficient on its own: the KEEP depends on the
    manifest existing, the manifest is all-or-nothing, and the per-component
    state is keyed on VEP release AND assembly. So a failed FASTA/REVEL check
    OR a VEP_RELEASE bump would re-enter the step inside what looks like an
    ordinary --monthly run. Requiring --with-vep makes the download a decision.
    """
    assert "--with-vep) WITH_VEP=1" in SCRIPT, "--with-vep must be parseable"
    assert re.search(r'WITH_VEP="\$\{WITH_VEP:-\}"', SCRIPT), (
        "WITH_VEP must also be settable from the environment — the Container "
        "App Job's command is fixed, so env is the only way to opt in there."
    )
    gate = re.search(
        r'if \[\[ -z "\$WITH_VEP" \]\].*?\nfi\n(?=run_sh vep)', SCRIPT, re.DOTALL,
    )
    assert gate, "the opt-in gate must sit immediately before `run_sh vep`"
    assert "vep" in gate.group(0) and "SKIP_LIST" in gate.group(0), (
        "the gate must add vep to SKIP_LIST when WITH_VEP is unset"
    )


def test_vep_is_never_monthly():
    """A calendar month must never re-download a version-locked 23 GB cache."""
    assert _call_bash_function("cadence_for", "vep") == "versioned"


def _lock_is_stale(tmpdir: Path, heartbeat_age_minutes: int | None,
                   stale_minutes: int = 15) -> bool:
    """Drive the real _lock_is_stale() against a temp lock dir.

    heartbeat_age_minutes=None writes no heartbeat file at all.
    """
    import os
    import time

    lock = tmpdir / ".build_all.lock"
    lock.mkdir(exist_ok=True)
    meta = lock / "heartbeat"
    if heartbeat_age_minutes is not None:
        meta.write_text("0\n")
        when = time.time() - heartbeat_age_minutes * 60
        os.utime(meta, (when, when))
    script = "\n".join([
        f'LOCK_META="{meta}"',
        f"LOCK_STALE_MINUTES={stale_minutes}",
        _extract_function("_mtime"),
        _extract_function("_lock_is_stale"),
        "_lock_is_stale && echo STALE || echo LIVE",
    ])
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True, timeout=30)
    return "STALE" in out.stdout


def test_live_heartbeat_is_not_reclaimed(tmp_path):
    """The guard's whole purpose: never steal the lock from a running build."""
    assert _lock_is_stale(tmp_path, heartbeat_age_minutes=0) is False
    assert _lock_is_stale(tmp_path, heartbeat_age_minutes=14) is False


def test_abandoned_heartbeat_is_reclaimed(tmp_path):
    """Past the window → the holder is gone (OOM/eviction); take the lock."""
    assert _lock_is_stale(tmp_path, heartbeat_age_minutes=16) is True
    assert _lock_is_stale(tmp_path, heartbeat_age_minutes=60 * 24) is True


def test_missing_heartbeat_is_reclaimed(tmp_path):
    """A lock with no heartbeat predates this mechanism, or died before its first
    write. An unattributable lock blocking every future run is the worse failure."""
    assert _lock_is_stale(tmp_path, heartbeat_age_minutes=None) is True


def test_lock_does_not_use_pid_or_boot_id_heuristics():
    """Regression guard on a bug in the first version of this fix.

    PIDs are meaningless across containers (separate PID namespaces) while
    /proc/sys/kernel/random/boot_id is SHARED by containers on one kernel — so
    "same boot_id and PID not alive → stale" wrongly reclaims a lock held by a
    live sibling replica, breaking the guard in precisely the two-replica case it
    exists for. Don't reintroduce it.
    """
    lock_section = SCRIPT[SCRIPT.index("LOCK_DIR="):SCRIPT.index("declare -a FAILED")]
    assert "boot_id" not in lock_section
    assert "kill -0" not in lock_section


def test_lock_is_released_on_term_and_int():
    """A graceful stop should free the lock at once, not leave the next run to
    wait out the staleness window."""
    for sig in ("EXIT", "TERM", "INT"):
        assert re.search(rf"trap '_lock_release[^']*' {sig}", SCRIPT), sig


def test_run_and_run_sh_share_one_implementation():
    """Both must delegate to _run_with, so the cadence/skip/KEEP logic exists once."""
    assert re.search(r"^_run_with\(\)", SCRIPT, re.MULTILINE), (
        "_run_with() not found — run() and run_sh() must share their bookkeeping."
    )
    assert "_run_with" in _extract_function("run"), "run() must delegate to _run_with"
    assert "_run_with" in _extract_function("run_sh"), "run_sh() must delegate to _run_with"


def test_run_sh_uses_bash_and_run_uses_python():
    assert re.search(r'_run_with\s+bash\s+"\$@"', _extract_function("run_sh"))
    assert re.search(r'_run_with\s+"\$PY"\s+"\$@"', _extract_function("run"))


def _drive_runner(runner_call: str, exit_code: int) -> tuple[str, str]:
    """Run one _run_with-based step against a stub command; return (stdout, summary).

    The 'script' is a bash function, so nothing real executes.
    """
    script = "\n".join([
        "DATA_DIR=/d; BACKEND_DATA=/b; PY=stub_py; RESUME=; MONTHLY=; SKIP_LIST=",
        "declare -a FAILED=(); declare -a DONE=(); declare -a SKIPPED=()",
        "declare -a PRESENT=(); declare -a REFRESHED=()",
        'skipped() { [[ ",$SKIP_LIST," == *",$1,"* ]]; }',
        'output_for() { echo ""; }',
        "cadence_for() { echo monthly; }",
        f"stub_py() {{ return {exit_code}; }}",
        f"bash() {{ return {exit_code}; }}",
        _extract_function("_run_with"),
        _extract_function("run"),
        _extract_function("run_sh"),
        runner_call,
        'echo "SUMMARY done=${DONE[*]:-} failed=${FAILED[*]:-}"',
    ])
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True, timeout=30)
    summary = [ln for ln in out.stdout.splitlines() if ln.startswith("SUMMARY")]
    return out.stdout, (summary[0] if summary else "")


def test_run_sh_records_success():
    _, summary = _drive_runner("run_sh mylabel /some/script.sh", exit_code=0)
    assert "done=mylabel" in summary
    assert "failed=mylabel" not in summary


def test_run_sh_records_failure_without_aborting():
    """A failed VEP install must not stop the rest of the build."""
    _, summary = _drive_runner("run_sh mylabel /some/script.sh", exit_code=1)
    assert "failed=mylabel" in summary


def test_run_still_uses_the_python_interpreter():
    """The refactor must not change run()'s behaviour for the 19 existing steps."""
    _, summary = _drive_runner("run uniprot /some/build.py", exit_code=0)
    assert "done=uniprot" in summary


def test_run_sh_honours_skip_list():
    script = "\n".join([
        "DATA_DIR=/d; BACKEND_DATA=/b; PY=stub_py; RESUME=; MONTHLY=; SKIP_LIST=vep",
        "declare -a FAILED=(); declare -a DONE=(); declare -a SKIPPED=()",
        "declare -a PRESENT=(); declare -a REFRESHED=()",
        'skipped() { [[ ",$SKIP_LIST," == *",$1,"* ]]; }',
        'output_for() { echo ""; }',
        "cadence_for() { echo versioned; }",
        "bash() { echo SHOULD_NOT_RUN; return 0; }",
        _extract_function("_run_with"),
        _extract_function("run_sh"),
        "run_sh vep /some/script.sh",
        'echo "SUMMARY skipped=${SKIPPED[*]:-}"',
    ])
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True, timeout=30)
    assert "SHOULD_NOT_RUN" not in out.stdout
    assert "skipped=vep" in out.stdout


def test_vep_is_wired_in_as_a_shell_step():
    assert re.search(r"^run_sh\s+vep\s", SCRIPT, re.MULTILINE), (
        "offline VEP must be a real build_all.sh step, not a NOTE telling the "
        "operator to run a separate job by hand."
    )


def test_vep_runs_last():
    """A multi-hour install must not delay the caches curators depend on, and a
    space failure must strand nothing after it."""
    labels = _run_labels()
    assert labels[-1] == "vep", f"vep should be the final step; order is {labels}"


def test_vep_output_is_the_all_or_nothing_manifest():
    """setup_offline_vep.sh writes this only once the cache, FASTA and REVEL have
    each been verified by annotating a real variant. A manifest that appeared
    after a PARTIAL install would make the `versioned` KEEP branch skip the step
    forever with REVEL missing."""
    assert _call_bash_function("output_for", "vep") == "/d/vep/.heartvar_vep_manifest.json"


def test_vep_is_versioned_not_monthly():
    """Version-locked to the vep binary in the image. A monthly rebuild would
    re-download ~27 GB of byte-identical data, which is the fastest way to get the
    monthly job switched off."""
    assert _call_bash_function("cadence_for", "vep") == "versioned"


def test_vep_is_skipped_not_failed_when_the_toolchain_is_absent():
    """On a dev machine INSTALL.pl is absent, and that is not a failure worth
    reporting — a local `bash scripts/build_all.sh` should stay green. In the
    deploy the builder image is built FROM the Ensembl VEP image, so it is always
    present there and a genuine failure still reports as one."""
    assert re.search(r"command -v INSTALL\.pl", SCRIPT)
    assert re.search(r'SKIP_LIST="\$\{SKIP_LIST:\+\$SKIP_LIST,\}vep"', SCRIPT)


def test_the_old_manual_note_is_gone():
    assert "or run scripts/setup_offline_vep.sh on a host that has Docker" not in SCRIPT


def test_free_space_is_logged():
    """build_all.sh recorded `du` (used) but never `df` (free), so the one number
    needed to answer "can the share take the VEP cache?" was absent from every
    build log and had to be asked of IT instead."""
    assert re.search(r"df -h", SCRIPT), "log free space next to the du -sh line"


def _acquire(tmpdir: Path, *, heartbeat_age_minutes=None, lock_exists=True,
             keep_alive=False, wait_seconds=2, timeout=40):
    """Drive the real _acquire_lock() against a temp lock dir."""
    import os
    import time

    lock = tmpdir / ".build_all.lock"
    meta = lock / "heartbeat"
    if lock_exists:
        lock.mkdir(exist_ok=True)
        if heartbeat_age_minutes is not None:
            meta.write_text("0\n")
            when = time.time() - heartbeat_age_minutes * 60
            os.utime(meta, (when, when))

    toucher = ""
    if keep_alive:
        toucher = (f'( for i in $(seq 1 60); do date -u +%s > "{meta}"; '
                   f"sleep 0.3; done ) & TOUCHER=$!\n")

    script = "\n".join([
        f'LOCK_DIR="{lock}"',
        f'LOCK_META="{meta}"',
        "LOCK_STALE_MINUTES=15",
        f"LOCK_WAIT_SECONDS={wait_seconds}",
        "LOCK_POLL_SECONDS=1",
        _extract_function("_mtime"),
        _extract_function("_lock_is_stale"),
        _extract_function("_acquire_lock"),
        toucher,
        "_acquire_lock && echo ACQUIRED || echo ABORTED",
        '[[ -n "${TOUCHER:-}" ]] && kill $TOUCHER 2>/dev/null; true',
    ])
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                         timeout=timeout)
    return out


def test_no_lock_is_acquired_immediately(tmp_path):
    out = _acquire(tmp_path, lock_exists=False)
    assert "ACQUIRED" in out.stdout, out.stderr


def test_dead_holder_is_reclaimed_not_waited_out(tmp_path):
    """Heartbeat past the staleness window: the holder is gone, take the lock."""
    out = _acquire(tmp_path, heartbeat_age_minutes=16)
    assert "ACQUIRED" in out.stdout, out.stderr
    assert "RECLAIMING" in out.stderr


def test_holder_killed_is_waited_out_then_reclaimed(tmp_path):
    """THE regression. A heartbeat that still looks live used to abort the build
    outright; it must now be waited out and reclaimed once it goes stale.

    Parked just UNDER the staleness threshold so the crossing happens within a
    second. Note the real-world consequence of this design: after a hard kill
    (SIGKILL/OOM — a clean stop releases the lock via the TERM trap) the
    replacement waits out the remaining staleness window, up to
    BUILD_ALL_LOCK_STALE_MINUTES, before reclaiming. That is far better than not
    building at all, but it is not instant.
    """
    out = _acquire(tmp_path, heartbeat_age_minutes=14.9, wait_seconds=30)
    assert "ACQUIRED" in out.stdout, f"{out.stdout}\n{out.stderr}"
    assert "RECLAIMING" in out.stderr, out.stderr


def test_a_genuinely_live_sibling_is_still_refused(tmp_path):
    """The guard's whole purpose. Concurrent writers corrupt the SQLite caches, so
    a holder that keeps refreshing its heartbeat must never have its lock stolen."""
    out = _acquire(tmp_path, heartbeat_age_minutes=0, keep_alive=True,
                   wait_seconds=3)
    assert "ABORTED" in out.stdout, f"{out.stdout}\n{out.stderr}"
    assert "replica-count = 1" in out.stderr


def test_lock_released_mid_wait_is_then_acquired(tmp_path):
    """A live holder that finishes normally hands the lock over rather than making
    the next run wait out the full window."""
    lock = tmp_path / ".build_all.lock"
    meta = lock / "heartbeat"
    lock.mkdir()
    meta.write_text("0\n")
    script = "\n".join([
        f'LOCK_DIR="{lock}"', f'LOCK_META="{meta}"',
        "LOCK_STALE_MINUTES=15", "LOCK_WAIT_SECONDS=20", "LOCK_POLL_SECONDS=1",
        _extract_function("_mtime"),
        _extract_function("_lock_is_stale"),
        _extract_function("_acquire_lock"),
        f'( for i in 1 2; do date -u +%s > "{meta}"; sleep 1; done; '
        f'rm -rf "{lock}" ) &',
        "_acquire_lock && echo ACQUIRED || echo ABORTED",
    ])
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                         timeout=40)
    assert "ACQUIRED" in out.stdout, f"{out.stdout}\n{out.stderr}"


def test_the_wait_window_outlasts_the_staleness_window():
    """If the wait were shorter than the staleness window, a dead holder would
    time out as "live" and we would be back to the original bug."""
    wait = re.search(r"LOCK_WAIT_SECONDS=\$\{LOCK_WAIT_SECONDS:-(\d+)\}", SCRIPT)
    stale = re.search(r"LOCK_STALE_MINUTES=\$\{BUILD_ALL_LOCK_STALE_MINUTES:-(\d+)\}",
                      SCRIPT)
    assert wait and stale, "lock timing knobs not found"
    assert int(wait.group(1)) > int(stale.group(1)) * 60, (
        f"wait {wait.group(1)}s must exceed staleness {stale.group(1)}m"
    )


def test_entrypoint_needs_no_retry_loop():
    """build_all.sh owns the waiting, so the entrypoint must not RETRY — a retry
    there would re-run a build against a genuinely live sibling replica.

    The invariant is "no unconditional second build", not "exactly one call". A
    second invocation is allowed when it is scoped with --only, because that is a
    targeted completion pass for one source rather than another whole build. The
    offline-VEP finisher is the case: without it a half-installed VEP (cache
    verified, FASTA/REVEL not) could never complete, since the opt-in gate makes
    every later monthly pass skip it.
    """
    entrypoint = (BUILD_ALL.parent / "entrypoint_builder.sh").read_text()
    invocations = re.findall(
        r'^\s*bash "\$SCRIPT_DIR/build_all\.sh"(.*)$', entrypoint, re.MULTILINE,
    )
    assert invocations, "entrypoint must invoke build_all.sh"
    unscoped = [a for a in invocations if "--only" not in a]
    assert len(unscoped) == 1, (
        "exactly one FULL build_all.sh invocation is allowed; found "
        f"{len(unscoped)}. Extra passes must be scoped with --only."
    )
    assert not re.search(r"^\s*(while|until)\b", entrypoint, re.MULTILINE), (
        "no retry loop here — build_all.sh owns the waiting, and retrying against "
        "a genuinely live sibling replica is what the lock exists to prevent"
    )


def test_entrypoint_refreshes_at_most_once_per_calendar_month():
    """A container start is not the same event as the monthly refresh.

    This entrypoint runs on every container start, which includes every webapp
    restart, so each deploy used to kick off a full ~19-minute rebuild of sources
    that had not moved — four times in one afternoon on 2026-08-25. The decision
    now comes from the mirror's own build_stamp.json, so an off-cycle start is a
    no-op. FORCE_BUILD stays available for a deliberate run.
    """
    entrypoint = (BUILD_ALL.parent / "entrypoint_builder.sh").read_text()
    assert "build_stamp.json" in entrypoint, (
        "the refresh-due decision must read the mirror's stamp, not assume that "
        "booting means a refresh is due"
    )
    assert "last_refresh_utc" in entrypoint
    assert "FORCE_BUILD" in entrypoint, "an off-cycle escape hatch must exist"
    assert re.search(r"sed -n 's/.*last_refresh_utc", entrypoint), (
        "parse the stamp without an interpreter dependency"
    )


def test_entrypoint_vep_finisher_never_starts_a_fresh_install():
    """The finisher must require EXISTING per-component cache state.

    That condition is the whole safety argument: cache-state present means an
    install stalled on a cheap component and finishing it costs ~1.5 GB, whereas
    cache-state absent means nothing is installed and proceeding would start an
    unasked-for ~23 GB download. Losing the second half of the condition turns a
    completion pass into exactly the automatic download the opt-in gate exists to
    prevent.
    """
    entrypoint = (BUILD_ALL.parent / "entrypoint_builder.sh").read_text()
    assert "--only vep --with-vep" in entrypoint, (
        "the completion pass must be scoped to vep and must opt in"
    )
    assert ".heartvar_vep_state/cache.json" in entrypoint, (
        "the finisher must gate on existing per-component cache state, or it "
        "becomes an automatic 23 GB download"
    )
    assert ".heartvar_vep_manifest.json" in entrypoint, (
        "the finisher must not run once the install is already complete"
    )


def test_check_only_and_force_vep_bypass_the_cadence_KEEP_for_vep():
    """THE 2026-08-28 BLOCKER, and it blocked BOTH directions.

    `_run_with`'s KEEP branch runs BEFORE the build script, so with the manifest
    present on the mount a dispatch reported

        >>> KEEP  vep  (versioned upstream — rebuild on a numbered release…)

    and setup_offline_vep.sh was never executed. Consequences:
      * CHECK_ONLY could not probe — its whole job is to report on what is there;
      * and the 26 GB re-install could not be STARTED BY ANY INPUT, because FORCE
        is read by setup_offline_vep.sh, which the KEEP prevented from running.
        The manifest was written over the WRONG cache flavour, so "present" was
        exactly what made it unfixable.

    Narrow on purpose: `vep` only, and on CHECK_ONLY/FORCE_VEP only. A generic
    FORCE would re-enter gnomad_freq (8 GB, multi-hour) inside what looks like a
    routine monthly build.
    """
    body = re.search(r"_run_with\(\) \{(.*?)\n\}", SCRIPT, re.S)
    assert body, "could not isolate _run_with"
    body = body.group(1)

    assert "FORCE_VEP" in body, (
        "nothing in _run_with reaches past the cadence KEEP, so a wrong VEP "
        "install cannot be replaced"
    )
    assert "CHECK_ONLY" in body, (
        "CHECK_ONLY does not bypass the KEEP, so the probe cannot probe"
    )
    assert re.search(r'"\$label" == vep', body), (
        "the KEEP bypass is not scoped to the vep label; a generic bypass would "
        "re-download every versioned source including gnomad_freq"
    )
    assert re.search(r'if \[\[ -z "\$force_this"', body), (
        "the KEEP branch does not consult the bypass, so setting it changes nothing"
    )


def test_a_generic_FORCE_does_not_reenter_every_versioned_source():
    """The guard on the guard. `FORCE` is what setup_offline_vep.sh reads for its
    own components; if build_all.sh ever keys the KEEP bypass on that name
    instead, one input starts an 8 GB gnomad_freq rebuild as a side effect."""
    body = re.search(r"_run_with\(\) \{(.*?)\n\}", SCRIPT, re.S).group(1)
    assert not re.search(r'\$\{FORCE:-\}', body), (
        "the KEEP bypass keys on the generic FORCE; use FORCE_VEP so the blast "
        "radius stays on the vep step"
    )


def test_an_unreadable_free_space_report_is_not_silent():
    """`df -h` returned NOTHING for the Azure Files mount on 2026-08-28, so the
    build printed an empty "data dir free:" section. A blank free-space report is
    indistinguishable from a healthy one, and free space is the guard that stops
    a 26 GB install failing halfway."""
    assert "free space UNKNOWN" in SCRIPT, (
        "a df that produces no output is reported as nothing at all"
    )


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
