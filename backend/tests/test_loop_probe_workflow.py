"""Tests for the manually-triggered loop-probe workflow.

WHY A WORKFLOW AT ALL. Arming the probe needs `az webapp config appsettings
set`, and reading the result needs the app's logs. We have NO Azure access —
deploy-doc item A0 (Reader on the resource group + log stream) is still an open
ask to IT. GitHub Actions, however, already authenticates to Azure via OIDC and
is trusted to run `az webapp stop/start` and `az webapp config appsettings set`
(see restart-webapp.yml). So the workflow is the only self-service route.

The contracts here are SAFETY contracts, not conveniences:

  1. Manual only. This restarts production; it must never fire on a push.
  2. It arms HEARTVAR_LOOP_PROBE.
  3. ⚠ IT ALWAYS DISARMS. Debug mode times every callback, so a probe left on is
     a permanent tax on production. The cleanup must run with `if: always()` so
     a failed curation, a timeout or a cancelled job cannot leave it armed.
  4. The curation it fires uses ai_mode=none — no Anthropic spend.
  5. No hardcoded resource group / webapp name / hostname. Those come from
     secrets or from `az`, exactly as the other workflows do.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

_WF = (Path(__file__).resolve().parents[2] / ".github" / "workflows"
       / "loop-probe.yml")
if not _WF.exists():
    pytest.skip("'.github/workflows/loop-probe.yml' is not in this checkout",
                allow_module_level=True)

WORKFLOW = (Path(__file__).resolve().parents[2]
            / ".github" / "workflows" / "loop-probe.yml")


@pytest.fixture(scope="module")
def wf() -> str:
    assert WORKFLOW.is_file(), f"missing workflow: {WORKFLOW}"
    return WORKFLOW.read_text(encoding="utf-8")


def test_manual_trigger_only(wf):
    """It restarts production — it must not be reachable from a push."""
    assert "workflow_dispatch:" in wf
    trigger_block = wf.split("jobs:")[0]
    assert not re.search(r"^\s*push:", trigger_block, re.MULTILINE), (
        "loop-probe must never run on push — it restarts the webapp"
    )
    assert not re.search(r"^\s*schedule:", trigger_block, re.MULTILINE)


def test_it_arms_the_probe(wf):
    assert "HEARTVAR_LOOP_PROBE=1" in wf


def test_it_always_disarms(wf):
    """The one contract that cannot be allowed to regress."""
    assert "appsettings delete" in wf, "must remove the setting, not just unset it"
    assert "HEARTVAR_LOOP_PROBE" in wf
    assert "if: always()" in wf, (
        "the disarm step must be if: always() — a failed or cancelled run must "
        "not leave asyncio debug mode on in production"
    )
    assert wf.index("appsettings delete") > wf.index("HEARTVAR_LOOP_PROBE=1")


def test_the_curation_spends_no_ai_credit(wf):
    """ai_mode=none is the public, deterministic path — no owner key is spent."""
    assert '"ai_mode": "none"' in wf or "'ai_mode': 'none'" in wf


def test_no_hardcoded_azure_identifiers(wf):
    assert "secrets.WEBAPP_NAME" in wf
    assert "secrets.RESOURCE_GROUP_NAME" in wf
    assert "azure/login" in wf
    assert "defaultHostName" in wf, (
        "derive the host from `az webapp show` rather than hardcoding a URL"
    )
    assert "heartvar.victorchang.edu.au" not in wf


def test_it_greps_for_the_verdict(wf):
    """The job log has to actually show the answer."""
    assert "took" in wf and "loop-probe" in wf


def test_the_verdict_thresholds_on_MAGNITUDE_not_existence(wf):
    """⚠ THE BUG THIS PINS, caught on the first real reading (2026-08-30).

    The first version printed "THE LOOP WAS BLOCKED" whenever ANY `took` line
    existed. The real reading had a 20.2 s stall whose largest blocked callback
    was 0.533 s — i.e. the loop was demonstrably NOT the problem — and the step
    still announced that it was. A diagnostic that reports the wrong conclusion
    is worse than one that reports nothing, because it gets acted on.

    So the verdict must compare the LARGEST `took` value against the stall it is
    trying to explain, and must always print that maximum.
    """
    assert "MAX_BLOCK" in wf, "the workflow must compute the largest block"
    assert "VERDICT_MIN_SECONDS" in wf, (
        "the workflow must compare the max against a magnitude floor, not "
        "merely assert that some `took` line exists"
    )
    assert "if grep -aE \"took [0-9.]+ seconds\" /tmp/logtail.txt; then" not in wf


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
