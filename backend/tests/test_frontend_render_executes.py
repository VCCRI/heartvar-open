"""Runs the frontend render path for real, in node.

WHY THIS FILE EXISTS. The same-site UI shipped with six unqualified
``esc(...)`` calls inside ``_renderEvidence``. ``esc`` is function-local in
eight other renderers and does not exist at script scope, so the page threw
``ReferenceError: esc is not defined`` before painting anything — a blank page,
no error boundary, no message. Every check that was supposed to guard it passed:

  * the eight assertions in ``test_same_site_markup.py`` are Python substring
    searches over the TEXT of heartvar.js, and the text was correct;
  * ``node --check`` parses but does not evaluate, and a ReferenceError is
    perfectly valid syntax.

A substring test cannot catch an unresolved identifier and no number of them
ever will. This harness evaluates the whole of heartvar.js in a node ``vm``
context with a minimal DOM and then CALLS the renderers, so an identifier that
does not resolve here does not resolve in a browser either.

It is not a substitute for a real browser pass — there is no layout, no CSS and
no paint here, so it cannot see the row-height and colour defects that a
Playwright run found. It catches the class of failure that makes the page blank.

Skips (rather than fails) when node is unavailable, so the suite still runs on
a machine without it.
"""
from __future__ import annotations

import os
import pathlib
import shutil
import subprocess

import pytest

_HARNESS = pathlib.Path(__file__).resolve().parent / "js" / "render_smoke.js"
_NODE = shutil.which("node")


def test_harness_file_is_present():
    assert _HARNESS.is_file(), f"missing render harness at {_HARNESS}"


@pytest.mark.skipif(_NODE is None, reason="node is not installed")
def test_render_path_executes_without_throwing():
    """_renderEvidence, _renderGeneContextTab and _buildLandscapeLollipop must
    each run against realistic payloads and return markup.

    The harness's own assertions cover the specific regressions: the
    numbering-caveat row rendering on non-MANE input, the codon label not
    carrying nested-parenthesised provenance, the PM5-ineligible chip not
    wearing the benign class, star-driven mark opacity, and the proband's own
    allele not appearing under "other alleles at this base".
    """
    proc = subprocess.run(
        [_NODE, str(_HARNESS)],
        capture_output=True, text=True, timeout=120,
    )
    if proc.returncode != 0:
        pytest.fail(
            "frontend render harness failed\n"
            f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
        )
    assert "OK" in proc.stdout, proc.stdout


@pytest.mark.skipif(_NODE is None, reason="node is not installed")
def test_harness_actually_detects_an_unresolved_identifier(tmp_path):
    """Guards the guard. If the sandbox ever leaked a global ``esc`` (or the
    harness stopped executing the render path at all), the test above would go
    green on the very bug it exists for. So: inject the original defect into a
    COPY of heartvar.js (via HEARTVAR_JS_PATH), point the harness at it, and
    require a failure. The tracked file is never written.
    """
    js = pathlib.Path(__file__).resolve().parent.parent.parent / "static" / "heartvar.js"
    src = js.read_text(encoding="utf-8")
    assert "_hvEsc(a.ref" in src, "the line the regression is injected into moved"
    broken = src.replace("_hvEsc(a.ref", "esc(a.ref", 1)
    copy = pathlib.Path(tmp_path) / "heartvar_broken.js"
    copy.write_text(broken, encoding="utf-8")
    proc = subprocess.run(
        [_NODE, str(_HARNESS)], capture_output=True, text=True, timeout=120,
        env={**os.environ, "HEARTVAR_JS_PATH": str(copy)},
    )
    assert js.read_text(encoding="utf-8") == src, (
        "the tracked static/heartvar.js was modified by this test"
    )
    assert proc.returncode != 0, (
        "the harness passed with `esc is not defined` reintroduced — it is no "
        "longer executing the render path, or the sandbox is leaking a global "
        "`esc`.\n" + proc.stdout
    )
    assert "esc is not defined" in (proc.stdout + proc.stderr), proc.stdout + proc.stderr


def test_no_bare_esc_call_outside_a_function_that_defines_one():
    """Static companion to the harness above.

    The harness only catches unresolved identifiers on code paths it actually
    executes. `esc` is function-local in eight renderers and does not exist at
    script scope, so a call in any OTHER function is a guaranteed
    ReferenceError the moment that branch runs — including branches no fixture
    reaches. This walks every `esc(` call site, finds its enclosing top-level
    function, and requires that function to declare its own `esc`.

    Cheap, total, and specific to the mistake that blanked the page.
    """
    import re

    js = pathlib.Path(__file__).resolve().parent.parent.parent / "static" / "heartvar.js"
    lines = js.read_text(encoding="utf-8").splitlines()

    fn_start = re.compile(r"^ {4,5}(?:async )?function ([A-Za-z_$][\w$]*)\s*\(")
    esc_def = re.compile(r"\bconst esc\s*=")
    esc_call = re.compile(r"(?<![A-Za-z0-9_$])esc\s*\(")

    bounds = [(i + 1, m.group(1))
              for i, line in enumerate(lines) if (m := fn_start.match(line))]

    def owner(lineno: int) -> str:
        current = "<script top level>"
        for start, name in bounds:
            if start <= lineno:
                current = name
            else:
                break
        return current

    definers = {owner(i + 1) for i, line in enumerate(lines) if esc_def.search(line)}
    offenders = [
        (i + 1, owner(i + 1), line.strip()[:80])
        for i, line in enumerate(lines)
        if esc_call.search(line) and not esc_def.search(line)
        and owner(i + 1) not in definers
    ]
    assert not offenders, (
        "bare esc(...) called in a function that does not define it — this is a "
        "ReferenceError at render time (the module-level escaper is _hvEsc):\n"
        + "\n".join(f"  heartvar.js:{ln} in {fn}: {src}" for ln, fn, src in offenders)
    )
