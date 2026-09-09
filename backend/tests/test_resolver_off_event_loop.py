"""A static guard: the HGVS resolver must never be called on the event loop.

WHY A SOURCE-TEXT TEST AND NOT A BEHAVIOURAL ONE. ``hgvs_resolver.resolve()`` is
synchronous and does real I/O against the SMB mount — the cdot SQLite, the
bgzipped toplevel FASTA, its .fai and .gzi. Called bare from ``async def
fetch_hgvs`` it blocks the whole event loop, which stalls EVERY other source in
the parallel gather rather than just the VEP one. That is the shape the
production logs kept showing: a dozen sources on three continents completing
inside the same few milliseconds, because none of them could make progress until
this one call returned.

``loop_probe.py:26`` named this exact call as suspect #1 on 2026-08-30, and it
still survived the PR #53–#58 sweep that moved 25 other call sites onto the
local-I/O pool. A blocking call is invisible in a unit test — the awaited result
is identical either way — so the only thing that catches a regression here is
looking at the source. The same reasoning as
``test_connect_timeout_guard.py``.

The rule: any call to ``resolve(`` in async application code goes through
``run_local``.
"""
from __future__ import annotations

import ast
import pathlib

import pytest

BACKEND = pathlib.Path(__file__).resolve().parent.parent
RESOLVE_NAMES = {"resolve"}


def _offending_calls(root: pathlib.Path) -> list[str]:
    """Bare ``…resolve(…)`` calls on hgvs_resolver that are not wrapped.

    A wrapped call looks like ``run_local(hgvs_resolver.resolve, …)`` — the
    function is passed as a VALUE, so it never appears as the callee of a Call
    node. That asymmetry is what makes this checkable at all: an offending call
    is syntactically a call, a correct one is syntactically an argument.
    """
    offenders: list[str] = []
    for path in sorted(root.rglob("*.py")):
        if "tests" in path.parts:
            continue
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except (OSError, SyntaxError):
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            # i.e. an attribute call whose attr is `resolve` on a name that
            if (isinstance(func, ast.Attribute)
                    and func.attr in RESOLVE_NAMES
                    and isinstance(func.value, ast.Name)
                    and "resolver" in func.value.id):
                offenders.append(
                    f"{path.relative_to(root)}:{node.lineno}: "
                    f"{func.value.id}.{func.attr}(...) called directly")
    return offenders


def test_the_resolver_is_never_called_on_the_event_loop():
    offenders = _offending_calls(BACKEND)
    assert not offenders, (
        "these calls run SMB-backed SQLite and FASTA reads on the event loop, "
        "which stalls every other source in the gather. Wrap as "
        "`await run_local(hgvs_resolver.resolve, ...)`:\n  "
        + "\n  ".join(offenders)
    )


def test_the_guard_can_actually_see_a_violation(tmp_path):
    """A guard that cannot fail is not a guard."""
    (tmp_path / "bad.py").write_text(
        "from . import hgvs_resolver\n"
        "async def f(h, t, g):\n"
        "    return hgvs_resolver.resolve(h, t, g)\n")
    offenders = _offending_calls(tmp_path)
    assert offenders, "the guard failed to flag a bare resolver call"
    assert "bad.py:3" in offenders[0]


def test_the_guard_accepts_the_wrapped_form(tmp_path):
    (tmp_path / "ok.py").write_text(
        "from . import hgvs_resolver\n"
        "from ..localio import run_local\n"
        "async def f(h, t, g):\n"
        "    return await run_local(hgvs_resolver.resolve, h, t, g)\n")
    assert _offending_calls(tmp_path) == []


def test_the_real_call_site_is_awaited_through_run_local():
    """Belt and braces on the one call site that matters, by name, so a
    refactor that renames the wrapper does not silently pass the AST check."""
    src = (BACKEND / "clients" / "vep_offline.py").read_text(encoding="utf-8")
    assert "run_local(\n        hgvs_resolver.resolve" in src or \
        "run_local(hgvs_resolver.resolve" in src, (
            "vep_offline.fetch_hgvs must dispatch the resolver to the "
            "local-I/O pool — see loop_probe.py:26"
        )


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
