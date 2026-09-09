"""No outbound HTTP call may inherit its connect budget from a read timeout.

⚠ THE BUG THIS PREVENTS, and it cost three deploys and two wrong fixes to find.
Every affected client passed a SCALAR timeout — medgen 20/25s, pubmed 20/25s,
uniprot 20s, gtex 30s, mgi 20/30s, ensembl 25/30s, spliceai 25s — and httpx
applies a scalar to ALL phases, connect included. On a cold container one attempt
stalled up to 20 s inside connect, so:

  * ``max_attempts`` never mattered — attempt 1 had not returned;
  * the retry deadline never mattered — it is checked BETWEEN attempts.

Which is why the first-curation cost sat at a near-constant 18.00 / 19.63 /
18.30 s across three deploys while the netwarm DNS total swung 27s -> 39s, and
why a DNS server change and a retry deadline both did nothing.

``request_with_retry`` now normalises a scalar centrally. But SEVEN call sites
bypass it and call httpx directly — including ``_ncbi_throttle``, which is the
path for pubmed, pubtator3, medgen and gene_literature, four of the affected
sources. A fix that misses those is not a fix, so this test enumerates every
direct call and fails on any bare scalar.

A read of 25 s is legitimate (pubtator3 routinely takes 2-5 s). A connect of 25 s
never is: measured across all 24 outbound hosts, connect was 1-350 ms.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

CLIENTS = Path(__file__).resolve().parents[1] / "clients"
_HTTP_METHODS = {"get", "post", "put", "request", "stream", "delete", "head"}
_WRAPPER = "_with_connect_cap"
_EXEMPT = {"_http_retry.py"}


def _offending_calls(root: Path | None = None) -> list[str]:
    """Direct httpx calls whose ``timeout=`` is a bare scalar or bare name."""
    bad: list[str] = []
    for path in sorted((root or CLIENTS).glob("*.py")):
        if root is None and path.name in _EXEMPT:
            continue
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except SyntaxError:  # pragma: no cover
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            fn = node.func
            if not (isinstance(fn, ast.Attribute) and fn.attr in _HTTP_METHODS):
                continue
            for kw in node.keywords:
                if kw.arg != "timeout":
                    continue
                val = kw.value
                if isinstance(val, ast.Call):
                    src = ast.unparse(val)
                    if _WRAPPER in src or "Timeout(" in src:
                        continue
                bad.append(
                    f"{path.name}:{node.lineno}  timeout={ast.unparse(val)}")
    return bad


def test_no_direct_http_call_passes_a_bare_scalar_timeout():
    offenders = _offending_calls()
    assert not offenders, (
        "these calls let CONNECT inherit a multi-second read budget — wrap the "
        "timeout in _with_connect_cap(...):\n  " + "\n  ".join(offenders)
    )


def test_the_guard_can_actually_see_a_violation(tmp_path):
    """A guard that cannot fail is not a guard."""
    (tmp_path / "x.py").write_text(
        "async def f(client):\n"
        "    return await client.get('https://x/', timeout=25.0)\n")
    assert _offending_calls(tmp_path), (
        "the guard failed to flag a bare scalar timeout"
    )


def test_the_guard_accepts_the_two_correct_forms(tmp_path):
    (tmp_path / "ok.py").write_text(
        "async def f(client, t):\n"
        "    await client.get('https://x/', timeout=_with_connect_cap(t))\n"
        "    await client.get('https://y/', timeout=httpx.Timeout(9.0, connect=3.0))\n")
    assert _offending_calls(tmp_path) == []

if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
