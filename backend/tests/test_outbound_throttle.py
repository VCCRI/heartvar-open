"""Unit tests for the process-wide outbound per-host throttle.

Guards the ban-prevention guarantee: the outbound rate to any single public
host is bounded regardless of inbound traffic, while a lone request and
distinct hosts are not slowed. Assertions are LOWER bounds on elapsed time
(the throttle enforces a *minimum* spacing — a slow machine can be slower but
never faster), so they are robust against CI timing jitter.

No pytest dependency — runnable with pytest or directly
(``python -m backend.tests.test_outbound_throttle``).
"""
from __future__ import annotations

import asyncio
import time

from backend.clients import _throttle
from backend.clients._throttle import throttle_host, _DEFAULT_MIN_INTERVAL


def test_single_call_is_not_blocked():
    async def go():
        t0 = time.monotonic()
        await throttle_host("https://unique-host-single.example/api")
        return time.monotonic() - t0
    assert asyncio.run(go()) < _DEFAULT_MIN_INTERVAL


def test_concurrent_calls_to_same_host_are_spaced():
    n = 4
    host = "https://unique-host-spaced.example/api"

    async def go():
        t0 = time.monotonic()
        await asyncio.gather(*[throttle_host(host) for _ in range(n)])
        return time.monotonic() - t0

    elapsed = asyncio.run(go())
    assert elapsed >= (n - 1) * _DEFAULT_MIN_INTERVAL * 0.9, elapsed


def test_distinct_hosts_do_not_block_each_other():
    async def go():
        t0 = time.monotonic()
        await asyncio.gather(
            throttle_host("https://unique-host-a.example/x"),
            throttle_host("https://unique-host-b.example/y"),
            throttle_host("https://unique-host-c.example/z"),
        )
        return time.monotonic() - t0
    assert asyncio.run(go()) < _DEFAULT_MIN_INTERVAL


def test_known_hosts_have_tighter_or_configured_intervals():
    assert _throttle._HOST_MIN_INTERVAL["rest.ensembl.org"] < _DEFAULT_MIN_INTERVAL
    assert _throttle._HOST_MIN_INTERVAL["gnomad.broadinstitute.org"] == _DEFAULT_MIN_INTERVAL


def test_blank_host_is_a_noop():
    asyncio.run(throttle_host("not-a-url"))


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
