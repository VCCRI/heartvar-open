"""Unit tests for backend.clients._cache.TTLCache.

Covers the four behaviours the external clients rely on: hit-within-TTL,
expiry-then-refetch, single-flight coalescing of concurrent identical calls,
should_cache=False not retaining a value (but still returning it), and
exception propagation without caching.

No pytest dependency for the logic — runnable directly
(``python -m backend.tests.test_cache``). Uses asyncio.run per test.
"""
from __future__ import annotations

import asyncio

import backend.claude as claude
from backend.clients._cache import TTLCache, _ok


def test_hit_within_ttl_calls_factory_once():
    cache = TTLCache("t", default_ttl=100)
    calls = {"n": 0}

    async def factory():
        calls["n"] += 1
        return {"ok": True, "v": calls["n"]}

    async def go():
        a = await cache.get_or_set("k", factory)
        b = await cache.get_or_set("k", factory)
        return a, b

    a, b = asyncio.run(go())
    assert a == b == {"ok": True, "v": 1}
    assert calls["n"] == 1
    assert cache.stats()["hits"] == 1


def test_expired_entry_refetches():
    cache = TTLCache("t", default_ttl=100)
    calls = {"n": 0}

    async def factory():
        calls["n"] += 1
        return calls["n"]

    async def go():
        first = await cache.get_or_set("k", factory, ttl=0)
        second = await cache.get_or_set("k", factory, ttl=0)
        return first, second

    first, second = asyncio.run(go())
    assert (first, second) == (1, 2)
    assert calls["n"] == 2


def test_single_flight_coalesces_concurrent_calls():
    cache = TTLCache("t", default_ttl=100)
    calls = {"n": 0}
    started = asyncio.Event()
    release = asyncio.Event()

    async def slow_factory():
        calls["n"] += 1
        started.set()
        await release.wait()
        return "value"

    async def go():
        leader = asyncio.create_task(cache.get_or_set("k", slow_factory))
        await started.wait()
        f1 = asyncio.create_task(cache.get_or_set("k", slow_factory))
        f2 = asyncio.create_task(cache.get_or_set("k", slow_factory))
        await asyncio.sleep(0)
        release.set()
        return await asyncio.gather(leader, f1, f2)

    results = asyncio.run(go())
    assert results == ["value", "value", "value"]
    assert calls["n"] == 1, "factory must run exactly once for coalesced calls"
    assert cache.stats()["coalesced"] == 2


def test_should_cache_false_returns_but_does_not_store():
    cache = TTLCache("t", default_ttl=100)
    calls = {"n": 0}

    async def factory():
        calls["n"] += 1
        return {"ok": False, "n": calls["n"]}

    async def go():
        a = await cache.get_or_set("k", factory, should_cache=_ok)
        b = await cache.get_or_set("k", factory, should_cache=_ok)
        return a, b

    a, b = asyncio.run(go())
    assert a == {"ok": False, "n": 1}
    assert b == {"ok": False, "n": 2}
    assert calls["n"] == 2
    assert cache.stats()["entries"] == 0


def test_exception_propagates_and_is_not_cached():
    cache = TTLCache("t", default_ttl=100)
    calls = {"n": 0}

    async def boom():
        calls["n"] += 1
        raise ValueError("upstream down")

    async def go():
        errors = 0
        for _ in range(2):
            try:
                await cache.get_or_set("k", boom)
            except ValueError:
                errors += 1
        return errors

    errors = asyncio.run(go())
    assert errors == 2
    assert calls["n"] == 2
    assert cache.stats()["entries"] == 0


def test_maxsize_eviction_bounds_store():
    cache = TTLCache("t", default_ttl=100, maxsize=3)

    async def go():
        for i in range(6):
            await cache.get_or_set(f"k{i}", _const(i), ttl=10 + i)

    asyncio.run(go())
    assert cache.stats()["entries"] == 3


def _const(value):
    async def factory():
        return value
    return factory


# patching the module-level attribute directly (and restoring it) rather than

def test_system_blocks_default_has_no_ttl_key():
    prev = claude._CACHE_TTL
    try:
        claude._CACHE_TTL = ""
        cc = claude._system_blocks()[0]["cache_control"]
        assert cc == {"type": "ephemeral"}
    finally:
        claude._CACHE_TTL = prev


def test_system_blocks_ttl_env_adds_ttl_key():
    prev = claude._CACHE_TTL
    try:
        claude._CACHE_TTL = "1h"
        cc = claude._system_blocks()[0]["cache_control"]
        assert cc == {"type": "ephemeral", "ttl": "1h"}
    finally:
        claude._CACHE_TTL = prev


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
