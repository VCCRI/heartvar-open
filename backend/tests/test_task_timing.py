"""The gather's per-source timing instrument, tested against the bug it exists
to prevent.

On 2026-08-30 the per-source timing log claimed fourteen sources each took
19.37 s, agreeing to within 10 ms, and three deploys were aimed at that number
before anyone checked what it measured. The number was orchestrator HARVEST
latency: ``evidence.py`` recorded elapsed time when ``asyncio.wait(...,
FIRST_COMPLETED)`` handed a finished task back, so a batch collected in one
pass all shared a single ``perf_counter()`` reading.

The first test reproduces that exact shape — three tasks of visibly different
durations, collected late in one batch — and asserts the two measurements
DISAGREE in the documented way. It computes the OLD number too and asserts that
it collapses: that is the control, without which "the timings look plausible"
would have passed against the broken code.

Stdlib only, deliberately, so the instrument can be validated without standing
up the app. See backend/task_timing.py.
"""
from __future__ import annotations

import asyncio
from time import perf_counter

import pytest

from backend.task_timing import track_task_timing

FAST_S = 0.02
MEDIUM_S = 0.10
SLOW_S = 0.40
COLLECTOR_STALL_S = 0.70

TOLERANCE_S = 0.25


async def _sleep_then(name: str, seconds: float) -> str:
    await asyncio.sleep(seconds)
    return name


def _harvest_the_broken_way(
    done: set[asyncio.Task],
    tasks: dict[str, asyncio.Task],
    shared_start: float,
) -> dict[str, float]:
    """The pre-fix measurement, kept as the control.

    One ``perf_counter()`` per harvested BATCH, against one shared start time —
    the two halves of the original bug (evidence.py:1463 and :1228).
    """
    harvest_elapsed: dict[str, float] = {}
    for finished in done:
        src = next(k for k, t in tasks.items() if t is finished)
        harvest_elapsed[src] = perf_counter() - shared_start
    return harvest_elapsed


def test_true_elapsed_survives_a_late_batch_harvest():
    """Three tasks, three durations, collected late in one pass.

    ``true_elapsed`` must report the three real durations. The harvest number
    must collapse them — that is what made 14 sources agree to within 10 ms.
    """
    async def scenario():
        started: dict[str, float] = {}
        true_elapsed: dict[str, float] = {}
        shared_start = perf_counter()

        tasks = {
            name: track_task_timing(
                name,
                asyncio.ensure_future(_sleep_then(name, seconds)),
                started,
                true_elapsed,
            )
            for name, seconds in (
                ("fast", FAST_S), ("medium", MEDIUM_S), ("slow", SLOW_S),
            )
        }

        await asyncio.sleep(COLLECTOR_STALL_S)
        done, pending = await asyncio.wait(
            tasks.values(), return_when=asyncio.FIRST_COMPLETED,
        )
        harvest_elapsed = _harvest_the_broken_way(done, tasks, shared_start)
        return started, true_elapsed, harvest_elapsed, done, pending

    started, true_elapsed, harvest_elapsed, done, pending = asyncio.run(scenario())

    assert not pending, "all three should have finished during the stall"
    assert len(done) == 3, "the point of the test is a multi-task batch"

    harvest_spread = max(harvest_elapsed.values()) - min(harvest_elapsed.values())
    assert harvest_spread < 0.05, (
        "this test only proves anything if the harvest numbers DO collapse; "
        f"spread was {harvest_spread:.3f}s — retune the timings above"
    )
    for src, value in harvest_elapsed.items():
        assert value >= COLLECTOR_STALL_S, (
            f"harvest time for {src} ({value:.3f}s) should be dominated by the "
            "collector stall, not by the task"
        )

    assert set(true_elapsed) == {"fast", "medium", "slow"}
    for src, expected in (("fast", FAST_S), ("medium", MEDIUM_S), ("slow", SLOW_S)):
        assert true_elapsed[src] == pytest.approx(expected, abs=TOLERANCE_S), (
            f"{src} reported {true_elapsed[src]:.3f}s for {expected:.3f}s of work"
        )
        assert true_elapsed[src] < COLLECTOR_STALL_S, (
            f"{src} reported {true_elapsed[src]:.3f}s — that is the collector's "
            "stall leaking into the source's own duration, i.e. the original bug"
        )

    assert true_elapsed["fast"] < true_elapsed["medium"] < true_elapsed["slow"], (
        "true_elapsed must ORDER the sources; ordering is what tells you which "
        f"one to chase: {true_elapsed}"
    )
    assert set(started) == {"fast", "medium", "slow"}


def test_each_task_gets_its_own_start_time():
    """Not one shared ``t_gather_start``. Sources created at different moments
    reported mutually contaminated numbers before the fix."""
    async def scenario():
        started: dict[str, float] = {}
        true_elapsed: dict[str, float] = {}

        early = track_task_timing(
            "early", asyncio.ensure_future(_sleep_then("early", FAST_S)),
            started, true_elapsed)
        await asyncio.sleep(MEDIUM_S)
        late = track_task_timing(
            "late", asyncio.ensure_future(_sleep_then("late", FAST_S)),
            started, true_elapsed)

        await asyncio.gather(early, late)
        return started, true_elapsed

    started, true_elapsed = asyncio.run(scenario())

    assert started["late"] - started["early"] == pytest.approx(
        MEDIUM_S, abs=TOLERANCE_S), (
        f"the two tasks were created ~{MEDIUM_S:.2f}s apart but their recorded "
        f"start times differ by {started['late'] - started['early']:.3f}s"
    )
    assert true_elapsed["late"] == pytest.approx(
        true_elapsed["early"], abs=TOLERANCE_S), (
        f"equal work reported unequal durations: {true_elapsed}"
    )


def test_a_cancelled_task_is_still_timed():
    """The curation deadline cancels outstanding tasks. How long a cancelled
    source ran is exactly the number you want when the deadline fires, so the
    callback must not be limited to tasks that returned a result."""
    async def scenario():
        started: dict[str, float] = {}
        true_elapsed: dict[str, float] = {}

        task = track_task_timing(
            "stuck", asyncio.ensure_future(_sleep_then("stuck", 30.0)),
            started, true_elapsed)
        await asyncio.sleep(MEDIUM_S)
        task.cancel()
        with pytest.raises(asyncio.CancelledError):
            await task
        await asyncio.sleep(0)
        return true_elapsed

    true_elapsed = asyncio.run(scenario())

    assert "stuck" in true_elapsed, "a cancelled task reported no duration"
    assert true_elapsed["stuck"] == pytest.approx(MEDIUM_S, abs=TOLERANCE_S)


def test_a_failed_task_is_still_timed():
    """A source that raises still tells you how long it took to fail."""
    async def _boom():
        await asyncio.sleep(FAST_S)
        raise RuntimeError("upstream is down")

    async def scenario():
        started: dict[str, float] = {}
        true_elapsed: dict[str, float] = {}
        task = track_task_timing(
            "boom", asyncio.ensure_future(_boom()), started, true_elapsed)
        with pytest.raises(RuntimeError):
            await task
        await asyncio.sleep(0)
        return true_elapsed

    true_elapsed = asyncio.run(scenario())

    assert true_elapsed["boom"] == pytest.approx(FAST_S, abs=TOLERANCE_S)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
