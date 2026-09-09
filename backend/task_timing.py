"""True per-task completion timing for the DB-gather, and why it is separate.

WHY THIS EXISTS. On 2026-08-30 the per-source timing log in ``evidence.py``
said fourteen sources — on three continents, over four different transports —
each took 19.37 s, agreeing to within 10 milliseconds. That is not a network
measurement. The line was::

    task_elapsed[src] = perf_counter() - task_started.get(src, t_gather_start)

and it ran when the orchestrator HARVESTED a finished task out of
``asyncio.wait(..., FIRST_COMPLETED)``, not when the task completed.
``asyncio.wait`` returns a SET, so every task harvested in one pass got the
same ``perf_counter()`` reading. Worse, every first-batch task shared one
``task_started`` value, and four sources were not in that list at all and fell
through to the ``.get()`` default. The number was orchestrator harvest latency
wearing a source's name.

THREE DEPLOYS were aimed at that number before anyone checked what it measured
(an Azure DNS server setting — reverted, made DNS worse; a retry deadline; a
connect-timeout cap). None moved the cost, and none could have: there was not
one ``attempt N/4 failed`` line in any log, so two of the three were tuning a
retry mechanism that never ran.

WHAT THIS DOES. ``track_task_timing`` stamps each task's own creation instant
and attaches a done callback, which the event loop runs the moment that task
completes. The result is the source's real duration, independent of when the
orchestrator gets round to collecting it.

BOTH NUMBERS ARE KEPT, and the gap between them is the whole diagnostic:

  * true small, harvest large  -> no source was ever slow; the delay is in the
    loop that collects them (or in whatever is back-pressuring it).
  * true also large            -> the work is real, and the question becomes
    which parts of it are avoidable.

This module deliberately imports nothing but the standard library, so the
instrument can be unit-tested without standing up the app — the failure it was
written for was trusting an instrument that had never been validated.
"""
from __future__ import annotations

import asyncio
from time import perf_counter


def track_task_timing(
    name: str,
    task: asyncio.Task,
    started: dict[str, float],
    true_elapsed: dict[str, float],
) -> asyncio.Task:
    """Record ``name``'s creation time in ``started`` and its TRUE duration in
    ``true_elapsed`` when the task completes. Returns the task, so it can wrap
    a ``create_task`` call in place.

    Call this at the moment of creation. Calling it later is only safe if no
    ``await`` has intervened, because ``t0`` is read here and not when the
    coroutine was scheduled.

    Cancelled and failed tasks are timed too — a done callback fires for every
    terminal state — which is what you want: a task killed by the curation
    deadline should still report how long it ran.
    """
    t0 = perf_counter()
    started[name] = t0
    task.add_done_callback(
        lambda _t, _s=name, _t0=t0: true_elapsed.__setitem__(_s, perf_counter() - _t0)
    )
    return task
