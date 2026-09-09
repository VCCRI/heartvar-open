"""The reveal animation must not get slower as the backend gets faster.

REPORTED FROM PRODUCTION: "the loading screen still shows for a few seconds
after all databases are ticked green." It did, and the cause was a perverse
coupling rather than any backend latency.

During loading, fakeCritSim ticks one ACMG box every SIM_TICK_MS (1500 ms) and
stops with SIM_RESERVE (4) boxes left. On reveal, rapidComplete used to tick
EVERY remaining box at a fixed SIM_FAST_MS (70 ms) and then hold 500 ms. So the
leftover count — and therefore the delay — was set by how FAST the gather was:

    gather 40 s -> ~4 boxes left  ->  4 x 70 + 500  = ~0.8 s
    gather  3 s -> ~25 boxes left -> 25 x 70 + 500  = ~2.3 s

Every second cut from the gather added ~0.5 s of animation here, so the work
done to speed up the curation was partly spent back on watching boxes tick.

A fixed TOTAL budget decouples them. These are static guards over the source,
in the same style as test_acmg_constants_parity's guard on this file, because
there is no browser in the suite to time.
"""
from __future__ import annotations

import pathlib
import re

import pytest

JS = (pathlib.Path(__file__).resolve().parent.parent.parent
      / "static" / "heartvar.js")


def _src() -> str:
    return JS.read_text(encoding="utf-8")


def _const(name: str) -> int:
    m = re.search(rf"const\s+{name}\s*=\s*(\d+)", _src())
    assert m, f"{name} not found in heartvar.js"
    return int(m.group(1))


def test_the_finale_has_a_total_budget():
    assert "SIM_FINALE_BUDGET_MS" in _src(), (
        "the rapid-complete finale has no total budget, so its duration scales "
        "with how many boxes the gather left unticked — i.e. it gets slower as "
        "the backend gets faster"
    )


def test_the_per_box_pace_is_derived_from_the_budget():
    """A budget constant that nothing divides by is decoration."""
    src = _src()
    assert re.search(r"SIM_FINALE_BUDGET_MS\s*/\s*this\.remaining\.length", src), (
        "the finale does not divide its budget by the remaining box count"
    )


def test_the_worst_case_finale_stays_under_half_a_second():
    """The whole point: the reveal must feel immediate no matter how quick the
    gather was. 27 boxes is the full grid."""
    budget = _const("SIM_FINALE_BUDGET_MS")
    hold = _const("SIM_DONE_HOLD_MS")
    worst = budget + hold
    assert worst <= 500, (
        f"a fast gather would still cost {worst} ms of animation after the "
        f"result is ready (budget {budget} + hold {hold})"
    )


def test_the_reveal_hold_is_not_the_old_half_second():
    """Guards the specific regression: a literal 500 ms timeout before the
    reveal callback."""
    src = _src()
    assert "setTimeout(cb || (() => {}), 500)" not in src, (
        "the reveal still holds a hard-coded 500 ms"
    )
    assert "setTimeout(cb || (() => {}), SIM_DONE_HOLD_MS)" in src


def test_something_still_animates_on_reveal():
    """The reserve exists so the reveal is not a dead jump. Budgeting the
    finale must not have removed it."""
    assert _const("SIM_RESERVE") >= 1, (
        "no boxes are reserved, so the reveal animates nothing"
    )
    assert re.search(r"Math\.max\(\s*\d+", _src()), (
        "no floor on the per-box pace — a big leftover would jump instantly"
    )


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
