"""The header markup invariant the responsive layout depends on.

Both headers wrap to two rows below 760px: brand and account chip on row 1, nav
links on row 2 (see the HEADER: TWO ROWS block at the end of static/heartvar.css).
That layout is driven by ``order`` and ``flex-basis`` on FLEX SIBLINGS, so it only
works while each account-chip slot is a sibling of its <nav>, not a child of it.

Worth pinning because the failure is silent and only visible on a phone. The chip
and the admin-only Admin tab were originally added inside a nav tuned for two
links; on a 390px screen the app nav then ran 118px past the viewport, putting the
Admin tab half off-screen and the account chip entirely out of reach, while the
landing header instead overlapped its logo with "Variant Curator" by 42px. Nothing
in CI catches that — there is no browser in the test run — so the structural
precondition is what gets tested here.

The Admin *link* deliberately stays inside the landing nav: it is a destination,
unlike sign-in/sign-out.
"""
from __future__ import annotations

import pathlib
from html.parser import HTMLParser

import pytest

_REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
_INDEX = _REPO_ROOT / "index.html"

_SLOTS = {
    "hv-auth-slot-landing": False,
    "hv-auth-slot-app": False,
    "hv-admin-slot-landing": True,
}


class _SlotNesting(HTMLParser):
    """Records, for each interesting id, whether it appeared inside a <nav>."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.nav_depth = 0
        self.found: dict[str, bool] = {}

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attr = dict(attrs)
        el_id = attr.get("id")
        if el_id in _SLOTS:
            self.found[el_id] = self.nav_depth > 0
        if tag == "nav":
            self.nav_depth += 1

    def handle_endtag(self, tag: str) -> None:
        if tag == "nav" and self.nav_depth > 0:
            self.nav_depth -= 1


def _nesting(html: str) -> dict[str, bool]:
    parser = _SlotNesting()
    parser.feed(html)
    return parser.found


@pytest.fixture(scope="module")
def nesting() -> dict[str, bool]:
    return _nesting(_INDEX.read_text(encoding="utf-8"))


@pytest.mark.parametrize("slot_id", sorted(_SLOTS))
def test_every_header_slot_is_present(nesting, slot_id):
    assert slot_id in nesting, f"{slot_id} is missing from index.html"


@pytest.mark.parametrize(
    "slot_id", sorted(i for i, in_nav in _SLOTS.items() if not in_nav)
)
def test_account_chip_slots_are_siblings_of_the_nav(nesting, slot_id):
    """Inside the <nav> the chip cannot be re-ordered onto row 1, so the phone
    header falls back to one overflowing row."""
    assert nesting[slot_id] is False, (
        f"{slot_id} is inside a <nav>; the two-row phone header orders it "
        f"against the nav as a flex sibling, so it must sit outside"
    )


@pytest.mark.parametrize("slot_id", sorted(i for i, in_nav in _SLOTS.items() if in_nav))
def test_navigation_slots_stay_inside_the_nav(nesting, slot_id):
    assert nesting[slot_id] is True, f"{slot_id} is a destination and belongs in the <nav>"


def test_the_two_row_header_breakpoint_is_documented_in_css():
    """The 760px breakpoint is load-bearing and not obvious: the rest of the
    phone layout switches at 620px. A signed-in admin overflowed a 621px viewport
    by 47px, because the 42px serif wordmark is still shown above 620px."""
    css = (_REPO_ROOT / "static" / "heartvar.css").read_text(encoding="utf-8")
    assert "@media (max-width: 760px)" in css
    for selector in (
        "#hv-auth-slot-app{ order:1",
        "#hv-landing #hv-auth-slot-landing{ order:1",
    ):
        assert selector in css, f"missing the row-1 ordering rule: {selector}"


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
