"""PP1 and BS4 are opposites and must never fire together.

  PP1 = "Cosegregation with disease in multiple affected family members in a
         gene definitively known to cause the disease."
  BS4 = "Lack of segregation in affected members of a family."

One segregation analysis cannot both support and refute cosegregation. Both are
AI-evaluated, so nothing stopped the model emitting both from one family history
— which scored +PP1 and -BS4 simultaneously and showed the curator two
contradictory statements about the same pedigree.

Resolution follows the existing PP3-vs-BP4 precedent in
``apply_cross_criterion_exclusions``: two opposite-direction criteria drawn from
the SAME evidence class cancel — keep NEITHER and surface the conflict, rather
than have the engine silently pick a side on data a human needs to adjudicate.

Ordering matters and is asserted here: the rule runs AFTER the PP1
literature-facts hardening (so a PP1 the facts gate already demoted cannot
create a phantom conflict) and BEFORE the PP5/BP6 reputational guards (so a
BS4 that has just been cancelled cannot go on to suppress PP5).
"""
from __future__ import annotations

from backend.acmg.hard_coded import apply_cross_criterion_exclusions


def _c(code: str, strength: str | None, status: str = "met",
       direction: str | None = None, evidence: str = "") -> dict:
    if direction is None:
        direction = "benign" if code.startswith("B") else "pathogenic"
    return {
        "code": code, "status": status, "criteria_strength": strength,
        "direction": direction, "evidence": evidence or f"{code} evidence",
        "source": "ai",
    }


def _run(criteria, **kw):
    out, forced = apply_cross_criterion_exclusions(criteria, **kw)
    return {c["code"]: c for c in out}, forced


def test_pp1_and_bs4_cancel_each_other():
    by, _ = _run(
        [
            _c("PP1", "PP1_Supporting",
               evidence="3 affected carriers, PMID 12345678"),
            _c("BS4", "BS4_Strong",
               evidence="2 affected relatives do not carry the variant"),
        ],
        has_curator_segregation=True,
    )
    assert by["PP1"]["status"] == "not_met", "PP1 survived a direct contradiction"
    assert by["BS4"]["status"] == "not_met", "BS4 survived a direct contradiction"


def test_cancelled_pair_contributes_zero_points():
    from backend.acmg.tiers import compute_points_total

    by, _ = _run(
        [_c("PP1", "PP1_Strong"), _c("BS4", "BS4_Strong")],
        has_curator_segregation=True,
    )
    assert compute_points_total(list(by.values())) == 0


def test_conflict_is_explained_to_the_curator():
    by, _ = _run(
        [_c("PP1", "PP1_Moderate"), _c("BS4", "BS4_Strong")],
        has_curator_segregation=True,
    )
    for code, other in (("PP1", "BS4"), ("BS4", "PP1")):
        ev = by[code]["evidence"]
        assert other in ev, f"{code}'s reason does not name {other}: {ev!r}"
        assert "segregation" in ev.lower()


def test_pp1_alone_is_untouched():
    by, _ = _run(
        [_c("PP1", "PP1_Supporting", evidence="4 affected carriers")],
        has_curator_segregation=True,
    )
    assert by["PP1"]["status"] == "met"
    assert by["PP1"]["criteria_strength"] == "PP1_Supporting"


def test_bs4_alone_is_untouched():
    by, _ = _run([_c("BS4", "BS4_Strong")], has_curator_segregation=True)
    assert by["BS4"]["status"] == "met"
    assert by["BS4"]["criteria_strength"] == "BS4_Strong"


def test_not_met_pp1_does_not_cancel_a_met_bs4():
    by, _ = _run(
        [_c("PP1", None, status="not_met"), _c("BS4", "BS4_Strong")],
        has_curator_segregation=True,
    )
    assert by["BS4"]["status"] == "met", "a not_met PP1 created a phantom conflict"


def test_not_met_bs4_does_not_cancel_a_met_pp1():
    by, _ = _run(
        [_c("PP1", "PP1_Supporting"), _c("BS4", None, status="not_met")],
        has_curator_segregation=True,
    )
    assert by["PP1"]["status"] == "met"


def test_pp1_demoted_by_the_facts_gate_leaves_bs4_standing():
    """A literature PP1 with no PMID is dropped by the facts gate first.

    That is not a segregation contradiction, so BS4 must survive intact.
    """
    by, _ = _run(
        [
            _c("PP1", "PP1_Strong", evidence="segregates in the family"),
            _c("BS4", "BS4_Strong"),
        ],
        has_curator_segregation=False,
    )
    assert by["PP1"]["status"] == "not_met"
    assert by["BS4"]["status"] == "met", (
        "BS4 was cancelled by a PP1 that the facts gate had already removed"
    )


def test_cancelled_bs4_no_longer_suppresses_pp5():
    """PP5 is withheld when strong benign evidence stands — but a BS4 that has
    just been cancelled as contradictory is not standing evidence."""
    by, _ = _run(
        [
            _c("PP1", "PP1_Strong"),
            _c("BS4", "BS4_Strong"),
            _c("PP5", "PP5_Strong", direction="pathogenic"),
        ],
        has_curator_segregation=True,
    )
    assert by["BS4"]["status"] == "not_met"
    assert by["PP5"]["status"] == "met", (
        "a cancelled BS4 still suppressed PP5"
    )


def test_gate_covers_the_deterministic_no_ai_path():
    """The double-fire is reachable WITHOUT the AI, from structured counts alone.

    ``no_ai.build_no_ai_criteria`` derives PP1 from ``seg_affected_carriers`` and
    BS4 from ``seg_affected_noncarriers`` independently, so a family with both
    produced both criteria deterministically. This is the reachable path, not
    just a hypothetical model slip, and app.py runs the same exclusions over it.
    """
    from backend.acmg.hard_coded import compute_hard_coded_criteria
    from backend.acmg.no_ai import build_no_ai_criteria, infer_supplementary_criteria

    cc = {
        "inheritance_input": "AD", "zygosity": "het",
        "seg_affected_carriers": 4, "seg_affected_noncarriers": 2,
        "seg_meioses": 6,
    }
    hc = compute_hard_coded_criteria({}, cc, gene="MYH7")
    crit = build_no_ai_criteria(hc, infer_supplementary_criteria({}, cc, "MYH7", hc))
    before = {c["code"]: c["status"] for c in crit}
    assert before["PP1"] == "met" and before["BS4"] == "met", (
        "precondition changed — this test needs both to fire pre-exclusion"
    )

    out, _ = apply_cross_criterion_exclusions(
        crit, "MYH7", "AD", has_curator_segregation=True,
    )
    after = {c["code"]: c["status"] for c in out}
    assert after["PP1"] == "not_met"
    assert after["BS4"] == "not_met"


def test_standing_bs4_still_suppresses_pp5():
    """The existing guard must keep working when there is no PP1 conflict."""
    by, _ = _run(
        [
            _c("BS4", "BS4_Strong"),
            _c("PP5", "PP5_Strong", direction="pathogenic"),
        ],
        has_curator_segregation=True,
    )
    assert by["BS4"]["status"] == "met"
    assert by["PP5"]["status"] == "not_met"
