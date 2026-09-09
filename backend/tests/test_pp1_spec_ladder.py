"""PP1's co-segregation ladder, graded against each gene's PUBLISHED thresholds.

Before 2026-09-07 `seg_affected_carriers` was capped at Supporting no matter how
large it was, so 11 affected relatives scored the same +1 as 2. The
Cardiomyopathy VCEP publishes its ladder in SEGREGATIONS ("≥7 segregations
(LOD 2.1) for STRONG") and a count of affected relatives carrying the variant IS
that count, so the cap discarded the evidence. PP1 came out Supporting on every
variant that carried it, where the specs' own thresholds would have given Strong
or Moderate, and variants sat at VUS instead of LP.

⚠ FOUR EARLIER ATTEMPTS AT THIS WERE REFUTED. The two failures these tests exist
to prevent:

  1. NON-MONOTONIC SCORING. `n = meioses if meioses else carriers` meant 7
     carriers with the meiosis field blank scored Strong while the same 7
     carriers with `meioses=1` scored zero — filling an optional field
     truthfully cost 4 points. Fixed by taking the max over the axes.
  2. UNIT CONFUSION ("FBN1's spec row is read two ways in one function").
     "segregations" / "meioses" / "affected family members" exclude the proband;
     "affected individuals" includes it. FBN1 is the only gene using the
     proband-inclusive phrasing, so its 2-3/4/≥5 is normalised to 1/3/4.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from backend.acmg.hard_coded import (
    _PP1_DEFAULT_LADDER,
    _pp1_spec_ladder,
    _pp1_strength_from_lod_meioses as ladder,
)

RANK = {None: 0, "Supporting": 1, "Moderate": 2, "Strong": 3}
SPEC = json.loads(
    (Path(__file__).resolve().parents[1] / "data" / "vcep_criteria_spec.json").read_text()
)


def test_the_published_ladders_are_exactly_two_groups():
    """🔴 PIN. 26 genes publish 3/5/7 and FBN1 alone publishes 2/4/5 (which is
    1/3/4 once the proband is removed). If a re-harvest changes this, the fix's
    whole premise needs re-reading rather than the numbers quietly moving."""
    groups: dict[tuple[int, int, int], list[str]] = {}
    for gene in SPEC["genes"]:
        groups.setdefault(_pp1_spec_ladder(gene), []).append(gene)
    assert {k: len(v) for k, v in groups.items()} == {(3, 5, 7): 26, (1, 3, 4): 1}
    assert groups[(1, 3, 4)] == ["FBN1"]


def test_fbn1_is_normalised_out_of_the_proband_inclusive_unit():
    """GN022 says "2-3 affected individuals" for Supporting. Two individuals is
    the proband plus ONE relative, i.e. one co-segregation — so the threshold in
    relatives is 1, not 2. Reading the 2 as relatives would make FBN1 look a
    tier stricter than it is."""
    assert _pp1_spec_ladder("FBN1") == (1, 3, 4)


def test_kcnq1_reads_its_family_member_count_not_its_schwartz_score():
    """GN112's row contains "Schwartz score ... >3" BEFORE "the proband + 3
    affected family members". A bare-number pattern reads the Schwartz
    threshold as a segregation count and happens to get 3 by luck at
    Supporting — and then the wrong number at Moderate and Strong."""
    assert _pp1_spec_ladder("KCNQ1") == (3, 5, 7)


def test_an_unspecd_gene_falls_back_to_the_svi_default(  ):
    assert _pp1_spec_ladder("TTN") == _PP1_DEFAULT_LADDER == (3, 5, 7)
    assert _pp1_spec_ladder(None) == _PP1_DEFAULT_LADDER


@pytest.mark.parametrize("carriers,expected", [
    (1, None), (2, "Supporting"), (3, "Supporting"), (4, "Supporting"),
    (5, "Moderate"), (6, "Moderate"), (7, "Strong"), (11, "Strong"),
])
def test_carriers_are_graded_on_the_ladder_not_capped(carriers, expected):
    assert ladder(None, 0, carriers, 0, gene="MYH7") == expected


def test_a_large_carrier_count_no_longer_scores_the_same_as_two():
    """The one-line statement of the bug."""
    assert ladder(None, 0, 2, 0, gene="MYH7") == "Supporting"
    assert ladder(None, 0, 11, 0, gene="MYH7") == "Strong"


def test_fbn1_grades_on_its_own_ladder():
    assert ladder(None, 0, 1, 0, gene="FBN1") == "Supporting"
    assert ladder(None, 0, 3, 0, gene="FBN1") == "Moderate"
    assert ladder(None, 0, 4, 0, gene="FBN1") == "Strong"
    assert ladder(None, 0, 4, 0, gene="MYH7") == "Supporting"


def test_the_exact_refuted_counterexample():
    """7 carriers, meiosis field blank vs truthfully set to 1. The blank form
    must not outscore the filled one."""
    blank = ladder(None, 0, 7, 0, gene="MYH7")
    filled = ladder(None, 1, 7, 0, gene="MYH7")
    assert RANK[filled] >= RANK[blank], f"blank={blank} filled={filled}"


def test_monotonic_across_the_whole_input_grid():
    """Raising either count, on any gene, must never lower the tier."""
    inversions = []
    for gene in ("MYH7", "FBN1", "PTPN11", "BMPR2", "KCNQ1", None):
        for unaff in (0, 1):
            for c in range(0, 13):
                for m in range(0, 13):
                    base = ladder(None, m, c, unaff, gene=gene)
                    for c2, m2 in ((c + 1, m), (c, m + 1)):
                        nxt = ladder(None, m2, c2, unaff, gene=gene)
                        if RANK[nxt] < RANK[base]:
                            inversions.append((gene, c, m, unaff, base, nxt))
    assert not inversions, f"{len(inversions)} inversions, first: {inversions[:3]}"


def test_still_requires_one_genotyped_affected_carrier():
    """Affected-but-untested relatives are PP4 context, not PP1."""
    assert ladder(None, 10, 0, 0, gene="MYH7") is None
    assert ladder(2.5, 10, 0, 0, gene="MYH7") is None


def test_the_two_carrier_supporting_floor_survives():
    """no_ai.py has always emitted Supporting at 2 carriers, below the 3/5/7
    ladder's own first rung. Dropping it would silently regress low-count rows."""
    assert ladder(None, 0, 2, 0, gene="MYH7") == "Supporting"


def test_explicit_lod_still_wins_where_it_is_the_strongest_axis():
    assert ladder(2.1, 0, 1, 0, gene="MYH7") == "Strong"
    assert ladder(1.5, 0, 1, 0, gene="MYH7") == "Moderate"
    assert ladder(0.9, 0, 1, 0, gene="MYH7") == "Supporting"
    assert ladder(0.5, 0, 1, 0, gene="MYH7") is None


def test_the_unaffected_carrier_temper_still_lowers_the_tier():
    """Contrary evidence on a DIFFERENT axis, so lowering is correct here and is
    not the monotonicity defect."""
    assert ladder(None, 0, 7, 0, gene="MYH7") == "Strong"
    assert ladder(None, 0, 7, 1, gene="MYH7") == "Moderate"


def test_meioses_alone_still_grades_when_no_carrier_count_is_richer():
    assert ladder(None, 7, 1, 0, gene="MYH7") == "Strong"
    assert ladder(None, 5, 1, 0, gene="MYH7") == "Moderate"
    assert ladder(None, 3, 1, 0, gene="MYH7") == "Supporting"


def test_junk_inputs_do_not_raise():
    for bad in ("", "x", None, True, -1, 3.7):
        ladder(bad, bad, bad, bad, gene="MYH7")
        ladder(None, bad, 5, bad, gene=bad if isinstance(bad, str) else None)
