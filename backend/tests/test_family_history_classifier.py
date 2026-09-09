"""classify_family_history — the Summary panel's Family-history bucket.

WHY THIS FILE EXISTS. The classifier's phrase patterns fixed the WORD ORDER, so
"affected sister" matched and "Mother and sister affected" did not: it returned
"unknown". A curator who types a real family history and sees "unknown" fairly
concludes the entry was ignored. Worse, the textbook PP1 phrasing "both affected
... both tested positive for this variant; 4 informative meioses" also returned
"unknown", because the segregation pattern required an "in" directly after
"tested positive".

Widened to pin BOTH directions,
because the risk in widening the positive side is that a NEGATED statement
starts reading as positive. Negation is checked first and was widened in step,
which is what the `no_family_history` cases here guard.

DISPLAY ONLY. This bucket label reaches no criterion, strength, point total or
tier — verified in production on 2026-09-02 with two otherwise-identical
`ai_mode=none` curations that emitted byte-identical criteria sets. The raw text
does reach the AI prompt, where prompt.py admits it as a PP1/BS4 fallback.
"""
from __future__ import annotations

import pytest

from backend.acmg.hard_coded import classify_family_history


@pytest.mark.parametrize(("text", "expected"), [
    ("", "unknown"),
    ("   ", "unknown"),
    ("Adopted, family history not available", "unknown"),

    ("No family history of cardiomyopathy", "no_family_history"),
    ("no sister affected", "no_family_history"),
    ("Parents unaffected; no family history of CHD", "no_family_history"),
    ("Both parents unaffected", "no_family_history"),
    ("Apparently sporadic", "no_family_history"),
    ("Family history is negative", "no_family_history"),
    ("Family history unremarkable", "no_family_history"),
    ("Nil", "no_family_history"),
    ("Proband is the only affected individual; parents unaffected",
     "no_family_history"),
    ("Denies any family history", "no_family_history"),

    ("Father carries the same variant", "segregation_data"),
    ("Cosegregates with disease", "segregation_data"),
    ("Two affected cousins genotyped and both carry this variant",
     "segregation_data"),
    ("LOD score of 2.3 across the pedigree", "segregation_data"),
    ("Father and paternal aunt both affected with HCM and both tested "
     "positive for this variant; 4 informative meioses.", "segregation_data"),
    ("Variant was also detected in the affected brother", "segregation_data"),

    ("Affected sister", "positive_family_history"),
    ("Positive family history", "positive_family_history"),
    ("Mother and sister affected", "positive_family_history"),
    ("Brother diagnosed at 30", "positive_family_history"),
    ("Maternal grandfather died suddenly aged 40", "positive_family_history"),
])
def test_bucket(text: str, expected: str) -> None:
    assert classify_family_history(text) == expected


def test_none_is_unknown_not_a_crash():
    assert classify_family_history(None) == "unknown"


def test_negation_is_checked_before_the_cooccurrence_rule():
    """The order-independent rule is the loosest thing in the classifier, so it
    must never see a text the negative patterns would have claimed. This is the
    specific regression that widening the positive side risks."""
    for text in ("no affected relatives",
                 "no sister affected",
                 "mother and sister unaffected",
                 "parents are unaffected"):
        assert classify_family_history(text) == "no_family_history", text
