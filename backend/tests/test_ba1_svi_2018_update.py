"""BA1 must follow the SVI's updated definition, not the bare 2015 one.

Ghosh et al., "Updated recommendation for the benign stand-alone ACMG/AMP
criterion", Genet Med 2018 (PMID 30311383) refines BA1 to:

    "Allele frequency is >0.05 in any general continental population dataset of
     at least 2,000 observed alleles and found in a gene without a gene- or
     variant-specific BA1 modification."

Two requirements the engine did not implement:

  1. THE >=2,000 OBSERVED ALLELE MINIMUM. A high AF computed from a handful of
     alleles is noise, and gnomAD AN genuinely does fall to the hundreds at
     poorly covered sites. BA1 forces a Benign classification and zeroes every
     pathogenic criterion, so firing it off a low-AN site is the worst-case
     false benign.

  2. THE NINE-VARIANT EXCEPTION LIST. ClinGen identified nine variants with
     population MAF >5% for which BA1 must NOT be applied because there is
     evidence of pathogenicity (four Pathogenic, five VUS).

Found in practice: HFE c.845G>A (p.Cys282Tyr) — hereditary haemochromatosis,
on ClinGen's exception list and classified PATHOGENIC — has a real gnomAD v4
exome FAF95 popmax of 0.071, so the engine fired BA1 and returned **Benign**.
HFE and ACAD9 are both on this project's CVD gene panel, so two of the nine are
reachable in production.
"""
from __future__ import annotations

import pytest

from backend.acmg.hard_coded import (
    apply_cross_criterion_exclusions,
    compute_hard_coded_criteria,
)
from backend.acmg.tiers import classification_for, compute_points_total


def _ev(faf95: float, an: int = 1_000_000, hgvsc: str | None = None):
    ev = {
        "gnomad": {"ok": True, "variant_found": True, "variant": {
            "exome": {"af": faf95, "an": an, "faf95": {"popmax": faf95}},
        }},
        "vep": {"ok": True, "most_severe_consequence": "missense_variant"},
    }
    if hgvsc:
        ev["vep"]["hgvsc"] = hgvsc
    return ev


def _ba1(ev, gene, inh="AD"):
    out = compute_hard_coded_criteria(ev, {"inheritance_input": inh}, gene=gene)
    return {c["code"]: c for c in out}["BA1"]


def _tier(ev, gene, inh="AD"):
    out = compute_hard_coded_criteria(ev, {"inheritance_input": inh}, gene=gene)
    out, forced = apply_cross_criterion_exclusions(out, gene, inh)
    return forced or classification_for(compute_points_total(out))


@pytest.mark.parametrize("an,should_fire", [
    (1_000_000, True),
    (2_000, True),
    (1_999, False),
    (400, False),
    (0, False),
])
def test_ba1_requires_at_least_2000_observed_alleles(an, should_fire):
    entry = _ba1(_ev(0.30, an=an), "ZZZZ9")
    assert (entry["status"] == "met") is should_fire, (
        f"AN={an}: BA1 status {entry['status']!r} — {entry.get('evidence')!r}"
    )


def test_low_an_reason_names_the_allele_floor():
    ev_text = _ba1(_ev(0.30, an=500), "ZZZZ9")["evidence"]
    assert "500" in ev_text and "2,000" in ev_text.replace("2000", "2,000")


def test_low_an_does_not_force_a_benign_classification():
    """The danger is BA1 zeroing pathogenic evidence off a noisy AF."""
    assert _tier(_ev(0.30, an=500), "ZZZZ9") != "Benign"


def test_missing_an_is_not_treated_as_meeting_the_floor():
    """AN absent → cannot confirm >=2,000 → fail closed."""
    ev = {
        "gnomad": {"ok": True, "variant_found": True, "variant": {
            "exome": {"af": 0.3, "faf95": {"popmax": 0.3}},
        }},
        "vep": {"ok": True, "most_severe_consequence": "missense_variant"},
    }
    assert _ba1(ev, "ZZZZ9")["status"] == "not_met"


EXCEPTIONS = [
    ("ACAD9", "NM_014049.5:c.-44_-41dupTAAG"),
    ("GJB2", "NM_004004.6:c.109G>A"),
    ("HFE", "NM_000410.4:c.187C>G"),
    ("HFE", "NM_000410.4:c.845G>A"),
    ("MEFV", "NM_000243.3:c.1105C>T"),
    ("MEFV", "NM_000243.3:c.1223G>A"),
    ("PIBF1", "NM_006346.3:c.1214G>A"),
    ("ACADS", "NM_000017.4:c.511C>T"),
    ("BTD", "NM_000060.4:c.1330G>C"),
]


@pytest.mark.parametrize("gene,hgvsc", EXCEPTIONS)
def test_ba1_withheld_for_the_nine_clingen_exception_variants(gene, hgvsc):
    entry = _ba1(_ev(0.30, hgvsc=hgvsc), gene)
    assert entry["status"] == "not_met", (
        f"BA1 fired on ClinGen exception {gene} {hgvsc}: "
        f"{entry.get('evidence')!r}"
    )


@pytest.mark.parametrize("gene,hgvsc", EXCEPTIONS)
def test_exception_variants_never_forced_benign(gene, hgvsc):
    assert _tier(_ev(0.30, hgvsc=hgvsc), gene) != "Benign"


def test_hfe_c282y_real_world_frequency_no_longer_reads_benign():
    """The case that surfaced this: real gnomAD v4 numbers for HFE C282Y."""
    ev = _ev(0.07099497, an=1_461_852, hgvsc="NM_000410.4:c.845G>A")
    ev["gnomad"]["variant"]["exome"]["ac_hom"] = 3248
    assert _ba1(ev, "HFE")["status"] == "not_met"
    assert _tier(ev, "HFE", "AR") != "Benign"


def test_exception_reason_explains_why():
    ev_text = _ba1(_ev(0.30, hgvsc="NM_000410.4:c.845G>A"), "HFE")["evidence"].lower()
    assert "exception" in ev_text or "excluded" in ev_text
    assert "clingen" in ev_text or "svi" in ev_text


def test_transcript_accession_prefix_is_ignored_when_matching():
    """Any transcript prefix (or none) must still match the published c. change."""
    for hgvsc in ("c.845G>A", "NM_000410.3:c.845G>A", "ENST00000357618.8:c.845G>A"):
        assert _ba1(_ev(0.30, hgvsc=hgvsc), "HFE")["status"] == "not_met", hgvsc


def test_match_is_case_insensitive():
    assert _ba1(_ev(0.30, hgvsc="NM_000410.4:C.845G>A"), "hfe")["status"] == "not_met"


def test_other_variants_in_an_exception_gene_still_earn_ba1():
    """Only the listed variants are excluded — not the whole gene."""
    entry = _ba1(_ev(0.30, hgvsc="NM_000410.4:c.500A>G"), "HFE")
    assert entry["status"] == "met", (
        "the exception list suppressed BA1 for an unlisted HFE variant"
    )


def test_same_hgvs_in_a_different_gene_still_earns_ba1():
    entry = _ba1(_ev(0.30, hgvsc="NM_999999.1:c.845G>A"), "ZZZZ9")
    assert entry["status"] == "met"


def _bs1(ev, gene, inh="AD"):
    out = compute_hard_coded_criteria(ev, {"inheritance_input": inh}, gene=gene)
    return {c["code"]: c for c in out}["BS1"]


def test_bs1_says_double_count_when_ba1_actually_fired():
    entry = _bs1(_ev(0.30, hgvsc="NM_1:c.9A>G"), "ZZZZ9")
    assert entry["status"] == "not_met"
    assert "fired" in entry["evidence"] or "double-count" in entry["evidence"]


@pytest.mark.parametrize("ev,gene", [
    (_ev(0.30, hgvsc="NM_000410.4:c.845G>A"), "HFE"),
    (_ev(0.30, an=480, hgvsc="NM_1:c.9A>G"), "ZZZZ9"),
])
def test_bs1_does_not_claim_ba1_coverage_when_ba1_was_withheld(ev, gene):
    """The old text read "covered by BA1" even when BA1 had declined."""
    entry = _bs1(ev, gene)
    assert entry["status"] == "not_met"
    text = entry["evidence"]
    assert "covered by BA1" not in text, f"stale claim survives: {text!r}"
    assert "WITHHELD" in text or "withheld" in text
    assert "frequency" in text.lower()


@pytest.mark.parametrize("ev,gene", [
    (_ev(0.30, hgvsc="NM_000410.4:c.845G>A"), "HFE"),
    (_ev(0.30, an=480, hgvsc="NM_1:c.9A>G"), "ZZZZ9"),
])
def test_no_frequency_benign_criterion_applies_when_ba1_is_withheld(ev, gene):
    """BS1 must not substitute for the criterion just declined."""
    out = compute_hard_coded_criteria(ev, {"inheritance_input": "AD"}, gene=gene)
    by = {c["code"]: c for c in out}
    assert by["BA1"]["status"] == "not_met"
    assert by["BS1"]["status"] == "not_met"


def test_bs1_still_fires_normally_below_the_ba1_ceiling():
    """The ordinary BS1 window is untouched by the ba1_met plumbing."""
    entry = _bs1(_ev(0.01, hgvsc="NM_1:c.9A>G"), "ZZZZ9")
    assert entry["status"] == "met"
    assert entry["criteria_strength"] == "BS1_Strong"


def test_ordinary_common_variant_unaffected():
    assert _ba1(_ev(0.30, hgvsc="NM_000257.4:c.1000A>G"), "ZZZZ9")["status"] == "met"
