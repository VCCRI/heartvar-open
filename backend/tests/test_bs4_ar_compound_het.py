"""BS4 does not apply to AR compound heterozygosity.

Biesecker et al. 2024 (ClinGen guidance for PP1/BS4 and PP4, AJHG 111:24-38),
verbatim: "it is critical to recognize that in autosomal-recessive inheritance
with compound heterozygosity in the affected proband, the existence of
non-segregations in relatives at a given locus provides little to no evidence
that those variants are not pathogenic. This is because only one of the two
alleles is necessarily benign, and the non-segregation data cannot distinguish in
a single family which of the variants is benign."

BS4 is -4.0 points, so leaving it in place here is a four-point false-benign
argument resting on evidence the guidance says carries none.

SCOPE MATTERS. The same paragraph KEEPS BS4 for autosomal-recessive inheritance
with homozygosity, so the exclusion must key on compound heterozygosity
specifically, not on AR inheritance.
"""
from __future__ import annotations

from backend.acmg.hard_coded import apply_cross_criterion_exclusions


def _bs4_met():
    return [{"code": "BS4", "status": "met", "criteria_strength": "BS4",
             "evidence": "two affected relatives do not carry the variant"}]


def _run(**kw):
    criteria, _forced = apply_cross_criterion_exclusions(_bs4_met(), **kw)
    return {c["code"]: c for c in criteria}["BS4"]


def test_bs4_withdrawn_for_ar_compound_heterozygote():
    out = _run(inheritance="AR", zygosity="het", in_trans_pathogenic="yes")
    assert out["status"] == "not_met"
    assert out["criteria_strength"] is None
    assert "compound heterozygote" in out["evidence"]


def test_bs4_retained_for_ar_homozygote():
    """Explicitly preserved by the same paragraph."""
    out = _run(inheritance="AR", zygosity="hom", in_trans_pathogenic="yes")
    assert out["status"] == "met"


def test_bs4_retained_for_autosomal_dominant():
    out = _run(inheritance="AD", zygosity="het", in_trans_pathogenic="yes")
    assert out["status"] == "met"


def test_bs4_fails_open_when_zygosity_is_unknown():
    """With no zygosity there is no way to tell a compound heterozygote from a
    homozygote. Suppressing would silently delete valid BS4 on every AR curation
    that left the field blank, so the guard must fail OPEN."""
    for kw in (
        {"inheritance": "AR"},
        {"inheritance": "AR", "zygosity": ""},
        {"inheritance": "AR", "zygosity": "het"},
        {"inheritance": "AR", "in_trans_pathogenic": "yes"},
    ):
        assert _run(**kw)["status"] == "met", kw


def test_bs4_survives_with_no_clinical_context_at_all():
    """The default-argument path — nothing supplied, nothing suppressed."""
    assert _run()["status"] == "met"
