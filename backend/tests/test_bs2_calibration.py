"""BS2 count-calibration and the removal of the curator-count arm.

ACMG/AMP 2015: BS2 = "Observed in a healthy adult individual for a recessive
(HOMOZYGOUS), dominant (HETEROZYGOUS), or X-linked (HEMIZYGOUS) disorder, with
full penetrance expected at an early age." Base strength Strong (-4).

Two things were wrong with the gnomAD arm:

1. NO STRENGTH LADDER. It fired full BS2_Strong (-4) at `gnomad_hom > 5` and
   nothing at all below that. Modern VCEP practice tiers BS2 by observation
   count. The ClinGen Pulmonary Arterial Hypertension / BMPR2 VCEP — an
   autosomal-dominant, incomplete-penetrance disease, the closest published
   analogue to the cardiac genes here — specifies BS2 as
   "≥ 3 counts: BS2_strong; ≥ 2 counts: BS2_supp"
   (Eichstaedt et al., Hum Mutat 2025, Table 1; CSpec BMPR2). So a variant seen
   homozygous in 3 healthy adults earned nothing under the old >5 rule, and one
   seen in 6 jumped straight to Strong with no intermediate tier.

2. IT GATED ON THE PROBAND'S ZYGOSITY. `zygosity == "hom"` was required before
   the gnomAD homozygote count was even consulted. BS2 is a statement about the
   healthy CONTROL's genotype, not the proband's — a heterozygous proband whose
   variant is homozygous in several healthy adults is exactly the observation
   BS2 describes, and it was silently dropped.

Heterozygous gnomAD counts remain unused: BMPR2 is explicit that BS2 "cannot be
used for heterozygotes due to incomplete penetrance", which is also why the
Cardiomyopathy / Long-QT / FBN1 VCEPs mark BS2 Not Applicable outright.

The curator-count arm (`unaffected_obs_count`) is GONE — its input box was
removed from the UI in the 9→3 segregation trim, so the field was dead.
"""
from __future__ import annotations

import pytest

from backend.acmg.hard_coded import compute_hard_coded_criteria


def _ev(hom=0, hemi=0):
    return {"gnomad": {"ok": True, "variant_found": True, "variant": {
        "exome": {"ac_hom": hom, "ac_hemi": hemi},
    }}}


def _bs2(ev, **cc):
    base = {"inheritance_input": "AD", "zygosity": "het"}
    base.update(cc)
    out = compute_hard_coded_criteria(ev, base, gene="ZZZZ9")
    return {c["code"]: c for c in out}["BS2"]


@pytest.mark.parametrize("hom,expected", [
    (0, None),
    (1, None),
    (2, None),
    (3, "BS2_Strong"),
    (4, "BS2_Strong"),
    (12, "BS2_Strong"),
])
def test_bs2_homozygote_count_ladder(hom, expected):
    entry = _bs2(_ev(hom=hom))
    if expected is None:
        assert entry["status"] == "not_met", (
            f"{hom} healthy homozygote(s) should not reach BS2"
        )
    else:
        assert entry["status"] == "met", f"{hom} homozygotes did not fire BS2"
        assert entry["criteria_strength"] == expected


@pytest.mark.parametrize("hemi,expected", [
    (1, None), (2, None), (3, "BS2_Strong"),
])
def test_bs2_hemizygote_count_ladder(hemi, expected):
    """X-linked disease: the healthy hemizygous male is the BS2 observation."""
    entry = _bs2(_ev(hemi=hemi), inheritance_input="XLR")
    if expected is None:
        assert entry["status"] == "not_met"
    else:
        assert entry["criteria_strength"] == expected


def test_bs2_takes_the_stronger_of_hom_and_hemi():
    entry = _bs2(_ev(hom=3, hemi=2))
    assert entry["criteria_strength"] == "BS2_Strong"


def test_bs2_never_fires_below_three_observations():
    """Deliberately MORE conservative than the published BMPR2 ladder.

    BMPR2 specifies "＞= 2 counts: BS2_supp". We withhold at 2 instead: BS2 is a
    benign criterion, a false-benign closes the diagnostic question, and the
    2-count tier is the weakest, newest place BS2 fires. Deviating from BMPR2
    only in the conservative direction is a deliberate project decision — the
    >=3 Strong bar and the removal of the proband-zygosity gate still stand.
    """
    for n in (0, 1, 2):
        assert _bs2(_ev(hom=n))["status"] == "not_met", f"{n} fired BS2"
        assert _bs2(_ev(hemi=n))["status"] == "not_met", f"{n} hemi fired BS2"


def test_no_bs2_supporting_tier_exists():
    for n in range(0, 15):
        s = _bs2(_ev(hom=n)).get("criteria_strength")
        assert s in (None, "BS2_Strong"), f"{n} homozygotes gave {s}"


@pytest.mark.parametrize("proband_zyg", ["het", "hom", "hemi", ""])
def test_bs2_does_not_depend_on_the_probands_own_zygosity(proband_zyg):
    """BS2 describes the healthy control's genotype, not the proband's."""
    entry = _bs2(_ev(hom=4), zygosity=proband_zyg)
    assert entry["status"] == "met", (
        f"proband zygosity {proband_zyg!r} suppressed a healthy-homozygote "
        f"observation: {entry.get('evidence')!r}"
    )
    assert entry["criteria_strength"] == "BS2_Strong"


def test_heterozygous_gnomad_counts_never_drive_bs2():
    """Incomplete penetrance — het controls are not BS2 evidence (BMPR2 VCEP)."""
    ev = {"gnomad": {"ok": True, "variant_found": True, "variant": {
        "exome": {"ac": 900, "ac_het": 900, "ac_hom": 0, "ac_hemi": 0},
    }}}
    assert _bs2(ev)["status"] == "not_met"


def test_unaffected_obs_count_is_no_longer_an_input():
    """The UI box is gone; a stray key must not resurrect the old -4 path."""
    entry = _bs2(_ev(), unaffected_obs_count=4)
    assert entry["status"] == "not_met", (
        "the deleted curator BS2 arm still fires from a leftover key"
    )


def test_curation_request_rejects_the_deleted_fields():
    from backend.models import CurationRequest

    req = CurationRequest(gene="MYH7", hgvs_c="c.1A>G", phenotype="HCM")
    for gone in ("unaffected_obs_count", "seg_unaffected_carriers"):
        assert not hasattr(req, gone), f"{gone} still on the request model"


def test_vcep_not_applicable_genes_still_never_fire_bs2():
    """Cardiomyopathy / FBN1 VCEPs remove BS2 for incomplete penetrance."""
    for gene in ("MYH7", "MYBPC3", "FBN1"):
        out = compute_hard_coded_criteria(
            _ev(hom=9), {"inheritance_input": "AD", "zygosity": "het"}, gene=gene,
        )
        entry = {c["code"]: c for c in out}["BS2"]
        assert entry["status"] == "not_met", f"BS2 fired for {gene}"


@pytest.mark.parametrize("hom", [2, 3, 9, 40])
def test_bs2_suppressed_for_mitochondrial_inheritance(hom):
    """mtDNA has heteroplasmy levels, not nuclear zygosity.

    gnomAD's mtDNA calls are heteroplasmy-fraction based, so a "homoplasmic"
    tally is not the healthy-homozygote count the BMPR2 ladder was calibrated
    on. Fail closed rather than reinterpret the number.
    """
    entry = _bs2(_ev(hom=hom), inheritance_input="MT")
    assert entry["status"] == "not_met", (
        f"BS2 fired on {hom} gnomAD counts for an MT variant: "
        f"{entry.get('evidence')!r}"
    )


def test_mt_suppression_explains_heteroplasmy():
    ev_text = _bs2(_ev(hom=5), inheritance_input="MT")["evidence"].lower()
    assert "heteroplasm" in ev_text, (
        f"MT reason should name heteroplasmy: {ev_text!r}"
    )


@pytest.mark.parametrize("inh", ["AD", "AR", "XLD", "XLR", ""])
def test_non_mt_inheritance_still_uses_the_ladder(inh):
    """Only MT is suppressed — every other mode keeps the calibrated ladder."""
    assert _bs2(_ev(hom=3), inheritance_input=inh)["criteria_strength"] == "BS2_Strong"


def test_missing_gnomad_does_not_fire_bs2():
    entry = _bs2({"gnomad": {"ok": False}})
    assert entry["status"] == "not_met"


def test_not_met_reason_quotes_the_threshold():
    ev_text = _bs2(_ev(hom=1))["evidence"]
    assert "1" in ev_text and "3" in ev_text, (
        f"reason should tell the curator how many are needed: {ev_text!r}"
    )
    assert "at least 3" in ev_text.lower()
    if "supporting" in ev_text.lower():
        assert "not applied" in ev_text.lower(), (
            f"reason mentions a Supporting tier without disclaiming it: {ev_text!r}"
        )


BA1_EXCEPTION_SAMPLES = [
    ("HFE", "NM_000410.4:c.845G>A"),
    ("HFE", "NM_000410.4:c.187C>G"),
    ("GJB2", "NM_004004.6:c.109G>A"),
    ("BTD", "NM_000060.4:c.1330G>C"),
    ("ACADS", "NM_000017.4:c.511C>T"),
    ("ACAD9", "NM_014049.5:c.-44_-41dupTAAG"),
    ("MEFV", "NM_000243.3:c.1105C>T"),
    ("MEFV", "NM_000243.3:c.1223G>A"),
    ("PIBF1", "NM_006346.3:c.1214G>A"),
]


@pytest.mark.parametrize("gene,hgvsc", BA1_EXCEPTION_SAMPLES)
def test_bs2_withheld_for_ba1_exception_variants(gene, hgvsc):
    """ACMG BS2 needs "full penetrance expected at an early age".

    All nine ClinGen BA1-exception variants are common REDUCED-PENETRANCE /
    hypomorphic alleles — that is why they are common yet not benign — so BS2's
    own precondition fails regardless of how many healthy homozygotes exist.
    """
    ev = {"gnomad": {"ok": True, "variant_found": True, "variant": {
              "exome": {"ac_hom": 3248, "an": 1_461_852}}},
          "vep": {"ok": True, "most_severe_consequence": "missense_variant",
                  "hgvsc": hgvsc}}
    out = compute_hard_coded_criteria(
        ev, {"inheritance_input": "AR", "zygosity": "hom"}, gene=gene,
    )
    entry = {c["code"]: c for c in out}["BS2"]
    assert entry["status"] == "not_met", (
        f"BS2 fired on BA1-exception variant {gene} {hgvsc}: "
        f"{entry.get('evidence')!r}"
    )
    assert "penetran" in entry["evidence"].lower()


def test_hfe_c282y_no_longer_reaches_likely_benign():
    """The end-to-end case: real gnomAD v4 numbers, ClinGen says Pathogenic."""
    from backend.acmg.hard_coded import apply_cross_criterion_exclusions
    from backend.acmg.tiers import classification_for, compute_points_total

    ev = {"gnomad": {"ok": True, "variant_found": True, "variant": {"exome": {
              "af": 0.05903, "an": 1_461_852, "ac_hom": 3248,
              "faf95": {"popmax": 0.07099497}}}},
          "vep": {"ok": True, "most_severe_consequence": "missense_variant",
                  "hgvsc": "NM_000410.4:c.845G>A"}}
    out = compute_hard_coded_criteria(
        ev, {"inheritance_input": "AR", "zygosity": "hom"}, gene="HFE",
    )
    by = {c["code"]: c for c in out}
    for code in ("BA1", "BS1", "BS2"):
        assert by[code]["status"] == "not_met", f"{code} still fired"
    out, forced = apply_cross_criterion_exclusions(out, "HFE", "AR")
    tier = forced or classification_for(compute_points_total(out))
    assert tier not in ("Benign", "Likely benign"), tier


def test_unlisted_variant_in_an_exception_gene_still_gets_bs2():
    """Only the nine listed variants are excluded — not the whole gene."""
    ev = {"gnomad": {"ok": True, "variant_found": True, "variant": {
              "exome": {"ac_hom": 9, "an": 1_461_852}}},
          "vep": {"ok": True, "most_severe_consequence": "missense_variant",
                  "hgvsc": "NM_000410.4:c.500A>G"}}
    out = compute_hard_coded_criteria(
        ev, {"inheritance_input": "AR", "zygosity": "hom"}, gene="HFE",
    )
    assert {c["code"]: c for c in out}["BS2"]["status"] == "met"
