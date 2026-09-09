"""Cardiomyopathy VCEP PM1 phenotype caveat.

THE BUG. The CM VCEP derived its PM1 hotspot codon ranges from **HCM case
cohorts** and says so in the spec text, verbatim, for MYH7 (GN002), MYBPC3
(GN095), TNNI3 (GN098) and TNNT2 (GN099) — and for no other gene:

    "Data from HCM case cohorts was used to derive these cluster regions.
     Therefore, this rule should NOT be applied when additional evidence for
     the variant supports that the variant causes a phenotype other than HCM
     (e.g., variant seen in multiple DCM cases)."

Neither `_gate_criteria_applicability` nor `apply_cross_criterion_exclusions`
received a phenotype, so the caveat was simply unimplemented — including on the
ASSERT side, where HeartVar manufactures PM1_Moderate (+2) from an
HCM-derived range for a proband whose phenotype is DCM. `req.hpo` was already in
scope at the call site and used two lines later by `carrier_status`.

SCOPE, stated rather than hidden. The spec's condition is evidence about the
VARIANT; the strongest signal HeartVar holds is the PROBAND's own phenotype.
Matching is exact on the four curated IDs in panelapp_hpo_map.json
(hcm = HP:0001639/HP:0005157, dcm = HP:0001644/HP:0006670), because the
available cardiac_hpo_descendants.json is an inverted term->category index whose
categories are not HCM/DCM-disjoint. Suppress-only: the failure mode is a
missed +2, never an over-call.
"""
from __future__ import annotations

from backend.acmg.hard_coded import (
    _gate_criteria_applicability,
    _pm1_non_hcm_phenotype,
)

_HCM = "HP:0001639"
_DCM = "HP:0001644"
_GENERIC_CM = "HP:0001638"

_IN_RANGE_VEP = {"vep": {"ok": True, "most_severe_consequence": "missense_variant",
                         "hgvsp": "NP_000248.2:p.Arg403Gln"}}


def _pm1_met():
    return [{"code": "PM1", "status": "met", "criteria_strength": "PM1_Moderate",
             "evidence": "residue 403 sits in the MYH7 head-domain hotspot"}]


def test_helper_fires_only_on_an_unambiguous_non_hcm_phenotype():
    assert _pm1_non_hcm_phenotype([_DCM])
    assert _pm1_non_hcm_phenotype(["hp:0006670"])
    assert _pm1_non_hcm_phenotype([_HCM]) is None
    assert _pm1_non_hcm_phenotype([_HCM, _DCM]) is None
    assert _pm1_non_hcm_phenotype([_GENERIC_CM]) is None
    assert _pm1_non_hcm_phenotype([]) is None
    assert _pm1_non_hcm_phenotype(None) is None


def test_pm1_suppressed_on_a_dcm_proband_in_the_four_caveat_genes():
    """The test that would have caught the bug. Pre-fix PM1 stayed met on a
    DCM proband in every one of these genes."""
    for gene in ("MYH7", "MYBPC3", "TNNI3", "TNNT2"):
        out = _gate_criteria_applicability(
            _pm1_met(), gene, _IN_RANGE_VEP, [_DCM])[0]
        assert out["status"] == "not_met", f"{gene}: PM1 must be suppressed"
        assert out["criteria_strength"] is None
        assert "phenotype other than HCM" in out["evidence"]


def test_pm1_kept_for_an_hcm_proband():
    out = _gate_criteria_applicability(
        _pm1_met(), "MYH7", _IN_RANGE_VEP, [_HCM])[0]
    assert out["status"] == "met"


def test_pm1_kept_when_no_phenotype_is_supplied():
    """Fail-open. No HPO terms is not the same fact as "the phenotype is not
    HCM", and treating it that way would delete PM1 on every curation that
    omits phenotype."""
    for hpo in (None, [], [_GENERIC_CM]):
        out = _gate_criteria_applicability(
            _pm1_met(), "MYH7", _IN_RANGE_VEP, hpo)[0]
        assert out["status"] == "met", f"hpo={hpo}"


def test_a_string_of_hpo_ids_works_exactly_like_a_list():
    assert _pm1_non_hcm_phenotype("HP:0001644")
    assert _pm1_non_hcm_phenotype(["HP:0001644"])
    assert _pm1_non_hcm_phenotype("HP:0001644,HP:0006670")
    assert _pm1_non_hcm_phenotype("HP:0001644; HP:0006670")


def test_free_text_phenotype_is_matched_too():
    """Curators type the disease, and the field has always allowed it."""
    assert _pm1_non_hcm_phenotype("dilated cardiomyopathy")
    assert _pm1_non_hcm_phenotype("DCM")
    assert _pm1_non_hcm_phenotype("HP:0001644, dilated cardiomyopathy")


def test_free_text_hcm_blocks_the_suppression():
    for text in ("hypertrophic cardiomyopathy", "HCM", "HP:0001639",
                 "HP:0001639, hypertrophic cardiomyopathy"):
        assert _pm1_non_hcm_phenotype(text) is None, text


def test_both_phenotypes_in_free_text_is_not_evidence_either():
    """Same fail-open rule as the HPO-ID path: a proband carrying both is not
    evidence for a non-HCM phenotype."""
    assert _pm1_non_hcm_phenotype("dilated and hypertrophic cardiomyopathy") is None
    assert _pm1_non_hcm_phenotype("DCM, HCM") is None


def test_an_unrelated_phenotype_suppresses_nothing():
    for text in ("RASopathy", "Noonan syndrome", "Marfan syndrome",
                 "long QT syndrome 1", "", None):
        assert _pm1_non_hcm_phenotype(text) is None, text


def test_the_gate_end_to_end_with_a_free_text_phenotype():
    """The whole point: a DCM proband loses PM1 on a caveat gene when the
    phenotype arrives the way the request actually carries it."""
    out = _gate_criteria_applicability(
        _pm1_met(), "MYH7", _IN_RANGE_VEP, "HP:0001644, dilated cardiomyopathy")[0]
    assert out["status"] == "not_met", out["evidence"]
    assert "phenotype other than HCM" in out["evidence"]


def test_the_caveat_is_keyed_to_the_spec_text_not_a_gene_list():
    """TPM1/ACTC1/MYL2/MYL3 are Cardiomyopathy-VCEP genes too, but their PM1 is
    marked Not Applicable and carries no hotspot text, so they must be handled
    by the applicability branch, not this one. KCNQ1's PM1 text carries no HCM
    caveat, so a DCM proband must NOT lose PM1 there.

    KCNQ1 needs its own inputs since its spec was implemented: PM2 met (GN112
    requires the variant to be rare before PM1 is considered) and a residue
    inside the pore helix, amino acids 300 to 320.
    """
    criteria = [
        {"code": "PM2", "status": "met", "criteria_strength": "PM2_Supporting",
         "evidence": "absent from gnomAD"},
        {"code": "PM1", "status": "met", "criteria_strength": "PM1_Moderate",
         "evidence": "in the pore helix"},
    ]
    pore_vep = {"vep": {"ok": True, "most_severe_consequence": "missense_variant",
                        "hgvsp": "NP_000209.2:p.Gly314Ser"}}
    out = [c for c in _gate_criteria_applicability(
        criteria, "KCNQ1", pore_vep, [_DCM]) if c["code"] == "PM1"][0]
    assert out["status"] == "met", out["evidence"]
    assert "phenotype other than HCM" not in out["evidence"]


def _assert_inputs():
    """PM2 met + rare missense in range is what makes the assert branch fire."""
    return [{"code": "PM2", "status": "met", "criteria_strength": "PM2_Supporting",
             "evidence": "absent from gnomAD"}]


def test_pm1_is_not_asserted_from_an_hcm_range_for_a_dcm_proband():
    """The more dangerous half: here HeartVar CREATES the +2 rather than
    passing the model's through."""
    out = _gate_criteria_applicability(
        _assert_inputs(), "MYH7", _IN_RANGE_VEP, [_DCM])
    pm1 = next((c for c in out if c["code"] == "PM1"), None)
    assert pm1 is None or pm1["status"] != "met", (
        "PM1_Moderate asserted from an HCM-derived hotspot range for a DCM "
        "proband — the exact case the CM VCEP rules out"
    )


def test_pm1_is_still_asserted_for_an_hcm_proband_in_range():
    """Guard against over-suppression: the assert branch must still work."""
    out = _gate_criteria_applicability(
        _assert_inputs(), "MYH7", _IN_RANGE_VEP, [_HCM])
    pm1 = next((c for c in out if c["code"] == "PM1"), None)
    assert pm1 is not None and pm1["status"] == "met"
    assert pm1["criteria_strength"] == "PM1_Moderate"
