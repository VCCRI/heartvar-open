"""Unit tests for backend.acmg.hard_coded.carrier_status / gene_inheritance_modes.

The deterministic, INFORMATIONAL carrier-status note resolves mode of inheritance
PER DISEASE (ClinGen Gene-Disease Validity + adequate-tier GenCC) and scopes the
recessive/dominant call to the disease that matches the proband's phenotype. So a
gene's off-target (e.g. non-cardiac) recessive association never bleeds into a
cardiac proband's interpretation, while a genuine recessive form of the matched
cardiac disease is surfaced as dual inheritance. With no phenotype (or no match)
it falls back to the gene-wide picture and discloses which conditions are dominant
vs recessive. It NEVER changes the ACMG variant classification.

The non-negotiable safety property: a single het P/LP whose phenotype matches a
disease with BOTH an adequate recessive AND dominant mode must NEVER read as a
bare "unaffected carrier" — it reads dual-risk.

Phenotypes are passed as free text (matched via the phrase map) so the tests don't
depend on the JAX HPO descendant cache being warm.
"""
from __future__ import annotations

from backend.app import carrier_status, gene_inheritance_modes
from backend.acmg.hard_coded import (
    _norm_gencc_moi, _norm_clingen_gv_moi, _norm_panelapp_moi,
)

LP = "Likely pathogenic"
P = "Pathogenic"

HCM = "hypertrophic cardiomyopathy"
DCM = "dilated cardiomyopathy"
CPVT = "catecholaminergic polymorphic ventricular tachycardia"
LQT = "long QT syndrome"
ARVC = "arrhythmogenic right ventricular cardiomyopathy"


def cc(**kw) -> dict:
    base = dict(
        inheritance_input="", zygosity="", in_trans_pathogenic="",
        denovo_status="", proband_sex="", chromosome="", zygosity_inferred=False,
    )
    base.update(kw)
    return base


def state(gene, classification, clin, hpo, ev=None):
    r = carrier_status(ev or {}, clin, gene, classification, hpo)
    return r["state"] if r else None


def test_norm_gencc_moi():
    assert _norm_gencc_moi("Autosomal recessive") == "AR"
    assert _norm_gencc_moi("Autosomal dominant") == "AD"
    assert _norm_gencc_moi("X-linked recessive") == "XLR"
    assert _norm_gencc_moi("X-linked") == "XL"
    assert _norm_gencc_moi("Semidominant") == "SD"
    assert _norm_gencc_moi("Mitochondrial") == "MT"
    assert _norm_gencc_moi("Unknown") is None
    assert _norm_gencc_moi("") is None
    assert _norm_gencc_moi("some new upstream label") is None


def test_norm_clingen_gv_moi():
    assert _norm_clingen_gv_moi("AR") == "AR"
    assert _norm_clingen_gv_moi("AD") == "AD"
    assert _norm_clingen_gv_moi("XL") == "XL"
    assert _norm_clingen_gv_moi("SD") == "SD"
    assert _norm_clingen_gv_moi("UD") is None
    assert _norm_clingen_gv_moi("???") is None


def test_norm_panelapp_moi():
    assert _norm_panelapp_moi("BIALLELIC, autosomal or pseudoautosomal") == ["AR"]
    assert _norm_panelapp_moi("MONOALLELIC, autosomal or pseudoautosomal") == ["AD"]
    assert set(_norm_panelapp_moi("BOTH monoallelic and biallelic")) == {"AR", "AD"}
    assert _norm_panelapp_moi(
        "X-LINKED: hemizygous in males, biallelic in females") == ["XL"]


_MYH7_EV = {"gencc": {"found": True, "submissions": [
    {"disease": "hypertrophic cardiomyopathy", "classification": "Definitive",
     "moi": "Autosomal dominant", "submitter": "ClinGen"},
    {"disease": "myopathy, myosin storage, autosomal recessive",
     "classification": "Strong", "moi": "Autosomal recessive", "submitter": "Labcorp"},
]}}


def test_cardiac_phenotype_scopes_out_noncardiac_recessive():
    assert carrier_status(_MYH7_EV, cc(zygosity="het"), "MYH7", P, HCM) is None


def test_no_phenotype_discloses_dual_and_names_conditions():
    r = carrier_status(_MYH7_EV, cc(zygosity="het"), "MYH7", P, "")
    assert r is not None and r["state"] == "dual_risk"
    assert r["severity"] == "info"
    assert "myosin storage" in r["detail"].lower()


def test_matched_cardiac_recessive_is_genuine_dual():
    ev = {"gencc": {"found": True, "submissions": [
        {"disease": "arrhythmogenic right ventricular dysplasia 8",
         "classification": "Definitive", "moi": "Autosomal recessive", "submitter": "G2P"},
    ]}}
    r = carrier_status(ev, cc(zygosity="het"), "DSP", P, ARVC)
    assert r is not None and r["state"] == "dual_risk"
    assert r["severity"] == "warning"


def test_tier_gate_blocks_benign_only():
    assert carrier_status({}, cc(zygosity="het"), "TECRL", "Benign", CPVT) is None
    assert carrier_status({}, cc(zygosity="het"), "TECRL", "Likely benign", CPVT) is None
    assert carrier_status({}, cc(zygosity="het"), "TECRL", "VUS", CPVT) is not None


def test_vus_carrier_is_conditional_and_info():
    r = carrier_status({}, cc(zygosity="het"), "TECRL", "VUS", CPVT)
    assert r["state"] == "carrier"
    assert r["severity"] == "info"
    assert "if reclassified pathogenic" in r["detail"].lower()
    assert "uncertain significance" in r["detail"].lower()


def test_vus_dual_gene_is_dual_risk_but_info_severity():
    r = carrier_status({}, cc(zygosity="het"), "CASQ2", "VUS", CPVT)
    assert r["state"] == "dual_risk"
    assert r["severity"] == "info"
    assert "if reclassified pathogenic" in r["detail"].lower()


def test_vus_dominant_only_gene_still_no_band():
    assert carrier_status({}, cc(zygosity="het"), "MYBPC3", "VUS", HCM) is None


def test_vus_xlinked_female_is_info_not_warning():
    r = carrier_status(
        {}, cc(zygosity="het", proband_sex="female", chromosome="X"), "GLA", "VUS",
        "Fabry disease")
    assert r["state"] == "xlinked_female_carrier"
    assert r["severity"] == "info"


def test_recessive_only_het_is_carrier():
    assert state("TECRL", LP, cc(zygosity="het"), CPVT) == "carrier"
    assert state("TRDN", P, cc(zygosity="het"), CPVT) == "carrier"


def test_carrier_band_always_carries_second_allele_caveat():
    r = carrier_status({}, cc(zygosity="het"), "TECRL", LP, CPVT)
    assert r["state"] == "carrier"
    assert "second allele" in r["detail"].lower()
    assert "does not change the acmg" in r["detail"].lower()


def test_dual_genes_are_dual_risk_never_carrier():
    for gene, ph in (("CASQ2", CPVT), ("KCNQ1", LQT)):
        r = carrier_status({}, cc(zygosity="het"), gene, LP, ph)
        assert r is not None and r["state"] == "dual_risk", (gene, r)
        assert r["severity"] == "warning"


def test_dominant_only_gene_gets_no_carrier_band():
    assert carrier_status({}, cc(zygosity="het"), "MYBPC3", P, HCM) is None
    assert carrier_status({}, cc(zygosity="het"), "DES", P, DCM) is None
    assert carrier_status({}, cc(zygosity="het"), "MYH7", P, HCM) is None


def test_homozygous_is_biallelic_affected():
    assert state("TECRL", P, cc(zygosity="hom"), CPVT) == "biallelic_affected"
    assert state("CASQ2", P, cc(zygosity="hom"), CPVT) == "biallelic_affected"


def test_in_trans_is_biallelic_affected():
    assert state("TECRL", P, cc(
        zygosity="het", inheritance_input="AR", in_trans_pathogenic="yes",
    ), CPVT) == "biallelic_affected"


def test_pm3_met_never_co_occurs_with_carrier():
    for clin in (
        cc(zygosity="het", inheritance_input="AR", in_trans_pathogenic="yes"),
        cc(zygosity="het", inheritance_input="AR", denovo_status="inherited_affected"),
    ):
        assert state("TECRL", P, clin, CPVT) == "biallelic_affected"
        assert state("CASQ2", P, clin, CPVT) == "biallelic_affected"


def test_xlinked_male_is_hemizygous_affected():
    assert state("GLA", P, cc(zygosity="hemi", proband_sex="male", chromosome="X"),
                 "Fabry disease") == "xlinked_male"
    assert state("GLA", P, cc(zygosity="het", proband_sex="male", chromosome="X"),
                 "Fabry disease") == "xlinked_male"


def test_xlinked_female_het_is_manifesting_carrier():
    r = carrier_status({}, cc(zygosity="het", proband_sex="female", chromosome="X"),
                       "GLA", P, "Fabry disease")
    assert r["state"] == "xlinked_female_carrier"
    assert r["severity"] == "warning"
    assert "x-inactivation" in r["detail"].lower()


def test_xlinked_unknown_sex_is_uncertain():
    assert state("GLA", P, cc(zygosity="het", chromosome="X"), "Fabry disease") == "uncertain"


def test_no_phenotype_dual_gene_discloses():
    r = carrier_status({}, cc(zygosity="het"), "CASQ2", LP, "")
    assert r["state"] == "dual_risk"
    assert r["severity"] == "info"


def test_no_phenotype_recessive_only_is_uncertain():
    assert state("TECRL", LP, cc(zygosity="het"), "") == "uncertain"


def test_unmatched_phenotype_falls_back_to_gene_wide():
    assert state("TECRL", LP, cc(zygosity="het"), HCM) == "uncertain"


def test_unknown_zygosity_is_uncertain():
    assert state("TECRL", LP, cc(zygosity=""), CPVT) == "uncertain"


def test_unknown_gene_is_uncertain_or_none():
    assert state("ZZZ_NOT_A_GENE", P, cc(zygosity="het"), HCM) == "uncertain"


def test_mito_and_y_are_not_applicable():
    assert carrier_status({}, cc(zygosity="het", inheritance_input="MT"), "MT-TL1", P, CPVT) is None
    assert carrier_status({}, cc(zygosity="het", chromosome="MT"), "SOMEGENE", P, CPVT) is None
    assert carrier_status({}, cc(zygosity="het", chromosome="Y"), "SOMEGENE", P, CPVT) is None


def test_rollup_canonical_genes():
    assert gene_inheritance_modes({}, "CASQ2")["has_dominant_signal"] is True
    assert gene_inheritance_modes({}, "CASQ2")["recessive_established"] is True
    assert gene_inheritance_modes({}, "TECRL")["has_dominant_signal"] is False
    assert gene_inheritance_modes({}, "TECRL")["recessive_established"] is True
    assert gene_inheritance_modes({}, "GLA")["xlinked"] is True


def test_facts_carry_per_disease_category():
    facts = gene_inheritance_modes(_MYH7_EV, "MYH7")["facts"]
    rec = [f for f in facts if f["mode"] == "AR"]
    assert rec and all(f["category"] is None for f in rec)
    assert any(f["category"] == "hcm" for f in facts if f["mode"] == "AD")


def test_gencc_evidence_supplements_rollup():
    ev = {"gencc": {"found": True, "submissions": [
        {"disease": "dilated cardiomyopathy", "classification": "Definitive",
         "moi": "Autosomal recessive", "submitter": "Ambry"},
    ]}}
    m = gene_inheritance_modes(ev, "FAKEREC1")
    assert m["recessive_established"] is True
    assert m["has_dominant_signal"] is False
    assert carrier_status(ev, cc(zygosity="het"), "FAKEREC1", LP, DCM)["state"] == "carrier"


def test_orphanet_submissions_excluded_from_rollup():
    ev = {"gencc": {"found": True, "submissions": [
        {"disease": "dilated cardiomyopathy", "classification": "Definitive",
         "moi": "Autosomal recessive", "submitter": "Orphanet"},
    ]}}
    assert gene_inheritance_modes(ev, "FAKEORP1")["recessive_established"] is False


def test_supportive_tier_excluded_from_rollup():
    ev = {"gencc": {"found": True, "submissions": [
        {"disease": "dilated cardiomyopathy", "classification": "Supportive",
         "moi": "Autosomal recessive", "submitter": "Labcorp"},
    ]}}
    assert gene_inheritance_modes(ev, "FAKESUP1")["recessive_established"] is False
