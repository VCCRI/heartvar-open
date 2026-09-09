"""Unit tests for the gene-specific VCEP frequency thresholds driving PM2/BS1.

Guards the coupled PM2-ceiling / BS1-floor table (backend/data/
vcep_frequency_thresholds.json) wired into compute_hard_coded_criteria. The
generic 1e-4 PM2 / prevalence-scaled BS1 over-fired PM2 (precision 66%, 11 FP
in the 2026-06-04 eRepo re-baseline) and missed BS1 on low-prevalence cardiac/
RASopathy genes. These tests pin the per-VCEP firing behaviour, the RASopathy
absence-only PM2 rule, the LZTR1 AR carve-out, and the generic fallback.

No pytest dependency — runnable with pytest or directly
(``python -m backend.tests.test_vcep_frequency_thresholds``).
"""
from __future__ import annotations

from backend.app import (
    _VCEP_FREQ,
    _pm2_threshold,
    _bs1_threshold,
    compute_hard_coded_criteria,
    _criterion_applicable,
    _gate_criteria_applicability,
    _parse_pm1_hotspots,
)


def _ev(af):
    """gnomAD evidence with FAF95 popmax == af (af=0 -> absent)."""
    return {"gnomad": {"ok": True, "variant": {"exome": {"faf95": {"popmax": af}}, "genome": {}}}}


def _ev_pops(pops, af=None):
    """gnomAD evidence carrying a PER-POPULATION breakdown.

    KCNQ1 (GN112), FBN1 (GN022) and BMPR2 (GN125) are the spec'd genes that
    compare a per-population POINT allele frequency rather than the filtering AF,
    so they need this shape — a FAF95-only fixture has nothing for them to read,
    which is the whole distinction the routing implements. (BMPR2 was added
    2026-09-01; before that this docstring correctly said "the only two".)

    ``pops`` is {population_id: (ac, an)}.
    """
    block = {"populations": [{"id": pid, "ac": ac, "an": an}
                             for pid, (ac, an) in pops.items()]}
    if af is not None:
        block["af"] = af
    return {"gnomad": {"ok": True, "variant": {"exome": block, "genome": None}}}


def _crit_ev(gene, ev, inh="AD"):
    out = compute_hard_coded_criteria(ev, {"inheritance_input": inh}, gene)
    return {e["code"]: e for e in out}


def _crit(gene, af, inh="AD"):
    """Run the hard-coded criteria and return {code: entry}."""
    out = compute_hard_coded_criteria(_ev(af), {"inheritance_input": inh}, gene)
    return {e["code"]: e for e in out}


def _met(c, code):
    return c[code]["status"] == "met"


def test_cardio_mybpc3_pm2_fires_when_rare():
    c = _crit("MYBPC3", 3e-5)
    assert _met(c, "PM2") and c["PM2"]["criteria_strength"] == "PM2_Supporting"
    assert not _met(c, "BS1")


def test_cardio_mybpc3_gap_zone_neither_fires():
    c = _crit("MYBPC3", 5e-5)
    assert not _met(c, "PM2"), "5e-5 is above MYBPC3 PM2 ceiling 4e-5"
    assert not _met(c, "BS1"), "5e-5 is below MYBPC3 BS1 floor 2e-4"


def test_cardio_mybpc3_bs1_fires_and_precludes_pm2():
    c = _crit("MYBPC3", 3e-4)
    assert _met(c, "BS1")
    assert not _met(c, "PM2"), "BS1 must preclude PM2"


def test_cardio_myh7_bs1_floor_lower_than_mybpc3():
    assert _met(_crit("MYH7", 1.5e-4), "BS1")
    assert not _met(_crit("MYBPC3", 1.5e-4), "BS1")


def test_rasopathy_pm2_absence_only_present_rare_does_not_fire():
    c = _crit("HRAS", 1e-5)
    assert not _met(c, "PM2"), "RASopathy PM2 is absence-only; 1e-5 is present"


def test_rasopathy_pm2_fires_when_truly_absent():
    c = _crit("HRAS", 0)
    assert _met(c, "PM2") and c["PM2"]["criteria_strength"] == "PM2_Supporting"


def test_rasopathy_bs1_fires_above_floor():
    assert _met(_crit("MAP2K2", 3e-4), "BS1")
    assert not _met(_crit("MAP2K2", 1e-4), "BS1")


def test_lztr1_ar_carveout():
    assert _met(_crit("LZTR1", 1e-5, inh="AR"), "PM2")
    assert not _met(_crit("LZTR1", 1e-5, inh="AD"), "PM2")


def test_rasopathy_added_genes_absence_only():
    for gene in ("NRAS", "RAF1", "RIT1", "SOS1", "SOS2", "PTPN11", "SHOC2", "RRAS2"):
        assert _met(_crit(gene, 0), "PM2"), f"{gene} PM2 should fire when absent"
        assert not _met(_crit(gene, 1e-5), "PM2"), f"{gene} PM2 is absence-only"
        assert _met(_crit(gene, 3e-4), "BS1"), f"{gene} BS1 fires above 2.5e-4"


def test_fbn1_low_pm2_ceiling_and_bs1_floor():
    assert _met(_crit_ev("FBN1", _ev_pops({"nfe": (1, 250_000)})), "PM2")
    assert not _met(_crit_ev("FBN1", _ev_pops({"nfe": (1, 100_000)})), "PM2")
    assert _met(_crit_ev("FBN1", _ev_pops({"nfe": (1, 10_000)})), "BS1")
    assert not _met(_crit_ev("FBN1", _ev_pops({"nfe": (1, 33_334)})), "BS1")


def test_fbn1_ignores_the_populations_its_spec_excludes():
    """Finnish, Ashkenazi Jewish and "Other" are barred by GN022. gnomAD v4
    renamed oth to remaining, so both spellings have to be handled or the
    exclusion silently stops working on the current dataset."""
    for barred in ("fin", "asj", "oth", "remaining"):
        ev = _ev_pops({barred: (10, 10_000), "nfe": (1, 1_000_000)})
        assert not _met(_crit_ev("FBN1", ev), "BS1"), barred


def test_fbn1_drops_a_population_below_its_allele_minimum():
    """GN022 requires >=2,000 studied alleles in the population used. An AF off a
    handful of alleles is noise, and this rule fires in the BENIGN direction."""
    thin = _ev_pops({"amr": (1, 500)})
    assert not _met(_crit_ev("FBN1", thin), "BS1")
    thick = _ev_pops({"amr": (1, 5_000)})
    assert _met(_crit_ev("FBN1", thick), "BS1")


def test_fbn1_r464c_is_a_real_v2_versus_v4_conflict():
    """The variant behind the decision, with its real gnomAD v4.1 counts. mid
    (Middle Eastern) is not one of the excluded populations and DID NOT EXIST in
    gnomAD v2, when this threshold was calibrated; 5,768 alleles also clears the
    2,000 minimum. So a faithful v4 implementation fires BS1 on a variant the
    expert panel calls Likely Pathogenic. Accepted deliberately — recorded here
    so it cannot become a surprise later."""
    ev = _ev_pops({"mid": (1, 5_768), "nfe": (2, 1_111_882), "fin": (0, 53_416)})
    c = _crit_ev("FBN1", ev)
    assert _met(c, "BS1")
    assert "mid" in c["BS1"]["evidence"], c["BS1"]["evidence"]


def test_kcnq1_pm2_and_bs1():
    assert _met(_crit_ev("KCNQ1", _ev_pops({"nfe": (1, 200_000)})), "PM2")
    assert not _met(_crit_ev("KCNQ1", _ev_pops({"nfe": (1, 20_000)})), "PM2")
    assert _met(_crit_ev("KCNQ1", _ev_pops({"nfe": (1, 2_000)})), "BS1")


def test_kcnq1_counts_only_its_five_continental_populations():
    for outside in ("fin", "mid", "ami", "asj", "remaining"):
        ev = _ev_pops({outside: (10, 10_000), "nfe": (1, 1_000_000)})
        assert not _met(_crit_ev("KCNQ1", ev), "BS1"), outside
    for allowed in ("afr", "eas", "nfe", "amr", "sas"):
        ev = _ev_pops({allowed: (10, 10_000)})
        assert _met(_crit_ev("KCNQ1", ev), "BS1"), allowed


def test_kcnq1_g589d_stays_pathogenic_on_the_specs_own_wording():
    """The Finnish founder variant, with its real gnomAD v4.1 counts. Finnish is
    27/52,964 = 5.10e-4, above the 4e-4 BS1 floor — and Finnish is not one of the
    five populations GN112 names. The highest allowed population is nfe at
    8.99e-7, 445x below the floor. No founder-exclusion special case needed: the
    spec's own wording does it."""
    ev = _ev_pops({"fin": (27, 52_964), "nfe": (1, 1_111_994)})
    c = _crit_ev("KCNQ1", ev)
    assert not _met(c, "BS1"), c["BS1"]["evidence"]
    assert not _met(c, "BA1")


def test_bmpr2_is_routed_to_its_own_thresholds_not_the_generic_default():
    """BMPR2 is 82 of the 852 cardiac-VCEP records and had NO entry here, so it
    took the generic path. Only BA1 was actually wrong, and badly: the generic
    ACMG 0.05 against GN125's 'above 1% in gnomAD', 5x too loose. BS1 and PM2
    coincided with the generic AD defaults."""
    assert _VCEP_FREQ["BMPR2"]["ba1"] == 0.01
    assert _pm2_threshold("BMPR2", "AD") == 0.0001
    assert _bs1_threshold({}, "BMPR2", "AD") == 0.001
    ev = _ev_pops({"nfe": (2_000, 100_000)})
    assert _met(_crit_ev("BMPR2", ev), "BA1"), "2% must reach BA1 under GN125"


def test_bmpr2_uses_a_point_af_with_a_1000_allele_minimum():
    """GN125 compares 'the subpopulation with the highest frequency and at least
    1,000 allele counts' — a point AF, not FAF95, making BMPR2 the third
    point_af gene after FBN1 and KCNQ1. Unlike those two it has no population
    allow-list and no exclusions, so the ONLY filter is the allele minimum."""
    spec = _VCEP_FREQ["BMPR2"]
    assert spec["metric"] == "point_af"
    assert spec["min_alleles"] == 1000
    assert "populations" not in spec and "exclude_populations" not in spec

    ev = _ev_pops({"amr": (1, 500), "nfe": (1, 800_000)})
    assert not _met(_crit_ev("BMPR2", ev), "BS1")
    ev = _ev_pops({"amr": (5, 2_000), "nfe": (1, 800_000)})
    assert _met(_crit_ev("BMPR2", ev), "BS1")


def test_ba1_can_fire_on_a_point_af_gene():
    """cc29ab6 (2026-08-26) added point_af routing but left BA1's allele-number
    gate reading the BLOCK-level AN via _gnomad_popmax_an. A point_af fixture
    carries AN per POPULATION and no block AN, so the gate saw None, failed
    closed, and BA1 could never fire on FBN1 or KCNQ1 — live for six days. BS1
    defers when BA1 is withheld, so such variants got NO benign criterion at
    all, a false negative in the toward-pathogenic direction.

    The AN now comes back from _vcep_freq_af alongside the AF, so it always
    describes the dataset that supplied the frequency."""
    ev = _ev_pops({"nfe": (5_000, 100_000)})
    for gene in ("BMPR2", "FBN1", "KCNQ1"):
        assert _met(_crit_ev(gene, ev), "BA1"), gene


def test_the_ba1_allele_number_gate_still_fails_closed_on_a_thin_population():
    """The other half of the fix. Restoring BA1 must not weaken the SVI's
    ">= 2,000 observed alleles" requirement — an AF off a thin cohort is noise,
    and BA1 is stand-alone benign, so a false positive here is a benign call on
    one variant in 1,500 alleles."""
    thin = _ev_pops({"nfe": (75, 1_500)})
    c = _crit_ev("BMPR2", thin)
    assert not _met(c, "BA1"), c["BA1"]["evidence"]
    assert "2,000" in c["BA1"]["evidence"], c["BA1"]["evidence"]


def test_faf_genes_are_untouched_by_the_routing():
    """The Cardiomyopathy and RASopathy VCEPs specify the filtering allele
    frequency, so FAF95 popmax is already conformant for 23 of the 25 spec'd
    genes. A blanket metric switch would break the 619 of 852 records it is right
    for in order to fix 118 — which is why the routing is per-VCEP."""
    for gene in ("MYH7", "MYBPC3", "PTPN11", "RAF1", "TNNT2"):
        assert _met(_crit(gene, 3e-4), "BS1"), gene
        assert _met(_crit(gene, 0), "PM2"), gene


def test_uncovered_gene_uses_generic_threshold():
    assert _pm2_threshold("FOO", "AD") == 0.0001
    assert _met(_crit("FOO", 5e-5), "PM2"), "5e-5 < generic 1e-4"
    assert not _met(_crit("FOO", 2e-4), "PM2"), "2e-4 > generic 1e-4"


def _ev_constraint(pli=None, loeuf=None):
    """gnomAD evidence carrying gene constraint (pLI / LOEUF=oe_lof_upper)."""
    return {"gnomad": {"ok": True, "gene": {"gnomad_constraint": {"pLI": pli, "oe_lof_upper": loeuf}}}}


def test_constraint_aware_tightens_uncovered_lof_intolerant_gene():
    ev = _ev_constraint(pli=0.99, loeuf=0.2)
    assert _pm2_threshold("NEWGENE", "AD", ev) == 0.00004
    assert _bs1_threshold(ev, "NEWGENE", "AD") == 0.0001


def test_constraint_aware_loeuf_only_triggers():
    assert _pm2_threshold("NEWGENE", "AD", _ev_constraint(loeuf=0.30)) == 0.00004
    assert _pm2_threshold("NEWGENE", "AD", _ev_constraint(loeuf=0.50)) == 0.0001


def test_unconstrained_uncovered_gene_keeps_generic():
    ev = _ev_constraint(pli=0.10, loeuf=1.20)
    assert _pm2_threshold("NEWGENE", "AD", ev) == 0.0001
    assert _bs1_threshold(ev, "NEWGENE", "AD") == 0.001


def test_constraint_not_applied_to_ar_or_vcep_genes():
    ev = _ev_constraint(pli=0.99, loeuf=0.1)
    assert _pm2_threshold("NEWGENE", "AR", ev) == 0.01
    assert _pm2_threshold("KRAS", "AD", ev) == 0.0


def _ev_noncoding(consequence, spliceai_ok=True, max_delta=0.0):
    """Evidence for a non-coding variant: gnomAD absent, given VEP consequence,
    SpliceAI present/absent with the given max delta, no missense scores."""
    return {
        "gnomad": {"ok": True, "variant": {"exome": {"faf95": {"popmax": 0.0}}, "genome": {}}},
        "vep": {"most_severe_consequence": consequence},
        "spliceai": {"ok": spliceai_ok, "max_delta": max_delta},
    }


def _hc(ev, gene="SOMEGENE", inh="AD"):
    return {e["code"]: e for e in compute_hard_coded_criteria(ev, {"inheritance_input": inh}, gene)}


def test_bp4_noncoding_low_spliceai_fires_supporting():
    c = _hc(_ev_noncoding("intron_variant", spliceai_ok=True, max_delta=0.02))
    assert c["BP4"]["status"] == "met"
    assert c["BP4"]["criteria_strength"] == "BP4_Supporting"


def test_bp4_noncoding_requires_spliceai_present():
    c = _hc(_ev_noncoding("intron_variant", spliceai_ok=False, max_delta=0.0))
    assert c["BP4"]["status"] == "not_met"


def test_bp4_noncoding_high_spliceai_does_not_fire():
    c = _hc(_ev_noncoding("intron_variant", spliceai_ok=True, max_delta=0.8))
    assert c["BP4"]["status"] == "not_met"


def test_bp4_noncoding_excludes_synonymous_and_canonical_splice():
    assert _hc(_ev_noncoding("synonymous_variant", max_delta=0.01))["BP4"]["status"] == "not_met"
    assert _hc(_ev_noncoding("splice_donor_variant", max_delta=0.01))["BP4"]["status"] == "not_met"


def _ev_hom(hom):
    """gnomAD evidence with `hom` homozygotes (for the generic BS2 hom rule)."""
    return {"gnomad": {"ok": True, "variant": {"exome": {"faf95": {"popmax": 0.0}, "ac_hom": hom}, "genome": {}}}}


def _bs2(ev, gene, zygosity="hom", inh="AD"):
    out = compute_hard_coded_criteria(ev, {"zygosity": zygosity, "inheritance_input": inh}, gene)
    return {e["code"]: e for e in out}["BS2"]


def test_bs2_not_applicable_for_cardiomyopathy_vcep():
    e = _bs2(_ev_hom(20), "MYH7")
    assert e["status"] == "not_met" and "Not Applicable" in e["evidence"]


def test_bs2_rasopathy_needs_family_phenotyping_not_gnomad():
    e = _bs2(_ev_hom(20), "KRAS")
    assert e["status"] == "not_met" and "family" in e["evidence"].lower()


def test_bs2_generic_hom_rule_preserved_for_uncovered_gene():
    assert _bs2(_ev_hom(20), "SOMEGENE")["criteria_strength"] == "BS2_Strong"
    assert _bs2(_ev_hom(3), "SOMEGENE")["criteria_strength"] == "BS2_Strong"
    assert _bs2(_ev_hom(2), "SOMEGENE")["status"] == "not_met"
    assert _bs2(_ev_hom(1), "SOMEGENE")["status"] == "not_met"


def test_criterion_applicable_flags():
    assert _criterion_applicable("ACTC1", "PP2") is False
    assert _criterion_applicable("TPM1", "PP2") is True
    assert _criterion_applicable("PTPN11", "BP1") is True
    assert _criterion_applicable("MYH7", "BP1") is False
    assert _criterion_applicable("SOMEGENE", "PP2") is True


def test_gate_demotes_not_applicable_pp2():
    crit = [{"code": "PP2", "status": "met", "criteria_strength": "PP2",
             "evidence": "missense in constrained gene"}]
    out = {c["code"]: c for c in _gate_criteria_applicability(crit, "ACTC1")}
    assert out["PP2"]["status"] == "not_met" and out["PP2"]["criteria_strength"] is None
    assert "Not Applicable" in out["PP2"]["evidence"]


def test_pvs1_cspec_override_excludes_na_gene_even_for_truncating():
    from backend.app import _gene_lof_mechanism
    chd = {"chdgene": {"listed": True, "inheritance": ["AD"]}}
    assert _gene_lof_mechanism(chd, "HRAS", "stop_gained") is False
    assert _gene_lof_mechanism(chd, "MYBPC3", "splice_acceptor_variant") is True
    assert _gene_lof_mechanism(chd, "SCN1A", "stop_gained") is True


def test_gate_demotes_bp2_on_kcnq1():
    """BP2 is DETERMINISTIC (_eval_bp2), which is why the data flag alone was
    not enough: _eval_bp2 takes no gene argument and never consults
    applicability, so KCNQ1's withdrawal ("Not applicable to KCNQ1 due to
    biallelic cases (Jervell and Lange-Nielsen syndrome)") only bites once BP2
    is in _APPLICABILITY_GATED_CODES. A KCNQ1 variant in trans with a
    pathogenic KCNQ1 allele is Jervell and Lange-Nielsen syndrome, not benign
    evidence, and BP2 is the only benign criterion that fires on that
    genotype — a lone -1 is Likely benign under the point system."""
    from backend.acmg.hard_coded import _APPLICABILITY_GATED_CODES
    assert "BP2" in _APPLICABILITY_GATED_CODES
    assert _criterion_applicable("KCNQ1", "BP2") is False
    crit = [{"code": "BP2", "status": "met", "criteria_strength": "BP2_Supporting",
             "evidence": "in trans with a pathogenic allele"}]
    out = {c["code"]: c for c in _gate_criteria_applicability(crit, "KCNQ1")}
    assert out["BP2"]["status"] == "not_met"
    assert out["BP2"]["criteria_strength"] is None
    assert "Not Applicable" in out["BP2"]["evidence"]
    from backend.acmg.hard_coded import _VCEP_CRIT
    na = sorted(g for g in _VCEP_CRIT
                if (_VCEP_CRIT[g].get("BP2") or {}).get("applicability")
                == "not_applicable")
    assert na == ["KCNQ1"], na
    keep = [{"code": "BP2", "status": "met", "criteria_strength": "BP2_Supporting",
             "evidence": "x"}]
    assert _gate_criteria_applicability(keep, "MYH7")[0]["status"] == "met"


def test_bs2_suppression_on_kcnq1_cites_the_criteria_spec_not_the_freq_table():
    """_eval_bs2 consulted vcep_frequency_thresholds.json BEFORE the criteria
    spec, so KCNQ1's BS2 was closed by a branch holding no BS2 evidence — the
    frequency table records PM2/BS1/BA1 cutoffs only. The answer was right by
    luck. Now the spec row decides, and the evidence string names the gene
    whose spec withdrew it."""
    e = _bs2(_ev_hom(20), "KCNQ1")
    assert e["status"] == "not_met"
    assert e["evidence"].startswith("KCNQ1 VCEP marks BS2 Not Applicable")
    assert _bs2(_ev_hom(20), "KRAS")["status"] == "not_met"
    assert "family" in _bs2(_ev_hom(20), "KRAS")["evidence"].lower()
    for gene in ("MYH7", "MYBPC3", "TNNI3", "TNNT2", "TPM1", "ACTC1", "MYL2",
                 "MYL3", "FBN1"):
        assert _bs2(_ev_hom(20), gene)["status"] == "not_met", gene
    assert _criterion_applicable("BMPR2", "BS2") is True
    assert _bs2(_ev_hom(20), "BMPR2")["status"] == "not_met"


def test_gate_keeps_applicable_and_unknown():
    crit = [{"code": "PP2", "status": "met", "criteria_strength": "PP2", "evidence": "x"}]
    assert _gate_criteria_applicability(crit, "TPM1")[0]["status"] == "met"
    crit2 = [{"code": "PP2", "status": "met", "criteria_strength": "PP2", "evidence": "x"}]
    assert _gate_criteria_applicability(crit2, "SOMEGENE")[0]["status"] == "met"


def _pm1_met():
    return [{"code": "PM1", "status": "met", "criteria_strength": "PM1",
             "evidence": "residue in EGF-like domain; P/LP density"}]


def test_pm1_demoted_on_lof_consequence():
    for cq in ("stop_gained", "frameshift_variant", "splice_donor_variant",
               "splice_acceptor_variant", "start_lost"):
        ev = {"vep": {"most_severe_consequence": cq}}
        out = _gate_criteria_applicability(_pm1_met(), "JAG1", ev)[0]
        assert out["status"] == "not_met", cq
        assert out["criteria_strength"] is None
        assert "loss-of-function" in out["evidence"]


def test_pm1_kept_on_missense_consequence():
    ev = {"vep": {"most_severe_consequence": "missense_variant"}}
    assert _gate_criteria_applicability(_pm1_met(), "JAG1", ev)[0]["status"] == "met"
    ev2 = {"vep": {"most_severe_consequence": "inframe_deletion"}}
    assert _gate_criteria_applicability(_pm1_met(), "JAG1", ev2)[0]["status"] == "met"


def test_pm1_lof_demotion_fails_open_without_consequence():
    assert _gate_criteria_applicability(_pm1_met(), "JAG1", {})[0]["status"] == "met"
    assert _gate_criteria_applicability(_pm1_met(), "JAG1", None)[0]["status"] == "met"


def test_pm1_hotspot_parser_captures_all_ranges_and_exons():
    rngs, exons = _parse_pm1_hotspots("exon 6, exon 11, P-loop [AA 459-474], CR3 [AA 594-627]")
    assert (459, 474) in rngs and (594, 627) in rngs
    assert exons == {6, 11}
    rngs2, exons2 = _parse_pm1_hotspots("Codons 485-502 and 1248-1266")
    assert rngs2 == [(485, 502), (1248, 1266)] and exons2 == set()


def test_pm1_asserted_for_exon_defined_hotspot():
    crit = [{"code": "PM2", "status": "met", "criteria_strength": "PM2_Supporting",
             "evidence": "rare"}]
    ev = {"vep": {"most_severe_consequence": "missense_variant", "exon": "6/18"}}
    out = {c["code"]: c for c in _gate_criteria_applicability(crit, "BRAF", ev)}
    assert out.get("PM1", {}).get("status") == "met", "PM1 should assert on exon-6 hotspot"
    assert out["PM1"]["criteria_strength"] == "PM1_Moderate"
    assert "exon 6" in out["PM1"]["evidence"]


def test_pm1_not_asserted_outside_hotspot_exon():
    crit = [{"code": "PM2", "status": "met", "criteria_strength": "PM2_Supporting",
             "evidence": "rare"}]
    ev = {"vep": {"most_severe_consequence": "missense_variant", "exon": "3/18"}}
    out = {c["code"]: c for c in _gate_criteria_applicability(crit, "BRAF", ev)}
    assert "PM1" not in out or out["PM1"]["status"] != "met", "exon 3 is not a BRAF hotspot"


def _pp4_met():
    return [{"code": "PP4", "status": "met", "criteria_strength": "PP4",
             "evidence": "green CHD panel + hpo_match; family has Alagille"}]


def test_pp4_suppressed_where_the_vcep_marks_it_not_applicable():
    """THE test that would have caught the bug, in the direction that matters.
    PP4 used to be demoted unconditionally, in a branch sitting AFTER this
    per-gene applicability check — so the flag below was dead code and the
    KCNQ1 GN112 PP4 rule could never fire. These four genes are the ones whose
    CSpec really does mark PP4 Not Applicable, and they must still suppress."""
    for gene in ("MYH7", "MYBPC3", "PTPN11", "TNNT2"):
        out = _gate_criteria_applicability(_pp4_met(), gene, {})[0]
        assert out["status"] == "not_met", f"{gene} must suppress PP4"
        assert out["criteria_strength"] is None
        assert "Not Applicable" in out["evidence"]


def test_pp4_survives_where_the_vcep_marks_it_applicable():
    """KCNQ1 GN112 spells out an affirmative PP4 rule (QT >480 ms AND a
    swimming-associated event / treadmill result / LQT1 T-wave morphology) and
    experts applied PP4 on 5 of 37 KCNQ1 eRepo records. It must be reachable."""
    for ev in (
        {"vep": {"most_severe_consequence": "stop_gained"}},
        {"vep": {"most_severe_consequence": "missense_variant"}},
        {},
        None,
    ):
        out = _gate_criteria_applicability(_pp4_met(), "KCNQ1", ev)[0]
        assert out["status"] == "met", (
            "KCNQ1 CSpec marks PP4 Applicable — it must not be demoted"
        )


def test_pp4_is_closed_on_a_gene_with_no_vcep_entry():
    """The JAG1 hole is CLOSED by the affirmative allow-list.

    This test previously asserted the opposite and said so: it pinned the COST of
    fix A3 removing the blanket PP4 demotion, and flagged that "a per-gene
    affirmative allow-list would close it". That allow-list now exists, so the
    assertion flips. JAG1 has no entry in vcep_criteria_applicability.json, so
    under default-open semantics the original JAG1 c.703C>T over-call (an
    isolated cardiac phenotype treated as specific for Alagille) was reachable
    with only the prompt's phenotype matrix in the way.
    """
    out = _gate_criteria_applicability(_pp4_met(), "JAG1", {})[0]
    assert out["status"] == "not_met"
    assert out["criteria_strength"] is None
    assert "affirmatively" in out["evidence"]


def test_pp4_still_fires_where_a_vcep_affirmatively_allows_it():
    """KCNQ1 GN112 is the only affirmative PP4 entry in the table. The allow-list
    must not re-create the blanket demotion fix A3 removed."""
    out = _gate_criteria_applicability(_pp4_met(), "KCNQ1", {})[0]
    assert out["status"] == "met"


def test_pp4_still_fires_on_fbn1_via_the_evidence_allow_list():
    """FBN1 GN022 has no PP4 key at all, so it is not affirmatively applicable —
    but experts applied PP4 on 63 of its 114 eRepo records. It is allow-listed on
    that evidence, which together with KCNQ1 reproduces cardiac expert practice
    exactly (68 of 852: FBN1 63 + KCNQ1 5, zero elsewhere)."""
    out = _gate_criteria_applicability(_pp4_met(), "FBN1", {})[0]
    assert out["status"] == "met"


def test_pp4_is_suppressed_on_a_gene_the_cspec_marks_not_applicable():
    for gene in ("MYH7", "MYBPC3", "PTPN11"):
        out = _gate_criteria_applicability(_pp4_met(), gene, {})[0]
        assert out["status"] == "not_met", gene


def _pvs1_met():
    return [{"code": "PVS1", "status": "met", "criteria_strength": "PVS1",
             "evidence": "null variant in established LoF gene"}]


def test_pvs1_suppressed_on_the_four_genes_the_mechanism_list_missed():
    """Of the 21 table genes marked PVS1 Not Applicable, _DOMINANT_NEGATIVE_GENES
    already blocks 17. These four leaked, and PVS1 is +8 — the largest single
    over-call the engine can make."""
    for gene in ("RRAS2", "SHOC2", "SOS1", "SOS2"):
        out = _gate_criteria_applicability(_pvs1_met(), gene, {})[0]
        assert out["status"] == "not_met", gene
        assert out["criteria_strength"] is None, gene


def test_pvs1_gating_does_not_touch_genes_outside_the_table():
    """_criterion_applicable defaults to applicable, so the 653 non-VCEP genes
    are unaffected. Suppressing there would gut PVS1 across most of the panel."""
    for gene in ("JAG1", "SOMEGENE", None):
        out = _gate_criteria_applicability(_pvs1_met(), gene, {})[0]
        assert out["status"] == "met", gene


def test_pvs1_still_fires_where_the_cspec_allows_it():
    for gene in ("MYBPC3", "KCNQ1"):
        out = _gate_criteria_applicability(_pvs1_met(), gene, {})[0]
        assert out["status"] == "met", gene


def test_pp2_demoted_on_non_missense():
    crit = [{"code": "PP2", "status": "met", "criteria_strength": "PP2",
             "evidence": "constrained gene"}]
    ev = {"vep": {"most_severe_consequence": "stop_gained"}}
    out = _gate_criteria_applicability(crit, "JAG1", ev)[0]
    assert out["status"] == "not_met" and "missense-only" in out["evidence"]
    crit2 = [{"code": "PP2", "status": "met", "criteria_strength": "PP2", "evidence": "x"}]
    ev2 = {"vep": {"most_severe_consequence": "missense_variant"}}
    assert _gate_criteria_applicability(crit2, "JAG1", ev2)[0]["status"] == "met"
    crit3 = [{"code": "PP2", "status": "met", "criteria_strength": "PP2", "evidence": "x"}]
    assert _gate_criteria_applicability(crit3, "JAG1", {})[0]["status"] == "met"


def _ev_lookup(ok=True, variant_found=False, indel_unresolved=False):
    gn = {"ok": ok, "variant_found": variant_found, "variant": None}
    if indel_unresolved:
        gn["indel_unresolved"] = True
    return {"gnomad": gn}


def _hc_gene(ev, gene="MYBPC3", inh="AD"):
    return {e["code"]: e for e in compute_hard_coded_criteria(ev, {"inheritance_input": inh}, gene)}


def test_pm2_fires_on_confirmed_absence():
    c = _hc_gene(_ev_lookup())
    assert _met(c, "PM2") and c["PM2"]["criteria_strength"] == "PM2_Supporting"


def test_pm2_not_fired_when_lookup_failed():
    c = _hc_gene(_ev_lookup(ok=False))
    assert not _met(c, "PM2"), "failed lookup must not be treated as absent"
    assert "unavailable" in c["PM2"]["evidence"].lower()


def test_pm2_not_fired_when_indel_unresolved():
    c = _hc_gene(_ev_lookup(indel_unresolved=True))
    assert not _met(c, "PM2"), "unresolved indel frequency must not over-call PM2"
    assert "indel" in c["PM2"]["evidence"].lower()


def test_threshold_ordering_invariant():
    for gene, spec in _VCEP_FREQ.items():
        pm2, bs1, ba1 = spec.get("pm2_max"), spec.get("bs1"), spec.get("ba1")
        assert pm2 is not None and bs1 is not None and ba1 is not None, gene
        assert pm2 <= bs1 < ba1, f"{gene}: expected pm2({pm2}) <= bs1({bs1}) < ba1({ba1})"


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
