"""Characterization / regression test for the deterministic ACMG scorer.

This test LOCKS IN the current, benchmark-validated behaviour of the
hard-coded ACMG scoring engine so a future refactor (or an unintended
threshold/gene-list edit) cannot silently change a score, a tier, or which
criteria fire. It is intentionally a *golden* test: the EXPECTED values below
were captured from the live scorer and written out explicitly. It passes now
and FAILS LOUDLY the moment deterministic scoring drifts.

It exercises ``compute_hard_coded_criteria`` (the 14 structurally-evaluated
ACMG codes) followed by ``apply_cross_criterion_exclusions`` (BA1 forced-Benign
demotion + PP5/BP6 guards), then ``compute_points_total`` + ``classification_for``
from the canonical ``backend.acmg.tiers`` module. These are the production code
paths, imported from the canonical modules (NOT the ``backend.app`` re-export
surface), so the test follows the symbols to their real home.

The cases span, at minimum:
  - the met path of every structurally-firable criterion:
    PVS1, PM2, BA1, BS1, BS2, PM3, PM4, BP2, BP3, BP7, PS2, PM6, PP5, BP6,
    PP3, BP4;
  - multiple gene mechanisms: AD / AR / XL, dominant-negative, RASopathy
    (gain-of-function), LoF-intolerant, and off-panel (no VCEP);
  - all five ACMG tiers: Benign, Likely benign, VUS, Likely pathogenic,
    Pathogenic.

Deterministic + fully offline: every input is an in-memory evidence dict; no
network, filesystem, or clock dependency. No external scratch files.

HARD RULE: this test must never be "fixed" by editing the EXPECTED values to
match new output. If it fails, that means scoring changed — investigate WHY,
and only update the goldens deliberately (with a benchmark re-validation).
"""
from __future__ import annotations

import copy

import pytest

from backend.acmg.hard_coded import (
    apply_cross_criterion_exclusions,
    compute_hard_coded_criteria,
)
from backend.acmg.tiers import classification_for, compute_points_total


def ev(
    *,
    gnomad_ok: bool = True,
    faf95: float | None = None,
    an: int | None = None,
    ac_hom: int = 0,
    ac_hemi: int | None = None,
    pli: float | None = None,
    loeuf: float | None = None,
    spliceai: float | None = None,
    consequence: str | None = None,
    exon: str | None = None,
    functionless_repeat: bool = False,
    chdgene_listed: bool = False,
    chdgene_inh: list[str] | None = None,
    clinvar_records: list[dict] | None = None,
    am_score: float | None = None,
    revel: float | None = None,
    cadd: float | None = None,
    phylop: float | None = None,
    confirmed_absent: bool = False,
) -> dict:
    """Build a compute_hard_coded_criteria-shaped evidence dict.

    ``confirmed_absent`` sets the gnomAD ``variant_found=False`` signal PM2
    needs to distinguish a true gnomAD absence from a failed lookup. In-silico
    scores live where the scorer reads them: REVEL/CADD under ``vep`` and
    AlphaMissense under ``alphamissense.am_pathogenicity``.
    """
    e: dict = {}
    variant: dict = {}
    exome: dict = {}
    if faf95 is not None:
        exome["faf95"] = {"popmax": faf95}
        exome["an"] = 1_461_852 if an is None else an
    if ac_hom:
        exome["ac_hom"] = ac_hom
    if ac_hemi is not None:
        exome["ac_hemi"] = ac_hemi
    if exome:
        variant["exome"] = exome
    gene_block: dict = {}
    if pli is not None or loeuf is not None:
        gene_block["gnomad_constraint"] = {"pLI": pli, "oe_lof_upper": loeuf}
    gnomad_block: dict = {
        "ok": gnomad_ok,
        "variant": variant or None,
        "gene": gene_block or None,
    }
    if confirmed_absent:
        gnomad_block["variant_found"] = False
        gnomad_block["indel_unresolved"] = False
    e["gnomad"] = gnomad_block
    if spliceai is not None:
        e["spliceai"] = {
            "ok": True,
            "max_delta": spliceai,
            "scores_per_transcript": [{"DS_DG": spliceai}],
        }
    vep: dict = {"ok": True}
    if consequence is not None:
        vep["most_severe_consequence"] = consequence
    if exon is not None:
        vep["exon"] = exon
    if functionless_repeat:
        vep["protein_start"] = 100
        e["uniprot"] = {
            "ok": True, "found": True, "domains": [],
            "features": [
                {"type": "Repeat", "start": 90, "end": 120,
                 "description": "Gln-rich"},
            ],
        }
    if revel is not None:
        vep["revel_score"] = revel
    if cadd is not None:
        vep["cadd_phred"] = cadd
    if phylop is not None:
        vep["phylop100way"] = phylop
    e["vep"] = vep
    if chdgene_listed:
        e["chdgene"] = {"listed": True, "inheritance": chdgene_inh or []}
    if clinvar_records is not None:
        e["clinvar"] = {
            "ok": True,
            "found": bool(clinvar_records),
            "records": clinvar_records,
        }
    if am_score is not None:
        e["alphamissense"] = {"ok": True, "am_pathogenicity": am_score}
    return e


def cv_rec(sig: str, stars: int, n_sub: int = 1,
           acc: str = "VCV000001", rs: str = "reviewed") -> dict:
    """A single ClinVar record in the shape the PP5/BP6 cluster reads."""
    return {
        "clinical_significance": sig,
        "stars": stars,
        "number_submitters": n_sub,
        "accession": acc,
        "review_status": rs,
    }


GOLDEN_CASES = [
    dict(id="ad_missense_constrained_lookup_absent",
         ev=ev(pli=0.99, consequence="missense_variant", exon="6/18"),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=0, tier="VUS", met=[]),

    dict(id="ad_missense_confirmed_absent_PM2",
         ev=ev(consequence="missense_variant", confirmed_absent=True),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=1, tier="VUS", met=[("PM2", "PM2_Supporting")]),

    dict(id="ar_offpanel_confirmed_absent_PM2",
         ev=ev(consequence="missense_variant", confirmed_absent=True),
         cc={"inheritance_input": "AR"}, gene="ZZZZ9",
         points=1, tier="VUS", met=[("PM2", "PM2_Supporting")]),

    dict(id="ad_missense_BA1_common",
         ev=ev(faf95=0.08, consequence="missense_variant"),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=-8, tier="Benign", met=[("BA1", "BA1")]),

    dict(id="ad_missense_BA1_vcep_default",
         ev=ev(faf95=0.002, consequence="missense_variant"),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=-8, tier="Benign", met=[("BA1", "BA1")]),

    dict(id="offpanel_BS1_only",
         ev=ev(faf95=0.002, consequence="missense_variant"),
         cc={"inheritance_input": "AD"}, gene="ZZZZ9",
         points=-4, tier="Likely benign", met=[("BS1", "BS1_Strong")]),

    dict(id="offpanel_BS2_gnomad_homozygotes",
         ev=ev(consequence="missense_variant", ac_hom=4),
         cc={"inheritance_input": "AD", "zygosity": "het"},
         gene="ZZZZ9",
         points=-4, tier="Likely benign", met=[("BS2", "BS2_Strong")]),

    dict(id="offpanel_BS2_two_homozygotes_withheld",
         ev=ev(consequence="missense_variant", ac_hom=2),
         cc={"inheritance_input": "AD", "zygosity": "het"},
         gene="ZZZZ9",
         points=0, tier="VUS", met=[]),

    dict(id="ad_frameshift_PVS1_lof_gene",
         ev=ev(pli=0.99, consequence="frameshift_variant", exon="3/18",
               chdgene_listed=True, chdgene_inh=["AD"]),
         cc={"inheritance_input": "AD"}, gene="MYBPC3",
         points=8, tier="Likely pathogenic", met=[("PVS1", "PVS1")]),

    dict(id="xl_frameshift_PVS1_lof_gene",
         ev=ev(pli=0.95, consequence="frameshift_variant", exon="2/10",
               chdgene_listed=True, chdgene_inh=["XL"]),
         cc={"inheritance_input": "XL"}, gene="FLNA",
         points=8, tier="Likely pathogenic", met=[("PVS1", "PVS1")]),

    dict(id="ad_splice_donor_PVS1_high_spliceai",
         ev=ev(spliceai=0.95, consequence="splice_donor_variant", pli=0.99,
               chdgene_listed=True, chdgene_inh=["AD"]),
         cc={"inheritance_input": "AD"}, gene="MYBPC3",
         points=8, tier="Likely pathogenic", met=[("PVS1", "PVS1")]),

    dict(id="domneg_missense_no_PVS1",
         ev=ev(consequence="missense_variant", exon="13/40"),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=0, tier="VUS", met=[]),

    dict(id="rasopathy_nonsense_PVS1",
         ev=ev(pli=0.99, consequence="stop_gained", exon="5/15"),
         cc={"inheritance_input": "AD"}, gene="HRAS",
         points=8, tier="Likely pathogenic", met=[("PVS1", "PVS1")]),

    dict(id="rasopathy_missense_no_fire",
         ev=ev(consequence="missense_variant"),
         cc={"inheritance_input": "AD"}, gene="PTPN11",
         points=0, tier="VUS", met=[]),

    dict(id="inframe_del_PM4",
         ev=ev(consequence="inframe_deletion"),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=2, tier="VUS", met=[("PM4", "PM4_Moderate")]),

    dict(id="stop_loss_PM4",
         ev=ev(consequence="stop_lost"),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=2, tier="VUS", met=[("PM4", "PM4_Moderate")]),

    dict(id="inframe_in_repeat_BP3",
         ev=ev(consequence="inframe_deletion", functionless_repeat=True),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=-1, tier="Likely benign", met=[("BP3", "BP3_Supporting")]),

    dict(id="synonymous_low_spliceai_BP7",
         ev=ev(spliceai=0.02, consequence="synonymous_variant", phylop=1.0),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=-1, tier="Likely benign", met=[("BP7", "BP7_Supporting")]),

    dict(id="insilico_pathogenic_PP3",
         ev=ev(consequence="missense_variant", am_score=0.99, revel=0.95,
               cadd=30),
         cc={"inheritance_input": "AD"}, gene="ZZZZ9",
         points=1, tier="VUS", met=[("PP3", "PP3_Supporting")]),

    dict(id="insilico_benign_BP4",
         ev=ev(consequence="missense_variant", am_score=0.05, revel=0.05,
               cadd=2),
         cc={"inheritance_input": "AD"}, gene="ZZZZ9",
         points=-1, tier="Likely benign", met=[("BP4", "BP4_Supporting")]),

    dict(id="denovo_confirmed_PS2",
         ev=ev(consequence="missense_variant"),
         cc={"inheritance_input": "AD", "denovo_status": "confirmed",
             "denovo_confirmed": True,
             "trio_status": "trio"},
         gene="MYH7",
         points=4, tier="VUS", met=[("PS2", "PS2_Strong")]),

    dict(id="denovo_assumed_PM6",
         ev=ev(consequence="missense_variant"),
         cc={"inheritance_input": "AD", "trio_status": "trio",
             "denovo_status": "unconfirmed"},
         gene="MYH7",
         points=2, tier="VUS", met=[("PM6", "PM6_Moderate")]),

    dict(id="ar_het_in_trans_PM3",
         ev=ev(faf95=0.0001, consequence="missense_variant"),
         cc={"inheritance_input": "AR", "zygosity": "het",
             "in_trans_pathogenic": "yes"},
         gene="ZZZZ9",
         points=3, tier="VUS",
         met=[("PM2", "PM2_Supporting"), ("PM3", "PM3_Moderate")]),

    dict(id="ad_het_in_trans_pathogenic_BP2",
         ev=ev(faf95=0.00005, consequence="missense_variant"),
         cc={"inheritance_input": "AD", "zygosity": "het",
             "in_trans_pathogenic": "yes"},
         gene="ZZZZ9",
         points=0, tier="VUS",
         met=[("BP2", "BP2_Supporting"), ("PM2", "PM2_Supporting")]),

    dict(id="clinvar_3star_PP5_strong",
         ev=ev(consequence="missense_variant",
               clinvar_records=[cv_rec("Pathogenic", 3, 4)]),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=0, tier="VUS", met=[("PP5", "PP5_Strong")]),

    dict(id="clinvar_2star_PP5_moderate",
         ev=ev(consequence="missense_variant",
               clinvar_records=[cv_rec("Likely pathogenic", 2, 3)]),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=0, tier="VUS", met=[("PP5", "PP5_Moderate")]),

    dict(id="clinvar_1star_PP5_supporting",
         ev=ev(consequence="missense_variant",
               clinvar_records=[cv_rec("Pathogenic", 1, 1)]),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=0, tier="VUS", met=[("PP5", "PP5_Supporting")]),

    dict(id="clinvar_0star_no_PP5",
         ev=ev(consequence="missense_variant",
               clinvar_records=[cv_rec("Pathogenic", 0, 1,
                                        rs="no assertion criteria")]),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=0, tier="VUS", met=[]),

    dict(id="clinvar_2star_BP6_moderate",
         ev=ev(consequence="missense_variant",
               clinvar_records=[cv_rec("Benign", 2, 3)]),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=0, tier="VUS", met=[("BP6", "BP6_Moderate")]),

    dict(id="clinvar_conflicting_no_PP5_BP6",
         ev=ev(consequence="missense_variant",
               clinvar_records=[cv_rec(
                   "Conflicting interpretations of pathogenicity", 1, 5)]),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=0, tier="VUS", met=[]),

    dict(id="clinvar_path_but_BA1_withholds_PP5",
         ev=ev(faf95=0.08, consequence="missense_variant",
               clinvar_records=[cv_rec("Pathogenic", 3, 4)]),
         cc={"inheritance_input": "AD"}, gene="MYH7",
         points=-8, tier="Benign", met=[("BA1", "BA1")]),

    dict(id="offpanel_frameshift_no_lof_signal",
         ev=ev(consequence="frameshift_variant", exon="4/20",
               chdgene_listed=True, chdgene_inh=["AD"]),
         cc={"inheritance_input": "AD"}, gene="GENE1",
         points=8, tier="Likely pathogenic", met=[("PVS1", "PVS1")]),

    dict(id="pathogenic_PVS1_PM2_PS2",
         ev=ev(consequence="frameshift_variant", confirmed_absent=True,
               exon="4/18", pli=0.99, chdgene_listed=True,
               chdgene_inh=["AD"]),
         cc={"inheritance_input": "AD", "denovo_status": "confirmed",
             "denovo_confirmed": True,
             "trio_status": "trio"},
         gene="MYBPC3",
         points=13, tier="Pathogenic",
         met=[("PM2", "PM2_Supporting"), ("PS2", "PS2_Strong"),
              ("PVS1", "PVS1")]),

    dict(id="pathogenic_PVS1_PP5_strong",
         ev=ev(consequence="frameshift_variant", exon="4/18", pli=0.99,
               chdgene_listed=True, chdgene_inh=["AD"],
               clinvar_records=[cv_rec("Pathogenic", 3, 4)]),
         cc={"inheritance_input": "AD"}, gene="MYBPC3",
         points=8, tier="Likely pathogenic",
         met=[("PP5", "PP5_Strong"), ("PVS1", "PVS1")]),

    dict(id="empty_inputs_safe",
         ev={}, cc={}, gene=None,
         points=0, tier="VUS", met=[]),
]


def _met_summary(criteria: list[dict]) -> list[tuple[str, str | None]]:
    """Sorted (code, criteria_strength) for every met criterion."""
    return sorted(
        (c["code"], c.get("criteria_strength"))
        for c in criteria
        if c.get("status") == "met"
    )


def _score(case: dict) -> tuple[int, str, list[tuple[str, str | None]]]:
    """Run the production deterministic pipeline for a case and return
    (points, classification, met-summary)."""
    inheritance = (case["cc"] or {}).get("inheritance_input")
    hard_coded = compute_hard_coded_criteria(
        copy.deepcopy(case["ev"]), copy.deepcopy(case["cc"]), case["gene"]
    )
    cleaned, forced = apply_cross_criterion_exclusions(
        copy.deepcopy(hard_coded), case["gene"], inheritance
    )
    points = compute_points_total(cleaned)
    tier = forced or classification_for(points)
    return points, tier, _met_summary(cleaned)


@pytest.mark.parametrize("case", GOLDEN_CASES, ids=[c["id"] for c in GOLDEN_CASES])
def test_scorer_regression(case: dict) -> None:
    """The scorer's points, tier, and met-criteria set must match the
    captured golden EXACTLY. A mismatch means deterministic scoring drifted."""
    points, tier, met = _score(case)
    assert points == case["points"], (
        f"{case['id']}: points {points} != golden {case['points']}"
    )
    assert tier == case["tier"], (
        f"{case['id']}: tier {tier!r} != golden {case['tier']!r}"
    )
    expected_met = sorted(tuple(x) for x in case["met"])
    assert met == expected_met, (
        f"{case['id']}: met {met} != golden {expected_met}"
    )


def test_all_five_tiers_covered() -> None:
    """Guard the SUITE itself: the golden battery must keep spanning all five
    ACMG tiers, so this regression net never silently narrows."""
    tiers = {c["tier"] for c in GOLDEN_CASES}
    assert tiers == {
        "Benign", "Likely benign", "VUS", "Likely pathogenic", "Pathogenic",
    }, f"tier coverage regressed: {sorted(tiers)}"


def test_all_target_criteria_covered() -> None:
    """Guard the SUITE itself: every criterion the regression net is meant to
    pin must appear (status met) in at least one golden case."""
    fired: set[str] = set()
    for case in GOLDEN_CASES:
        for code, _strength in case["met"]:
            fired.add(code)
    required = {
        "PVS1", "PM2", "BA1", "BS1", "BS2", "PP3", "BP4", "BP7", "PM4",
        "BP3", "PS2", "PM6", "PM3", "BP2", "PP5", "BP6",
    }
    missing = required - fired
    assert not missing, f"criteria no longer covered by a golden case: {missing}"
