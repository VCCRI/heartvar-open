"""Gate-correctness audit for the deterministic ACMG criteria.

Companion to test_bp2_in_trans_spec.py. Both pin the same class of defect: a
criterion whose *gating inputs* admit a clinical state the criterion's ACMG
definition excludes. The UI renders the correct ACMG definition next to each
code, so a mis-gated criterion is invisible to the curator until they question
an individual call.

Covers the two gates found mis-specified in the 2026-08-14 audit:

PS2 on a duo
    ACMG/AMP 2015: PS2 = "De novo (both maternity AND paternity confirmed) in a
    patient with the disease and no family history." A *duo* tests ONE parent,
    so it cannot confirm both. `denovo_confirmed` admitted
    trio_status in ("trio", "duo"), so duo + "confirmed de novo" — both
    selectable in the UI — scored PS2_Strong (+4) and printed the evidence
    string "confirmed by parental testing (full trio)", which is false. The
    correct code for de novo asserted without both parents confirmed is PM6
    (+2), which the same inputs were also denying.

BS2 on recessive obligate carriers
    ACMG/AMP 2015: BS2 = "Observed in a healthy adult individual for a
    recessive (HOMOZYGOUS), dominant (HETEROZYGOUS), or X-linked (HEMIZYGOUS)
    disorder, with full penetrance expected at an early age." The genotype must
    be the disease-causing one for that inheritance mode. The curator path fired
    BS2_Strong (-4) off `unaffected_obs_count` alone, ignoring inheritance — so
    for a recessive disorder the proband's own unaffected heterozygous carrier
    parents, who are obligate and entirely expected, scored -4 and drove the
    call to Likely benign.
"""
from __future__ import annotations

from backend.acmg.hard_coded import compute_hard_coded_criteria
from backend.acmg.tiers import compute_points_total


def _crit(cc: dict, gene: str = "ZZZZ9") -> dict:
    out = compute_hard_coded_criteria({}, cc, gene=gene)
    return {c["code"]: c for c in out}


def _denovo_cc(trio: str, status: str = "confirmed") -> dict:
    """Mirror app._build_clinical_context's denovo_confirmed derivation."""
    from backend.app import _build_clinical_context
    from backend.models import CurationRequest

    req = CurationRequest(
        gene="MYH7", hgvs_c="c.1A>G", phenotype="HCM",
        inheritance_input="AD", zygosity="het",
        trio_status=trio, denovo_status=status,
    )
    return _build_clinical_context(req, {})


def test_duo_cannot_confirm_de_novo_so_ps2_must_not_fire():
    """One parent tested cannot confirm both maternity and paternity."""
    crit = _crit(_denovo_cc("duo"), gene="MYH7")
    assert crit["PS2"]["status"] == "not_met", (
        "PS2_Strong fired off a duo — ACMG requires BOTH parents confirmed: "
        f"{crit['PS2'].get('evidence')!r}"
    )


def test_duo_confirmed_de_novo_earns_NO_de_novo_criterion():
    """🔴 CORRECTED 2026-09-07. This asserted PM6_Moderate, on the reasoning
    that "the evidence isn't lost — it downgrades to the assumed-de-novo code".

    That was half right. It correctly stopped a duo earning PS2, but PM6 is not
    a safe landing place either: PM6 is "ASSUMED de novo, but without
    confirmation of paternity and maternity", and both codes presuppose the
    variant is ABSENT FROM BOTH PARENTS. The confirmed/assumed axis is about
    whether the tested people are genetically verified as the parents (ClinGen
    SVI de novo v1.1: "confirmed parental relationships versus assumed parental
    relationships status"), NOT about how many parents were tested.

    With one parent untested the variant may simply have been inherited from
    the untested parent, so there is no de-novo observation to discount. ACMG
    PS2, verbatim: "Confirmation of paternity only is insufficient."

    The curator is told this rather than left with a silent zero."""
    crit = _crit(_denovo_cc("duo"), gene="MYH7")
    assert crit["PS2"]["status"] == "not_met"
    assert crit["PM6"]["status"] == "not_met"
    ev = crit["PM6"]["evidence"].lower()
    assert "one parent" in ev and "both parents" in ev, crit["PM6"]["evidence"]


def test_full_trio_confirmed_still_earns_ps2():
    """The legitimate PS2 path is untouched."""
    crit = _crit(_denovo_cc("trio"), gene="MYH7")
    assert crit["PS2"]["status"] == "met"
    assert crit["PS2"]["criteria_strength"] == "PS2_Strong"
    assert crit["PM6"]["status"] == "not_met"


def test_ps2_evidence_string_never_claims_full_trio_for_a_duo():
    """The old message asserted "(full trio)" regardless of trio_status."""
    for trio in ("duo", "trio"):
        crit = _crit(_denovo_cc(trio), gene="MYH7")
        for code in ("PS2", "PM6"):
            entry = crit[code]
            if entry["status"] == "met" and "full trio" in entry["evidence"]:
                assert trio == "trio", (
                    f"{code} claimed 'full trio' for trio_status={trio!r}"
                )


def test_duo_scores_nothing_and_a_trio_scores_four():
    """Corrected 2026-09-07 from "duo scores two, not four". A duo is not a
    de-novo observation at all, so it contributes nothing — see
    test_duo_confirmed_de_novo_earns_NO_de_novo_criterion for the reasoning.
    Direction of the correction is DOWN, removing +2 that was never earned."""
    pts = compute_points_total(list(_crit(_denovo_cc("duo"), gene="MYH7").values()))
    trio_pts = compute_points_total(
        list(_crit(_denovo_cc("trio"), gene="MYH7").values())
    )
    assert pts == 0, f"duo should contribute nothing, got {pts}"
    assert trio_pts == 4, f"trio should contribute PS2 (+4), got {trio_pts}"


def test_unconfirmed_trio_earns_pm6_and_unconfirmed_duo_does_not():
    """Corrected with the test above: a stated TRIO with unverified parentage
    is exactly PM6; a DUO is not a de-novo observation at all."""
    crit = _crit(_denovo_cc("trio", status="unconfirmed"), gene="MYH7")
    assert crit["PM6"]["status"] == "met"
    assert crit["PS2"]["status"] == "not_met"

    crit = _crit(_denovo_cc("duo", status="unconfirmed"), gene="MYH7")
    assert crit["PM6"]["status"] == "not_met"
    assert crit["PS2"]["status"] == "not_met"


def test_pm4_fires_for_stop_loss_even_inside_a_repeat_region():
    from backend.acmg.hard_coded import _eval_pm4_bp3

    d = {c["code"]: c for c in _eval_pm4_bp3({"in_repeat_region": True}, "stop_lost")}
    assert d["PM4"]["status"] == "met", (
        "stop-loss denied PM4 by the repeat-region gate: "
        f"{d['PM4'].get('evidence')!r}"
    )
    assert d["BP3"]["status"] == "not_met", "BP3 is for in-frame indels only"


def test_pm4_not_met_message_never_contradicts_the_consequence():
    from backend.acmg.hard_coded import _eval_pm4_bp3

    for cons in ("stop_lost", "inframe_deletion", "missense_variant"):
        for rep in (True, False):
            d = {c["code"]: c for c in _eval_pm4_bp3({"in_repeat_region": rep}, cons)}
            pm4 = d["PM4"]
            if pm4["status"] == "not_met" and "Not an in-frame indel or stop-loss" in pm4["evidence"]:
                assert "inframe" not in cons and "stop_lost" not in cons, (
                    f"PM4 said {cons!r} is not an in-frame indel or stop-loss"
                )


def test_bp3_still_owns_inframe_indels_in_functionless_repeats():
    """BP3 keeps the in-frame-indel arm — see test_bp3_repeat_function.py for
    the full "without a known function" gate this now runs through."""
    from backend.acmg.hard_coded import _eval_pm4_bp3

    ev = {"uniprot": {"ok": True, "found": True, "domains": [], "features": [
        {"type": "Repeat", "start": 90, "end": 120, "description": "Gln-rich"},
    ]}}
    d = {c["code"]: c for c in
         _eval_pm4_bp3({"protein_start": 100}, "inframe_deletion", ev)}
    assert d["BP3"]["status"] == "met"
    assert d["PM4"]["status"] == "not_met"


def test_legacy_insilico_consensus_never_claims_a_strength_it_cannot_reach():
    """The consensus pool holds 4 tools; the old code escalated at >= 5.

    Dead code that reads as an escalation path. Pin the pool size so that if a
    fifth predictor is ever added, whoever adds it must decide deliberately
    whether the legacy consensus escalates (SVI says prefer ONE calibrated tool
    over a correlated-tool vote, and the cardiac VCEPs cap in-silico at
    Supporting — so the answer is probably still no).
    """
    from backend.acmg.hard_coded import _eval_pp3_bp4

    vep = {"revel_score": 0.99, "cadd_phred": 40}
    ev = {"alphamissense": {"am_pathogenicity": 0.99},
          "spliceai": {"ok": True}}
    out, _, _ = _eval_pp3_bp4(ev, vep, "inframe_deletion", 0.99)
    pp3 = {c["code"]: c for c in out}["PP3"]
    assert pp3["status"] == "met", "all four tools agreed yet PP3 did not fire"
    assert pp3["criteria_strength"] == "PP3_Supporting", (
        "the legacy consensus escalated above Supporting — SVI prefers one "
        "calibrated tool over a correlated-tool vote, and the cardiac VCEPs cap "
        f"in-silico at Supporting; got {pp3['criteria_strength']}"
    )


def test_bs2_gnomad_homozygote_path_unaffected_by_the_inheritance_gate():
    """The population hom/hemi path is a different, already-zygosity-aware arm."""
    ev = {"gnomad": {"ok": True, "variant_found": True,
                     "variant": {"joint": {"ac_hom": 12}}}}
    out = compute_hard_coded_criteria(
        ev, {"inheritance_input": "AR", "zygosity": "hom"}, gene="ZZZZ9",
    )
    bs2 = {c["code"]: c for c in out}["BS2"]
    if bs2["status"] == "not_met":
        assert "carrier" not in bs2["evidence"].lower()
