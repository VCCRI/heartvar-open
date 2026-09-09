"""Unit tests for backend.app.mechanism_consistency_flag.

The deterministic variant-type ↔ disease-mechanism note must fire for exactly
ONE curated, ~zero-false-positive direction: a canonical loss-of-function (null)
variant in a gene whose ClinGen VCEP marks PVS1 "Not Applicable" (loss of
function is not the established disease mechanism). It must stay SILENT for:
  * missense / non-null consequences (even in those genes),
  * established-LoF genes that have pathogenic missense (KCNQ1, MYBPC3),
  * genes with no VCEP applicability entry (default-open — no warning without
    explicit expert evidence).

This guards against the dangerous false-warning class the feature was designed
to avoid (e.g. mis-flagging a frameshift in a haploinsufficient gene, or any
missense as "inconsistent").

No pytest dependency — runnable with pytest or directly.
"""
from __future__ import annotations

from backend.acmg.hard_coded import (
    _GENE_DOSAGE,
    _criterion_applicable,
    gene_mechanism,
)
from backend.app import mechanism_consistency_flag


def _ev(consequence: str, mis_z: float | None = None) -> dict:
    ev: dict = {"vep": {"most_severe_consequence": consequence}}
    if mis_z is not None:
        ev["gnomad"] = {"gene": {"gnomad_constraint": {"mis_z": mis_z}}}
    return ev


def test_fires_for_nonsense_in_rasopathy_gof_gene():
    flag = mechanism_consistency_flag(_ev("stop_gained"), "HRAS")
    assert flag is not None
    assert flag["kind"] == "lof_in_non_lof_gene"
    assert flag["severity"] == "warning"
    assert flag["gene"] == "HRAS"
    assert "PVS1" in flag["detail"]
    assert "loss of function" in flag["detail"].lower()


def test_fires_for_frameshift_in_dominant_negative_sarcomeric_gene():
    flag = mechanism_consistency_flag(_ev("frameshift_variant"), "MYH7")
    assert flag is not None
    assert flag["gene"] == "MYH7"


def test_fires_for_canonical_splice_in_not_applicable_gene():
    assert mechanism_consistency_flag(_ev("splice_donor_variant"), "TPM1") is not None
    assert mechanism_consistency_flag(_ev("splice_acceptor_variant"), "TNNT2") is not None


def test_gene_symbol_case_insensitive():
    assert mechanism_consistency_flag(_ev("stop_gained"), "hras") is not None


def test_silent_for_missense_even_in_not_applicable_gene():
    assert mechanism_consistency_flag(_ev("missense_variant"), "HRAS") is None
    assert mechanism_consistency_flag(_ev("missense_variant"), "MYH7") is None


def test_silent_for_synonymous_and_intronic():
    assert mechanism_consistency_flag(_ev("synonymous_variant"), "HRAS") is None
    assert mechanism_consistency_flag(_ev("intron_variant"), "MYH7") is None


def test_silent_for_lof_gene_with_pvs1_applicable():
    assert mechanism_consistency_flag(_ev("stop_gained"), "KCNQ1") is None


def test_silent_for_gene_absent_from_vcep_table():
    assert mechanism_consistency_flag(_ev("frameshift_variant"), "MYBPC3") is None
    assert mechanism_consistency_flag(_ev("stop_gained"), "TTN") is None


def test_mixed_mechanism_gene_gets_info_note_not_warning():
    for gene in ("SCN5A", "PTPN11"):
        for csq in ("stop_gained", "missense_variant"):
            flag = mechanism_consistency_flag(_ev(csq), gene)
            assert flag is not None, (gene, csq)
            assert flag["kind"] == "mixed_mechanism"
            assert flag["severity"] == "info"
    assert mechanism_consistency_flag(_ev("intron_variant"), "SCN5A") is None


def test_curated_only_dn_gene_null_variant_is_info_not_warning():
    # ClinGen dosage record — a null variant gets the softer, attributed note.
    flag = mechanism_consistency_flag(_ev("stop_gained"), "MYH6")
    assert flag is not None
    assert flag["kind"] == "lof_in_dn_gene"
    assert flag["severity"] == "info"


def test_soft_note_for_missense_in_haploinsufficiency_gene():
    flag = mechanism_consistency_flag(_ev("missense_variant", mis_z=-1.4), "GATA4")
    assert flag is not None
    assert flag["kind"] == "missense_in_lof_gene"
    assert flag["severity"] == "info"
    assert "not a benign call" in flag["detail"].lower()


def test_soft_note_suppressed_when_gene_is_missense_constrained():
    assert mechanism_consistency_flag(_ev("missense_variant", mis_z=8.2), "GATA4") is None


def test_soft_note_suppressed_when_vcep_keeps_pp2_applicable():
    """The pp2_applicable suppression inside direction (3).

    SENTINEL CHOICE IS LOAD-BEARING and the assertion is worthless without it.
    Direction (3) is only entered for a gene with mechanism
    "haploinsufficiency", confidence "established" and ClinGen HI == 3; for
    anything else mechanism_consistency_flag returns None several branches
    earlier, so `is None` would pass with the suppression deleted outright.
    TPM1 is exactly that trap — ClinGen HI = 0 and mechanism gof_or_dn — so the
    three preconditions are asserted here rather than assumed.

    KCNQ1 used to be this sentinel and no longer can be: GN112 withdraws PP2.
    See test_soft_note_fires_for_kcnq1_because_gn112_withdraws_pp2, which is the
    negative control on the same branch.
    """
    for gene in ("BMPR2", "FBN1"):
        mech = gene_mechanism(_ev("missense_variant"), gene)
        assert mech["mechanism"] == "haploinsufficiency", (gene, mech)
        assert mech["confidence"] == "established", (gene, mech)
        assert str(_GENE_DOSAGE[gene]["hi_score"]) == "3", gene
        assert _criterion_applicable(gene, "PP2") is True, gene
        assert mechanism_consistency_flag(
            _ev("missense_variant"), gene) is None, gene


def test_soft_note_fires_for_kcnq1_because_gn112_withdraws_pp2():
    """Negative control for the test above, on the SAME branch and with the
    same preconditions — which is what makes the pair non-vacuous.

    KCNQ1 is ClinGen HI=3 / established haploinsufficiency like BMPR2 and FBN1,
    but GN112 withdraws PP2: "Not applicable due to presence of benign
    variation throughout the KCNQ1 gene (since the missense constraint Z-score
    in gnomAD is 1.83, lower than 3)." So the suppression must NOT apply and the
    soft mechanism-fit note fires. Informational only — this function never
    changes a score."""
    mech = gene_mechanism(_ev("missense_variant"), "KCNQ1")
    assert mech["mechanism"] == "haploinsufficiency", mech
    assert mech["confidence"] == "established", mech
    assert str(_GENE_DOSAGE["KCNQ1"]["hi_score"]) == "3"
    assert _criterion_applicable("KCNQ1", "PP2") is False
    flag = mechanism_consistency_flag(_ev("missense_variant"), "KCNQ1")
    assert flag is not None
    assert flag["kind"] == "missense_in_lof_gene"
    assert flag["severity"] == "info"
    assert "not a benign call" in flag["detail"].lower()


def test_null_variant_in_haploinsufficiency_gene_is_consistent():
    assert mechanism_consistency_flag(_ev("stop_gained"), "GATA4") is None


def test_silent_for_undetermined_gene():
    assert mechanism_consistency_flag(_ev("missense_variant"), "RYR2") is None
    assert mechanism_consistency_flag(_ev("stop_gained"), "RYR2") is None


def test_silent_for_missing_gene_or_consequence():
    assert mechanism_consistency_flag(_ev("stop_gained"), "") is None
    assert mechanism_consistency_flag(_ev(""), "HRAS") is None
    assert mechanism_consistency_flag({}, "HRAS") is None
    assert mechanism_consistency_flag({"vep": {}}, "HRAS") is None


if __name__ == "__main__":  # pragma: no cover
    import sys

    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"ok  {fn.__name__}")
    print(f"\n{len(fns)} tests passed")
    sys.exit(0)
