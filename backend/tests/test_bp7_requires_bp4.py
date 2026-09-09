"""BP7 requires BP4 on the genes whose VCEP mandates the conjunction.

Most of the VCEP-spec'd genes state it. The RASopathy specs use:

    "This rule is also applicable for intronic positions (except canonical
     splice sites) or non-coding variants and should be used in conjunction
     with BP4."

and BMPR2 makes it a precondition outright:

    "Applicable after assignment of BP4 for no adverse splicing predictions,
     and inclusive of exonic and intronic variants."

⚠ The 8 Cardiomyopathy genes are DELIBERATELY EXCLUDED. Their wording is
permissive, not mandatory — "Rule CAN be combined with BP4 to make a variant
likely benign per Richards et al." — so requiring BP4 there would invent a rule
the panel did not write.

Most BP7 false positives are on mandate genes with no BP4, and removing them
moves no variant's tier in either arm (BP7 is worth -1 and those variants sit
well inside Benign). So this is a criterion-precision fix with no tier cost.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from backend.acmg.hard_coded import _bp7_requires_bp4, apply_cross_criterion_exclusions


def _crit(code, status="met", strength=None, evidence="x"):
    return {"code": code, "status": status,
            "criteria_strength": strength or code, "evidence": evidence}


def _run(codes, gene):
    cleaned, _ = apply_cross_criterion_exclusions(
        [_crit(c) for c in codes], gene=gene, inheritance="AD")
    return {c["code"]: c for c in cleaned}


def test_the_rasopathy_genes_and_bmpr2_mandate_bp4():
    for gene in ("BRAF", "KRAS", "MAP2K1", "MAP2K2", "PTPN11", "RAF1", "SHOC2",
                 "SOS1", "SOS2", "HRAS", "NRAS", "RIT1", "RRAS2", "LZTR1",
                 "MRAS", "PPP1CB", "BMPR2"):
        assert _bp7_requires_bp4(gene) is True, gene


def test_the_cardiomyopathy_genes_only_permit_it():
    """"Rule CAN be combined with BP4" is permission, not a requirement."""
    for gene in ("ACTC1", "MYBPC3", "MYH7", "MYL2", "MYL3", "TNNI3", "TNNT2",
                 "TPM1"):
        assert _bp7_requires_bp4(gene) is False, gene


def test_junk_and_unknown_genes_fail_open():
    for bad in ("", "NOTAGENE", None, "  "):
        assert _bp7_requires_bp4(bad) is False


def test_bp7_without_bp4_is_demoted_on_a_mandate_gene():
    got = _run(["BP7"], "MAP2K1")
    assert got["BP7"]["status"] == "not_met", got["BP7"]
    assert "BP4" in got["BP7"]["evidence"]


def test_bp7_with_bp4_survives_on_a_mandate_gene():
    got = _run(["BP7", "BP4"], "MAP2K1")
    assert got["BP7"]["status"] == "met", got["BP7"]


def test_bp7_alone_survives_on_a_permissive_gene():
    got = _run(["BP7"], "MYH7")
    assert got["BP7"]["status"] == "met", got["BP7"]


def test_bp7_alone_survives_on_a_gene_outside_the_panel():
    got = _run(["BP7"], "GATA4")
    assert got["BP7"]["status"] == "met", got["BP7"]


def test_bp4_is_never_touched_by_this_rule():
    got = _run(["BP7", "BP4"], "MAP2K1")
    assert got["BP4"]["status"] == "met"
    got = _run(["BP4"], "MAP2K1")
    assert got["BP4"]["status"] == "met"
