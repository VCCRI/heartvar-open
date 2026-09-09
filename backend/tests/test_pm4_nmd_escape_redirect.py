"""PM4 for NMD-escaping truncating variants, where the VCEP redirects to it.

Seven cardiomyopathy genes withdraw PVS1 and point at PM4 by name. TNNT2's PVS1
record, verbatim from the harvest:

    "Not currently applicable to TNNT2. See PM4 for truncating variants that do
     NOT undergo NMD."

The same sentence appears on ACTC1, MYL2, MYL3, TNNI3 and TPM1. MYH7 carries no
cross-reference at all ("Not applicable for MYH7."), and reaches the route
through the PM4 side, whose Moderate rung states the general rule:

    "For genes where PVS1 is not applicable (i.e., where there is no evidence
     that pLOF variants cause disease), consider using this rule at MODERATE or
     SUPPORTING strength for truncating variants that do NOT undergo nonsense
     mediated decay (NMD)."

MYBPC3 is DELIBERATELY NOT one of them, and this file originally had it the
wrong way round in both directions. MYBPC3 is the one Cardiomyopathy gene where
LOF is an established mechanism, so its PVS1 stays `applicable`; the "see PM4"
phrase in its record comes from shared boilerplate inside the PVS1 strength
text. Letting the route fire there would have scored PVS1 (downgraded for NMD
escape) AND PM4 for one truncation.

HeartVar's PM4 implemented only the two ACMG-2015 arms — in-frame indels and
stop-loss — so the redirect was dead: PM4 never fired on the case an expert
applies it to. TNNT2 c.890G>A (p.Trp297Ter) is that shape — a last-exon
nonsense variant that escapes NMD, where PM4 is applicable and HeartVar was
not applying it.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from backend.acmg.hard_coded import _eval_pm4_bp3, _pm4_nmd_redirect_gene


def _pm4(entries):
    return next(e for e in entries if e["code"] == "PM4")


REDIRECT_GENES = ("ACTC1", "MYH7", "MYL2", "MYL3", "TNNI3", "TNNT2", "TPM1")


def test_the_seven_redirect_genes_are_recognised():
    for gene in REDIRECT_GENES:
        assert _pm4_nmd_redirect_gene(gene) is True, gene


def test_genes_without_the_redirect_are_not_swept_in():
    for gene in ("PTPN11", "FBN1", "KCNQ1", "BMPR2", "MYBPC3"):
        assert _pm4_nmd_redirect_gene(gene) is False, gene


def test_the_redirect_set_is_exactly_the_pvs1_withdrawn_cardiomyopathy_genes():
    """A whole-registry pin. The route may only fire where the VCEP has
    withdrawn PVS1, so no gene can score PVS1 and PM4 for one truncation."""
    import json
    from backend.acmg.hard_coded import _criterion_applicable
    path = Path(__file__).resolve().parents[1] / "data" / "vcep_criteria_spec.json"
    genes = json.loads(path.read_text())["genes"]
    fires = sorted(g for g in genes if _pm4_nmd_redirect_gene(g))
    assert fires == sorted(REDIRECT_GENES), fires
    for gene in fires:
        assert _criterion_applicable(gene, "PVS1") is False, gene


def test_pvs1_and_pm4_can_never_both_take_one_truncation():
    import json
    from backend.acmg.hard_coded import _criterion_applicable
    path = Path(__file__).resolve().parents[1] / "data" / "vcep_criteria_spec.json"
    for gene in json.loads(path.read_text())["genes"]:
        assert not (_pm4_nmd_redirect_gene(gene)
                    and _criterion_applicable(gene, "PVS1")), gene


def test_junk_genes_never_raise():
    for bad in ("", "NOTAGENE", None, "  "):
        assert _pm4_nmd_redirect_gene(bad) is False


def test_nmd_escaping_nonsense_earns_pm4_on_a_redirect_gene():
    """TNNT2 c.890G>A p.Trp297Ter — last exon, escapes NMD, expert applied PM4."""
    vep = {"ok": True, "nmd_escape": True}
    got = _pm4(_eval_pm4_bp3(vep, "stop_gained", {}, gene="TNNT2"))
    assert got["status"] == "met", got
    assert got["criteria_strength"] == "PM4_Moderate", got["criteria_strength"]
    assert "NMD" in got["evidence"]


def test_frameshift_escaping_nmd_also_qualifies():
    vep = {"ok": True, "nmd_escape": True}
    got = _pm4(_eval_pm4_bp3(vep, "frameshift_variant", {}, gene="TNNI3"))
    assert got["status"] == "met", got


def test_a_truncating_variant_that_UNDERGOES_nmd_gets_nothing():
    """The redirect is explicitly for variants that do NOT undergo NMD. One that
    does is PVS1's business on genes where PVS1 applies, and nothing here."""
    vep = {"ok": True, "nmd_escape": False}
    got = _pm4(_eval_pm4_bp3(vep, "stop_gained", {}, gene="TNNT2"))
    assert got["status"] == "not_met", got


def test_the_route_is_scoped_to_redirect_genes():
    vep = {"ok": True, "nmd_escape": True}
    for gene in ("PTPN11", "FBN1", "MYBPC3"):
        got = _pm4(_eval_pm4_bp3(vep, "stop_gained", {}, gene=gene))
        assert got["status"] == "not_met", (gene, got)


def test_missing_nmd_signal_fails_closed():
    """No nmd_escape key means we do not know, and an unknown must not assert
    a pathogenic criterion."""
    got = _pm4(_eval_pm4_bp3({"ok": True}, "stop_gained", {}, gene="TNNT2"))
    assert got["status"] == "not_met", got


def test_stop_loss_still_earns_pm4_everywhere():
    got = _pm4(_eval_pm4_bp3({"ok": True}, "stop_lost", {}, gene="PTPN11"))
    assert got["status"] == "met" and got["criteria_strength"] == "PM4_Moderate"


def test_inframe_indel_still_earns_pm4():
    got = _pm4(_eval_pm4_bp3({"ok": True}, "inframe_deletion", {}, gene="MYBPC3"))
    assert got["status"] == "met", got


def test_missense_still_gets_nothing():
    got = _pm4(_eval_pm4_bp3({"ok": True, "nmd_escape": True},
                             "missense_variant", {}, gene="TNNT2"))
    assert got["status"] == "not_met", got
