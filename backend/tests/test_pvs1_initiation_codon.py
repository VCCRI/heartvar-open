"""PVS1 on an initiation-codon (`start_lost`) variant follows the SVI tree.

HISTORY, in two steps.

STEP 1 (fixed earlier). `start_lost` is a member of `_LOF_CONSEQUENCE_TYPES` but
had no branch of its own in `_eval_pvs1`, so it fell through to the generic
null-variant arm and scored bare `PVS1` = **+8, Very Strong**. A 6-point error in
the pathogenic direction on every start_lost in an established LoF gene.

STEP 2 (this file). The replacement capped it at a flat `PVS1_Moderate`, which is
still wrong: Abou Tayoun 2018 specifies a THREE-WAY tree, quoted verbatim —

    "the SVI Working Group generally does not recommend assigning PVS1 or
    PVS1_Strong for start loss variants. If alternative functional gene
    transcripts (i.e., found in transcript or expression databases) use an
    alternative start codon, then we recommend not applying PVS1 at any strength
    level for an initiation codon variant. If there are no alternative start
    codons in the transcript set for a gene, then we recommend applying
    PVS1_Moderate for a start loss variant if one or more pathogenic variant(s)
    have been reported 5' of the next downstream putative in-frame start codon
    (Methionine). On the other hand, if no pathogenic variant(s) occur upstream
    of the new Methionine then PVS1_Supporting should be applied."

So Moderate is the UPGRADE and Supporting is the DEFAULT. The flat cap over-called
by one point on every start-loss with no reported pathogenic variant in that
window — which is the common case, since it requires a *previously reported*
pathogenic variant in a short 5' region.

DIRECTION OF EVERY FAILURE PATH. No sequence, no ClinVar answer, no data at all →
PVS1_Supporting, never Moderate. Moderate is the branch that requires positive
evidence, so "we could not check" must not buy an extra point.

Decision, 2026-08-25: follow the literature rather than observed panel
practice. Across all 12,702 eRepo records, 36 of the 56 start-loss curations
applied PVS1 or PVS1_Strong — what this paper says not to do — so conforming
moves HeartVar AWAY from expert practice here. That is accepted deliberately and
will be cited in the manuscript.
"""
from __future__ import annotations

from backend.acmg.hard_coded import _eval_pvs1, _next_inframe_met

_LOF_EV = {
    "gnomad": {"gene": {"gnomad_constraint": {"pLI": 0.99, "oe_lof_upper": 0.20}}},
    "chdgene": {"listed": False},
    "gencc": {},
}


def _pvs1(consequence, gene="MYBPC3", vep_extra=None, start_loss=None):
    vep = {"gene_symbol": gene, "most_severe_consequence": consequence}
    vep.update(vep_extra or {})
    ev = dict(_LOF_EV)
    if start_loss is not None:
        ev["pvs1_start_loss"] = start_loss
    return _eval_pvs1(ev, vep, consequence, 0.0)


def test_next_inframe_met_finds_the_first_downstream_m():
    assert _next_inframe_met("MADEQMKL") == 6


def test_next_inframe_met_ignores_the_initiator():
    """The initiator M at residue 1 is the one being LOST — returning it would
    make the 5' window empty and silently force Supporting for the wrong reason."""
    assert _next_inframe_met("MKLST") is None


def test_next_inframe_met_handles_missing_sequence():
    for bad in (None, "", "M", "   "):
        assert _next_inframe_met(bad) is None


def test_moderate_when_a_pathogenic_variant_lies_5_prime_of_the_next_met():
    out = _pvs1("start_lost", start_loss={
        "next_inframe_met": 40, "upstream_plp": 3,
        "top_variants": [{"name": "NM_x.1(MYBPC3):c.55C>T (p.Arg19Ter)"}],
    })
    assert out["status"] == "met"
    assert out["criteria_strength"] == "PVS1_Moderate"
    assert "residue 40" in out["evidence"]


def test_supporting_when_no_pathogenic_variant_upstream_of_the_new_met():
    out = _pvs1("start_lost", start_loss={
        "next_inframe_met": 40, "upstream_plp": 0,
    })
    assert out["criteria_strength"] == "PVS1_Supporting", (
        "no reported pathogenic variant 5' of the next in-frame Met is the "
        "paper's PVS1_Supporting branch, not Moderate"
    )


def test_supporting_when_the_sequence_is_unavailable():
    """The pre-rebuild state: caches built before the UniProt `sequence` field
    was requested cannot resolve the next Met. That must not read as Moderate."""
    out = _pvs1("start_lost", start_loss={
        "next_inframe_met": None, "upstream_plp": 0,
        "reason": "UniProt sequence unavailable",
    })
    assert out["criteria_strength"] == "PVS1_Supporting"
    assert "could not be established" in out["evidence"]


def test_supporting_when_no_payload_at_all():
    out = _pvs1("start_lost")
    assert out["status"] == "met"
    assert out["criteria_strength"] == "PVS1_Supporting", (
        "absent evidence must land on the paper's default, not the upgrade"
    )


def test_a_missing_payload_can_never_produce_moderate_or_higher():
    """The whole safety property in one assertion."""
    for payload in (None, {}, {"next_inframe_met": None},
                    {"next_inframe_met": 40}, {"upstream_plp": 0},
                    {"next_inframe_met": None, "upstream_plp": 5}):
        out = _pvs1("start_lost", start_loss=payload)
        assert out["criteria_strength"] == "PVS1_Supporting", payload


def test_alternative_start_transcript_withholds_pvs1_completely():
    out = _pvs1("start_lost", start_loss={"alt_start_transcript": True})
    assert out["status"] == "not_met"
    assert out["criteria_strength"] is None


def test_branch_one_is_disclosed_as_unassessed_when_it_is_not_evaluated():
    """HeartVar does not compute branch 1, and its absence can only leave PVS1
    too HIGH — so the evidence string must say so rather than imply it passed."""
    out = _pvs1("start_lost", start_loss={"next_inframe_met": 40,
                                          "upstream_plp": 0})
    assert "not assessed" in out["evidence"].lower()
    assert "alternative" in out["evidence"].lower()


def test_start_lost_never_scores_very_strong_or_strong():
    for payload in (None, {"next_inframe_met": 40, "upstream_plp": 99}):
        out = _pvs1("start_lost", start_loss=payload)
        assert out["criteria_strength"] not in (
            "PVS1", "PVS1_VeryStrong", "PVS1_Strong",
        ), "the paper rules out PVS1 and PVS1_Strong for start-loss outright"


def test_other_null_consequences_keep_full_pvs1():
    for csq in ("stop_gained", "frameshift_variant", "transcript_ablation"):
        out = _pvs1(csq)
        assert out["status"] == "met", csq
        assert out["criteria_strength"] == "PVS1", (
            f"{csq} must still score full PVS1; got {out['criteria_strength']}"
        )


def test_a_compound_consequence_with_a_truncation_is_not_capped():
    """Guard on the branch condition. If VEP hands back a consequence string
    naming both start_lost and a truncating event, the truncation argument is
    the stronger one and must win — the arm keys on start_lost being the ONLY
    LoF term present."""
    out = _pvs1("start_lost,stop_gained")
    assert out["criteria_strength"] == "PVS1"


def test_start_lost_still_needs_an_established_lof_mechanism():
    ev = {"gnomad": {"gene": {"gnomad_constraint": {"pLI": 0.01,
                                                    "oe_lof_upper": 1.4}}},
          "chdgene": {"listed": False}, "gencc": {}}
    out = _eval_pvs1(ev, {"gene_symbol": "MYBPC3"}, "start_lost", 0.0)
    assert out["status"] == "not_met"
