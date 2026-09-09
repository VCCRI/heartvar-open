"""FBN1 (GN022) and KCNQ1 (GN112) PM1, implemented from their published specs.

WHY THESE TWO. Of the 25 genes with a CSpec, 23 already had PM1 handled correctly
— 6 suppressed outright, 17 by numeric hotspot ranges. FBN1 and KCNQ1 were the
two where a published definition exists and HeartVar substituted its own
domain-P/LP-density heuristic instead. For the other 652 panel genes there IS no
spec, so the heuristic is a house rule there and legitimate.

KCNQ1's rule is a range plus a precondition. FBN1's is an enumeration over the
protein's domain architecture that publishes TWO strengths — which is exactly
what the old hand-transcribed table could not represent, since FBN1's PM1 arrived
there as one merged cell joining the Strong rule to the Moderate one. Implementing
from that copy would have scored a Strong rule at Moderate.

The FBN1 fixture is synthetic, not fetched: cysteines are placed at known indices
so the inter-cysteine interval each test residue falls in is arithmetic rather
than something to look up. UniProt's real naming is reproduced exactly, though —
"EGF-like N; calcium-binding" for the 43 cbEGF domains, plain "EGF-like N" for the
4 others, "TB N" for the 9 — because that naming is what makes the classes
decidable at all.
"""
from __future__ import annotations

import pytest

from backend.acmg.hard_coded import (
    _cys_index_within,
    _fbn1_domain_class,
    _gate_criteria_applicability,
    _hgvsp_ref_alt_aa,
    _parse_pm1_hotspots,
    _pm1_fbn1,
    _crit_text,
)

_CB_CYS = [102, 108, 114, 120, 126, 132]
_DOMAINS = [
    {"name": "EGF-like 4; calcium-binding", "start": 40, "end": 90},
    {"name": "EGF-like 5; calcium-binding", "start": 100, "end": 140},
    {"name": "EGF-like 3", "start": 200, "end": 240},
    {"name": "TB 1", "start": 300, "end": 340},
]


def _sequence() -> str:
    seq = ["A"] * 400
    for i in _CB_CYS:
        seq[i - 1] = "C"
    for i in (42, 48, 54, 60):
        seq[i - 1] = "C"
    seq[110 - 1] = "G"
    seq[117 - 1] = "G"
    seq[112 - 1] = "N"
    seq[202 - 1] = "C"
    seq[302 - 1] = "C"
    seq[210 - 1] = "R"
    return "".join(seq)


def _ev(domains=None, sequence=None) -> dict:
    return {"uniprot": {"domains": _DOMAINS if domains is None else domains,
                        "sequence": _sequence() if sequence is None else sequence}}


def _p(ref3: str, pos: int, alt3: str) -> str:
    return f"NP_000129.3:p.{ref3}{pos}{alt3}"


def test_uniprots_naming_distinguishes_cbegf_from_plain_egf():
    """The whole implementation rests on this. UniProt describes FBN1's 43
    calcium-binding EGF-like domains as "EGF-like N; calcium-binding" and its 4
    others as plain "EGF-like N"; the spec gives them different STRENGTHS, so
    conflating them would score a Strong rule at Moderate."""
    assert _fbn1_domain_class("EGF-like 5; calcium-binding") == "cbEGF"
    assert _fbn1_domain_class("EGF-like 3") == "EGF-like"
    assert _fbn1_domain_class("TB 1") == "TB"
    assert _fbn1_domain_class("Hybrid 2") == "hybrid"
    assert _fbn1_domain_class("Fibrillin repeat") is None
    assert _fbn1_domain_class(None) is None


def test_inter_cysteine_interval_arithmetic():
    seq = _sequence()
    assert _cys_index_within(seq, 100, 140, 110) == 2
    assert _cys_index_within(seq, 100, 140, 117) == 3
    assert _cys_index_within(seq, 100, 140, 102) is None
    assert _cys_index_within(None, 100, 140, 110) is None
    assert _cys_index_within(seq, 100, 140, 999) is None


def test_cysteine_in_a_cbegf_domain_is_strong():
    verdict, why = _pm1_fbn1(_ev(), _p("Cys", 102, "Arg"))
    assert verdict == "Strong", why
    assert "cbEGF" in why


def test_cysteine_in_a_plain_egf_domain_is_moderate():
    verdict, why = _pm1_fbn1(_ev(), _p("Cys", 202, "Arg"))
    assert verdict == "Moderate", why


def test_cysteine_in_a_tb_domain_is_moderate():
    verdict, why = _pm1_fbn1(_ev(), _p("Cys", 302, "Arg"))
    assert verdict == "Moderate", why


def test_a_cysteine_creating_variant_is_moderate_anywhere():
    """The spec attaches no domain qualifier to "Cys-creating variants"."""
    verdict, why = _pm1_fbn1(_ev(), _p("Arg", 380, "Cys"))
    assert verdict == "Moderate", why
    assert "creating" in why


def test_critical_glycine_between_cys2_and_cys3_is_moderate():
    verdict, why = _pm1_fbn1(_ev(), _p("Gly", 110, "Ser"))
    assert verdict == "Moderate", why
    assert "Cys2 and Cys3" in why


def test_glycine_between_cys3_and_cys4_needs_an_upstream_cbegf_domain():
    """The spec conditions this one on an upstream cbEGF domain existing, so the
    condition has to be checked rather than assumed."""
    verdict, why = _pm1_fbn1(_ev(), _p("Gly", 117, "Ser"))
    assert verdict == "Moderate", why
    assert "upstream" in why

    only_this = [d for d in _DOMAINS if d["start"] != 40]
    verdict2, _ = _pm1_fbn1(_ev(domains=only_this), _p("Gly", 117, "Ser"))
    assert verdict2 is None


def test_g_to_a_is_barred_even_at_a_critical_glycine():
    """Residue 110 is a critical Gly between Cys2 and Cys3, so the Moderate rule
    matches — and the caveat still forbids PM1. That precedence IS the caveat."""
    assert _pm1_fbn1(_ev(), _p("Gly", 110, "Ser"))[0] == "Moderate"
    verdict, why = _pm1_fbn1(_ev(), _p("Gly", 110, "Ala"))
    assert verdict == "barred", why
    assert "tolerated" in why


def test_n_to_s_in_a_cbegf_domain_is_barred():
    verdict, why = _pm1_fbn1(_ev(), _p("Asn", 112, "Ser"))
    assert verdict == "barred", why


def test_n_to_s_outside_cbegf_is_not_barred_by_that_caveat():
    """The caveat is written about the cbEGF consensus sequence, so it must not
    leak into other domains."""
    verdict, _ = _pm1_fbn1(_ev(), _p("Asn", 210, "Ser"))
    assert verdict != "barred"


def test_a_non_cysteine_in_a_tb_domain_is_not_covered():
    """GN022 scopes PM1 in TB / plain-EGF domains to cysteines, so anything else
    there is outside the enumeration."""
    verdict, why = _pm1_fbn1(_ev(), _p("Arg", 310, "Trp"))
    assert verdict == "no_rule", why


def test_outside_every_domain_is_not_covered():
    verdict, why = _pm1_fbn1(_ev(), _p("Arg", 5, "Trp"))
    assert verdict == "no_rule", why


def test_undecidable_inside_cbegf_is_left_alone_not_suppressed():
    """FAIL-OPEN, deliberately. Two Moderate rules are NOT implemented — the
    (D/N)-X-(D/N)-(E/Q) consensus substitution and the invariant
    calcium-binding/hydroxylation residues — and both live inside cbEGF domains.
    Suppressing an unmatched cbEGF residue would delete PM1 for variants those
    rules cover, so an unmatched residue there returns None."""
    verdict, why = _pm1_fbn1(_ev(), _p("Arg", 105, "Trp"))
    assert verdict is None, why


def test_missing_domain_annotation_is_undecidable_not_uncovered():
    """No UniProt annotation is not evidence that the enumeration excludes the
    variant — a hybrid-domain cysteine would look identical."""
    assert _pm1_fbn1(_ev(domains=[]), _p("Cys", 102, "Arg"))[0] is None


def test_a_frameshift_hgvsp_is_undecidable():
    assert _hgvsp_ref_alt_aa("NP_000129.3:p.Arg545fs") == (None, None)
    assert _pm1_fbn1(_ev(), "NP_000129.3:p.Arg545fs")[0] is None


def _criteria(*codes):
    return [{"code": c, "status": "met", "criteria_strength": f"{c}_Moderate",
             "evidence": "from the model"} for c in codes]


def _vep(hgvsp: str) -> dict:
    ev = _ev()
    ev["vep"] = {"ok": True, "most_severe_consequence": "missense_variant",
                 "hgvsp": hgvsp}
    return ev


def _by_code(out, code):
    return [c for c in out if c["code"] == code][0]


def test_gate_sets_pm1_strong_and_bars_pm5_and_ps1():
    """GN022's Strong rule carries "PM5/PS1 should not be used when this argument
    applies" — the cysteine's importance IS the argument, so a different
    substitution at the same cysteine is the same evidence counted twice."""
    out = _gate_criteria_applicability(
        _criteria("PM1", "PM5", "PS1"), "FBN1", _vep(_p("Cys", 102, "Arg")))
    assert _by_code(out, "PM1")["criteria_strength"] == "PM1_Strong"
    assert _by_code(out, "PM5")["status"] == "not_met"
    assert _by_code(out, "PS1")["status"] == "not_met"


def test_gate_leaves_pm5_alone_when_pm1_is_only_moderate():
    """The caveat is attached to the Strong argument, not to PM1 generally."""
    out = _gate_criteria_applicability(
        _criteria("PM1", "PM5"), "FBN1", _vep(_p("Cys", 302, "Arg")))
    assert _by_code(out, "PM1")["criteria_strength"] == "PM1_Moderate"
    assert _by_code(out, "PM5")["status"] == "met"


def test_gate_suppresses_pm1_outside_the_enumeration():
    out = _gate_criteria_applicability(
        _criteria("PM1"), "FBN1", _vep(_p("Arg", 5, "Trp")))
    assert _by_code(out, "PM1")["status"] == "not_met"


def test_kcnq1_pore_helix_range_now_parses():
    """GN112 writes its hotspot as "amino acids 300 to 320", which the dash-only
    pattern could not see — so KCNQ1's PM1 had no range and the residue sub-gate
    silently did nothing for it."""
    ranges, _ = _parse_pm1_hotspots(_crit_text("KCNQ1", "PM1"))
    assert (300, 320) in ranges


def test_kcnq1_pm1_requires_pm2():
    """GN112: "The variant must be rare (meeting PM2_Supporting) in order to be
    considered for PM1.\""""
    vep = {"vep": {"ok": True, "most_severe_consequence": "missense_variant",
                   "hgvsp": "NP_000209.2:p.Gly314Ser"}}
    without = _gate_criteria_applicability(_criteria("PM1"), "KCNQ1", vep)
    assert _by_code(without, "PM1")["status"] == "not_met"
    assert "rare" in _by_code(without, "PM1")["evidence"]

    with_pm2 = _gate_criteria_applicability(
        _criteria("PM2", "PM1"), "KCNQ1", vep)
    assert _by_code(with_pm2, "PM1")["status"] == "met"


def test_kcnq1_pm1_suppressed_outside_the_pore_helix():
    vep = {"vep": {"ok": True, "most_severe_consequence": "missense_variant",
                   "hgvsp": "NP_000209.2:p.Arg591His"}}
    out = _gate_criteria_applicability(_criteria("PM2", "PM1"), "KCNQ1", vep)
    assert _by_code(out, "PM1")["status"] == "not_met"
    assert "300, 320" in _by_code(out, "PM1")["evidence"].replace("(", "").replace(")", "")


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
