"""Unit tests for the deterministic (server-owned) PM5 criterion.

PM5 was moved from AI-evaluated to server-owned (2026-07-12): benchmarking
showed the AI applied the "PM5 MUST-APPLY" rule on only ~4-5 of ~22 qualifying
variants even though get_pm5_evidence had already retrieved the ≥2★
same-residue candidates. `_clinvar_pm5_criterion` now derives PM5 directly from
that retrieval. These tests pin every branch: the ≥2★ gate, the RASopathy
≥2-distinct-change guard, and the not-applicable / no-candidate paths.
"""
from __future__ import annotations

from backend.acmg.hard_coded import _clinvar_pm5_criterion, compute_hard_coded_criteria


def _cand(name, cs="Pathogenic", stars=2):
    return {"name": name, "clinical_significance": cs, "stars": stars, "tier": "P"}


def _ev(candidates, *, ok=True, not_applicable=False, pos=403, gene="MYH7"):
    d = {"ok": ok, "gene": gene, "protein_position": pos,
         "candidates": candidates, "count": len(candidates),
         "count_two_star": sum(1 for c in candidates if (c.get("stars") or 0) >= 2)}
    if not_applicable:
        d["not_applicable"] = True
    return d


def test_pm5_met_on_two_star_candidate():
    ev = _ev([_cand("NM_000257.4(MYH7):c.1207C>T (p.Arg403Trp)", stars=3)])
    crit = _clinvar_pm5_criterion(ev, "MYH7")
    assert crit["status"] == "met"
    assert crit["criteria_strength"] == "PM5_Moderate"
    assert crit["source"] == "hard_coded"
    assert "residue 403" in crit["evidence"]


def test_pm5_not_met_when_only_one_star():
    ev = _ev([_cand("NM_000257.4(MYH7):c.1207C>T (p.Arg403Trp)", stars=1)])
    crit = _clinvar_pm5_criterion(ev, "MYH7")
    assert crit["status"] == "not_met"
    assert crit["criteria_strength"] is None
    assert "below the PM5 confidence bar" in crit["evidence"]


def test_pm5_not_met_when_no_candidate():
    crit = _clinvar_pm5_criterion(_ev([]), "MYH7")
    assert crit["status"] == "not_met"
    assert "No other P/LP missense change" in crit["evidence"]


def test_pm5_not_applicable_passthrough():
    for ev in ({"ok": False}, _ev([], not_applicable=True), None):
        crit = _clinvar_pm5_criterion(ev, "MYH7")
        assert crit["status"] == "not_met"
        assert crit["criteria_strength"] is None


def test_pm5_rasopathy_moderate_needs_only_one_change():
    """CORRECTED 2026-09-01. This test previously asserted that a single ≥2★
    change must NOT fire PM5 on a RASopathy gene. That was a misreading: the
    ≥2-distinct-changes requirement belongs to the RASopathy STRONG row ("≥2
    different [likely] pathogenic residues changes at the same codon observed in
    ≥5 probands"), while the MODERATE row of all 15 RASopathy genes reads simply
    "1 [likely] pathogenic residue change at the same codon".

    Gating Moderate on ≥2 suppressed legitimate calls — an UNDER-call and a
    divergence from the published spec, which is the class of defect that
    matters even when it moves us closer to the experts."""
    one = _ev([_cand("NM_002834.5(PTPN11):c.922A>G (p.Asn308Asp)", stars=3)],
              gene="PTPN11", pos=308)
    crit = _clinvar_pm5_criterion(one, "PTPN11")
    assert crit["status"] == "met", crit["evidence"]
    assert crit["criteria_strength"] == "PM5_Moderate"
    assert "distinct alt-AA change" in crit["evidence"]
    assert "proband count HeartVar does not hold" in crit["evidence"]

    two = _ev([_cand("NM_002834.5(PTPN11):c.922A>G (p.Asn308Asp)", stars=3),
               _cand("NM_002834.5(PTPN11):c.923A>C (p.Asn308Thr)", stars=2)],
              gene="PTPN11", pos=308)
    got = _clinvar_pm5_criterion(two, "PTPN11")
    assert got["status"] == "met" and got["criteria_strength"] == "PM5_Moderate"


def test_pm5_strength_follows_the_comparators_classification_on_cm_genes():
    """The 8 Cardiomyopathy genes tie PM5's strength to the comparator:
    Moderate when the same-codon variant is classified PATHOGENIC, Supporting
    when it is LIKELY pathogenic. The builder hardcoded Moderate for every gene,
    so an LP-only comparator scored +2 where the spec allows +1."""
    lp_only = _ev([_cand("NM_000257.4(MYH7):c.1207C>T (p.Arg403Trp)", stars=2,
                         cs="Likely pathogenic")],
                  gene="MYH7", pos=403)
    got = _clinvar_pm5_criterion(lp_only, "MYH7")
    assert got["status"] == "met"
    assert got["criteria_strength"] == "PM5_Supporting", got["evidence"]
    assert "Likely pathogenic" in got["evidence"]

    full_p = _ev([_cand("NM_000257.4(MYH7):c.1207C>T (p.Arg403Trp)", stars=2,
                        cs="Pathogenic")],
                 gene="MYH7", pos=403)
    assert _clinvar_pm5_criterion(full_p, "MYH7")["criteria_strength"] == "PM5_Moderate"

    ras_lp = _ev([_cand("NM_002834.5(PTPN11):c.922A>G (p.Asn308Asp)", stars=3,
                        cs="Likely pathogenic")],
                 gene="PTPN11", pos=308)
    assert _clinvar_pm5_criterion(ras_lp, "PTPN11")["criteria_strength"] == "PM5_Moderate"


def test_compute_hard_coded_emits_pm5():
    ev = {"clinvar_pm5_candidates": _ev(
        [_cand("NM_000257.4(MYH7):c.1207C>T (p.Arg403Trp)", stars=3)])}
    out = compute_hard_coded_criteria(ev, {"inheritance_input": "AD"}, "MYH7")
    pm5 = [c for c in out if c["code"] == "PM5"]
    assert len(pm5) == 1
    assert pm5[0]["status"] == "met"
    assert pm5[0]["source"] == "hard_coded"


if __name__ == "__main__":
    import sys
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn(); print(f"PASS {name}")
            except AssertionError as exc:
                failures += 1; print(f"FAIL {name}: {exc}")
    sys.exit(1 if failures else 0)
