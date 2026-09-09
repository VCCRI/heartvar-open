"""Unit tests for the ClinGen-SVI calibrated REVEL PP3/BP4 logic and the BP7
conservation (phyloP) gate added in the conservation + calibrated-REVEL change.

Three layers, all fully offline (no network / DB / clock):

  * ``_calibrated_revel_call`` — the pure REVEL→strength mapping, at every
    band boundary (Pejaver 2022 thresholds).
  * ``_eval_pp3_bp4`` — REVEL-primary behaviour for missense, plus the
    consensus fallback when REVEL is indeterminate or absent, and the
    untouched non-coding/splice arms.
  * ``_eval_bp7`` — the "nucleotide not highly conserved" (phyloP100way > 2.0 =
    highly conserved, MM-VCEP) gate,
    including the missing-phyloP passthrough.

Runnable directly (``python -m backend.tests.test_calibrated_insilico``).
"""
from __future__ import annotations

from backend.acmg.hard_coded import (
    _calibrated_revel_call,
    _clamp_insilico_strength,
    _eval_pp3_bp4,
    _eval_bp7,
    _insilico_strength_ceiling,
    apply_cross_criterion_exclusions,
)


def _crit(code, strength):
    return {"code": code, "status": "met", "criteria_strength": strength,
            "direction": "benign" if code.startswith("B") else "pathogenic",
            "evidence": ""}


def _by_code(entries: list[dict]) -> dict[str, dict]:
    return {e["code"]: e for e in entries}


def _vep(**kw) -> dict:
    v = {"ok": True, "most_severe_consequence": "missense_variant"}
    v.update(kw)
    return v


def _ev(am: float | None = None, spliceai_ok: bool = False) -> dict:
    e: dict = {}
    if am is not None:
        e["alphamissense"] = {"ok": True, "am_pathogenicity": am}
    if spliceai_ok:
        e["spliceai"] = {"ok": True}
    return e


def test_revel_pathogenic_bands():
    assert _calibrated_revel_call(0.999) == ("PP3", "PP3_Strong")
    assert _calibrated_revel_call(0.932) == ("PP3", "PP3_Strong")
    assert _calibrated_revel_call(0.9319) == ("PP3", "PP3_Moderate")
    assert _calibrated_revel_call(0.773) == ("PP3", "PP3_Moderate")
    assert _calibrated_revel_call(0.7729) == ("PP3", "PP3_Supporting")
    assert _calibrated_revel_call(0.644) == ("PP3", "PP3_Supporting")


def test_revel_benign_bands():
    assert _calibrated_revel_call(0.290) == ("BP4", "BP4_Supporting")
    assert _calibrated_revel_call(0.183) == ("BP4", "BP4_Moderate")
    assert _calibrated_revel_call(0.05) == ("BP4", "BP4_Moderate")
    assert _calibrated_revel_call(0.016) == ("BP4", "BP4_Strong")
    assert _calibrated_revel_call(0.010) == ("BP4", "BP4_Strong")
    assert _calibrated_revel_call(0.003) == ("BP4", "BP4_VeryStrong")
    assert _calibrated_revel_call(0.001) == ("BP4", "BP4_VeryStrong")


def test_revel_indeterminate_and_missing_return_none():
    assert _calibrated_revel_call(0.30) is None
    assert _calibrated_revel_call(0.50) is None
    assert _calibrated_revel_call(0.6439) is None
    assert _calibrated_revel_call(None) is None
    assert _calibrated_revel_call(".") is None


def test_pp3_revel_primary_strength_tiers():
    for revel, strength in [(0.95, "PP3_Strong"), (0.80, "PP3_Moderate"),
                            (0.65, "PP3_Supporting")]:
        out, _, _ = _eval_pp3_bp4(_ev(), _vep(revel_score=revel),
                                  "missense_variant", 0.0)
        d = _by_code(out)
        assert d["PP3"]["status"] == "met"
        assert d["PP3"]["criteria_strength"] == strength
        assert d["BP4"]["status"] == "not_met"


def test_bp4_revel_primary_strength_tiers():
    for revel, strength in [(0.25, "BP4_Supporting"), (0.10, "BP4_Moderate"),
                            (0.012, "BP4_Strong"), (0.002, "BP4_VeryStrong")]:
        out, _, _ = _eval_pp3_bp4(_ev(), _vep(revel_score=revel),
                                  "missense_variant", 0.0)
        d = _by_code(out)
        assert d["BP4"]["status"] == "met"
        assert d["BP4"]["criteria_strength"] == strength
        assert d["PP3"]["status"] == "not_met"


def test_revel_decisive_overrides_disagreeing_other_tools():
    out, _, _ = _eval_pp3_bp4(_ev(am=0.99), _vep(revel_score=0.04, cadd_phred=35),
                              "missense_variant", 0.0)
    d = _by_code(out)
    assert d["BP4"]["status"] == "met"
    assert d["BP4"]["criteria_strength"] == "BP4_Moderate"
    assert d["PP3"]["status"] == "not_met"


def test_revel_shortcircuit_gated_to_missense():
    out, _, _ = _eval_pp3_bp4(_ev(), _vep(revel_score=0.95),
                              "splice_donor_variant", 0.0)
    d = _by_code(out)
    assert not (d["PP3"]["status"] == "met"
                and d["PP3"]["criteria_strength"] == "PP3_Strong")


def test_benign_revel_does_not_bury_splice_signal():
    out, _, _ = _eval_pp3_bp4(_ev(spliceai_ok=True), _vep(revel_score=0.01),
                              "missense_variant", 0.95)
    d = _by_code(out)
    assert d["BP4"]["status"] == "not_met"


def test_indeterminate_revel_falls_back_to_consensus_benign():
    out, _, _ = _eval_pp3_bp4(_ev(am=0.05), _vep(revel_score=0.50, cadd_phred=2),
                              "missense_variant", 0.0)
    d = _by_code(out)
    assert d["BP4"]["status"] == "met"
    assert d["BP4"]["criteria_strength"] == "BP4_Supporting"


def test_indeterminate_revel_falls_back_to_consensus_pathogenic():
    out, _, _ = _eval_pp3_bp4(_ev(am=0.99), _vep(revel_score=0.50, cadd_phred=30),
                              "missense_variant", 0.0)
    d = _by_code(out)
    assert d["PP3"]["status"] == "met"
    assert d["PP3"]["criteria_strength"] == "PP3_Supporting"


def test_absent_revel_falls_back_to_consensus():
    out, _, _ = _eval_pp3_bp4(_ev(am=0.99), _vep(cadd_phred=30),
                              "missense_variant", 0.0)
    d = _by_code(out)
    assert d["PP3"]["status"] == "met"
    assert d["PP3"]["criteria_strength"] == "PP3_Supporting"


def test_absent_revel_no_tools_both_not_met():
    out, _, _ = _eval_pp3_bp4(_ev(), _vep(), "missense_variant", 0.0)
    d = _by_code(out)
    assert d["PP3"]["status"] == "not_met"
    assert d["BP4"]["status"] == "not_met"


def test_noncoding_splice_arm_unaffected_by_revel_logic():
    out, ncs, clear = _eval_pp3_bp4(_ev(spliceai_ok=True), _vep(
        most_severe_consequence="intron_variant"), "intron_variant", 0.02)
    d = _by_code(out)
    assert d["BP4"]["status"] == "met"
    assert d["BP4"]["criteria_strength"] == "BP4_Supporting"
    assert ncs is True and clear is True


def test_bp7_fires_when_not_conserved():
    e = _eval_bp7("synonymous_variant", 0.02, sa_ok=True,
                  noncoding_splice=False, splice_clear=True, phylop=0.05)
    assert e["status"] == "met"
    assert e["criteria_strength"] == "BP7_Supporting"


def test_bp7_fires_at_mmvcep_boundary():
    e = _eval_bp7("synonymous_variant", 0.02, sa_ok=True,
                  noncoding_splice=False, splice_clear=True, phylop=2.0)
    assert e["status"] == "met"


def test_bp7_blocked_when_highly_conserved():
    e = _eval_bp7("synonymous_variant", 0.02, sa_ok=True,
                  noncoding_splice=False, splice_clear=True, phylop=3.0)
    assert e["status"] == "not_met"
    assert "conserved" in e["evidence"].lower()


def test_bp7_missing_phylop_withheld_fail_closed():
    e = _eval_bp7("synonymous_variant", 0.02, sa_ok=True,
                  noncoding_splice=False, splice_clear=True, phylop=None)
    assert e["status"] == "not_met"
    assert "could not be determined" in e["evidence"].lower()


def test_bp7_splice_impact_blocks_regardless_of_conservation():
    e = _eval_bp7("synonymous_variant", 0.8, sa_ok=True,
                  noncoding_splice=False, splice_clear=False, phylop=0.05)
    assert e["status"] == "not_met"


def test_clamp_insilico_strength():
    assert _clamp_insilico_strength("PP3_Strong", "Supporting") == "PP3_Supporting"
    assert _clamp_insilico_strength("BP4_VeryStrong", "Supporting") == "BP4_Supporting"
    assert _clamp_insilico_strength("PP3_Moderate", "Supporting") == "PP3_Supporting"
    assert _clamp_insilico_strength("PP3_Supporting", "Supporting") == "PP3_Supporting"
    assert _clamp_insilico_strength("PP3_Strong", None) == "PP3_Strong"


def test_strength_ceiling_supporting_engine_wide():
    for g in ("MAP2K2", "MYBPC3", "FBN1", "KCNQ1", "KCNH2", "MYH7",
              "ATM", "APOB", "ZZZZ9", ""):
        assert _insilico_strength_ceiling(g, {}) == "Supporting", g


def test_pp3_capped_to_supporting_for_all_genes():
    vep = _vep(revel_score=0.966)
    for g in ("MYBPC3", "MAP2K2", "ATM", "ZZZZ9"):
        out, _, _ = _eval_pp3_bp4(_ev(), vep, "missense_variant", 0.0,
                                  _insilico_strength_ceiling(g, {}))
        assert _by_code(out)["PP3"]["criteria_strength"] == "PP3_Supporting", g


def test_bp4_capped_for_gof_gene():
    vep = _vep(revel_score=0.001)
    out, _, _ = _eval_pp3_bp4(_ev(), vep, "missense_variant", 0.0,
                              _insilico_strength_ceiling("MRAS", {}))
    assert _by_code(out)["BP4"]["criteria_strength"] == "BP4_Supporting"


def test_pp3_pm1_cap_demotes_pm1_when_combined_exceeds_strong():
    cleaned, _ = apply_cross_criterion_exclusions(
        [_crit("PP3", "PP3_Strong"), _crit("PM1", "PM1_Moderate")], "ATM", "AD")
    d = _by_code(cleaned)
    assert d["PP3"]["criteria_strength"] == "PP3_Strong"
    assert d["PM1"]["status"] == "not_met"


def test_pp3_pm1_no_cap_when_within_strong():
    cleaned, _ = apply_cross_criterion_exclusions(
        [_crit("PP3", "PP3_Moderate"), _crit("PM1", "PM1_Moderate")], "ATM", "AD")
    d = _by_code(cleaned)
    assert d["PP3"]["criteria_strength"] == "PP3_Moderate"
    assert d["PM1"]["status"] == "met"
    assert d["PM1"]["criteria_strength"] == "PM1_Moderate"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"ok  {fn.__name__}")
    print(f"\n{len(fns)} tests passed")


def test_alphamissense_benign_cutoff_matches_the_published_am_class_boundary():
    """The AlphaMissense benign cutoff is 0.34, verified against the published
    table rather than from memory.

    In data/AlphaMissense_hg38.tsv.gz the authors' own `am_class` column has
    `likely_benign` topping out at 0.34 and `ambiguous` beginning at 0.34.
    hard_coded.py had 0.337 — a typo that had drifted away from the three places
    already saying 0.340 (clients/alphamissense.py:7, prompt.py:1286,
    static/heartvar.js:6645).

    While it stood, a score in [0.337, 0.340) — 0.193% of missense rows — was
    shown to the curator as likely-benign, described to the MODEL as
    BP4-supporting, and then excluded from tools_benign by the engine. Three
    components disagreeing about one variant, which is why this is pinned."""
    from backend.acmg import hard_coded as H
    import inspect, re
    src = inspect.getsource(H._eval_pp3_bp4)
    m = re.search(r'"alphamissense":\s*([0-9.]+),?\s*\n\s*\}', src)
    vals = re.findall(r'"alphamissense":\s*([0-9.]+)', src)
    assert vals, "could not find the alphamissense thresholds"
    assert "0.564" in vals, f"pathogenic side changed: {vals}"
    assert "0.340" in vals or "0.34" in vals, (
        f"benign side must be 0.34 per the published am_class boundary, got {vals}"
    )
    assert "0.337" not in vals, "the 0.337 typo is back"

    out, _, _ = _eval_pp3_bp4(_ev(am=0.338), _vep(revel_score=0.50, cadd_phred=2),
                              "missense_variant", 0.0)
    assert _by_code(out)["BP4"]["status"] == "met"
