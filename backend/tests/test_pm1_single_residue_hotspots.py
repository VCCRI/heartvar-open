"""_parse_pm1_hotspots — enumerated SINGLE residues count, not just N-N spans.

WHY THIS FILE EXISTS. PTPN11 GN043's PM1_Moderate enumerates a mix of single
residues and ranges:

    "Applicable only to critical and well-established functional domains
     available in the supplementary table (Directly interacting residues between
     N-SH2 and PTPN domains [AA 4, AA 7-9, AA 58-63, AA 69-77, AA 247, AA 251,
     AA 255, AA 256, AA 258, AA 261, AA 265, AA 278-281, AA 284]). Not
     applicable to specific amino acid residues (see PM5)."

The parser read only N-N spans, so it returned four ranges and dropped all nine
singles. `_gate_criteria_applicability`'s PM1 suppress sub-gate then removed a
PM1 the spec licenses, which is why PTPN11 NM_002834.4:c.781C>T (p.Leu261Phe),
expert Pathogenic, came back VUS without PM1 — codon 261 is in the spec's own
list.

The trailing exclusion clause is why this sat parked, and it is settled by
measurement rather than by reading: thirteen genes carry that exact clause and
TWELVE of them enumerate no single residues at all, so on those twelve it cannot
mean "exclude the enumerated singles". See test_the_exclusion_clause_reading.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from backend.acmg.hard_coded import _parse_pm1_hotspots

SPEC = json.loads(
    (Path(__file__).resolve().parent.parent / "data" / "vcep_criteria_spec.json").read_text()
)["genes"]

PTPN11_SINGLES = [4, 247, 251, 255, 256, 258, 261, 265, 284]
PTPN11_RANGES = [(7, 9), (58, 63), (69, 77), (278, 281)]


def _pm1_text(gene: str) -> str:
    pm1 = SPEC[gene].get("PM1") or {}
    return json.dumps(pm1.get("strengths") or {}) + " " + str(pm1.get("acmg_summary") or "")


def _applicable_pm1_genes() -> list[str]:
    return [g for g in sorted(SPEC)
            if (SPEC[g].get("PM1") or {}).get("applicability") != "not_applicable"]


def test_ptpn11_keeps_its_ranges_and_gains_its_nine_singles():
    ranges, _ = _parse_pm1_hotspots(_pm1_text("PTPN11"))
    assert all(r in ranges for r in PTPN11_RANGES), ranges
    assert sorted(lo for lo, hi in ranges if lo == hi) == PTPN11_SINGLES, ranges


@pytest.mark.parametrize("residue", PTPN11_SINGLES)
def test_each_enumerated_ptpn11_residue_is_in_scope(residue: int):
    """Codon 261 is the one that produced the VUS; the other eight are in the
    same list and must behave identically."""
    ranges, _ = _parse_pm1_hotspots(_pm1_text("PTPN11"))
    assert any(lo <= residue <= hi for lo, hi in ranges), residue


@pytest.mark.parametrize("residue", [1, 100, 200, 300, 500])
def test_residues_outside_the_enumeration_stay_out(residue: int):
    """The fix must widen PM1 to the spec's list, not to everything. PM1 is +2
    pathogenic, so an over-wide parse is the dangerous direction."""
    ranges, _ = _parse_pm1_hotspots(_pm1_text("PTPN11"))
    assert not any(lo <= residue <= hi for lo, hi in ranges), residue


def test_ptpn11_is_the_only_gene_that_gains_a_single_residue():
    """THE BLAST RADIUS, pinned. The parked concern was that "all 13 RASopathy
    genes carry the same clause, so a parser change hits every one of them".
    Measured: it does not. Only PTPN11 enumerates single residues at all."""
    gainers = {}
    for gene in _applicable_pm1_genes():
        ranges, _ = _parse_pm1_hotspots(_pm1_text(gene))
        singles = sorted(lo for lo, hi in ranges if lo == hi)
        if singles:
            gainers[gene] = singles
    assert gainers == {"PTPN11": PTPN11_SINGLES}, gainers


def test_the_exclusion_clause_reading():
    """The evidence that settles the interpretation, kept executable.

    Twelve of the thirteen genes carrying "not applicable to specific amino acid
    residues (see PM5)" enumerate NO single residues, so on them the clause
    cannot mean "exclude the enumerated singles". It is a gene-independent
    instruction that residue-specific arguments belong to PM5.
    """
    clause = "not applicable to specific amino acid"
    carriers = [g for g in _applicable_pm1_genes() if clause in _pm1_text(g).lower()]
    assert len(carriers) >= 11, carriers
    with_singles = [g for g in carriers
                    if any(lo == hi for lo, hi in _parse_pm1_hotspots(_pm1_text(g))[0])]
    assert with_singles == ["PTPN11"], with_singles


def test_a_word_written_span_is_not_also_read_as_a_single():
    """KCNQ1 writes "amino acids 300 to 320". Marking only DASHED spans as
    consumed let "amino acids 300" match as a single and produced a spurious
    (300, 300) — subsumed by (300, 320) so it changed nothing, but a mis-parse.
    """
    ranges, _ = _parse_pm1_hotspots(_pm1_text("KCNQ1"))
    assert (300, 320) in ranges, ranges
    assert not [r for r in ranges if r[0] == r[1]], ranges


def test_stray_numbers_in_rule_prose_are_not_read_as_hotspots():
    """The single pattern is anchored on an amino-acid/codon/residue word for
    the same reason the span patterns are: rule text carries proband counts and
    odds-ratio bounds that must not become hotspots."""
    ranges, _ = _parse_pm1_hotspots(
        "Observed in 12 unrelated probands; odds ratio 8 (95% CI 3-19). "
        "Applicable to AA 55 and codons 90-95."
    )
    singles = sorted(lo for lo, hi in ranges if lo == hi)
    assert singles == [55], ranges
    assert (90, 95) in ranges, ranges
    assert not any(lo == 12 or lo == 8 for lo, hi in ranges), ranges
