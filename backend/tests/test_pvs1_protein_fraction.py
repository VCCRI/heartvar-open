"""PVS1's ">10% of the protein" rule must measure the PROTEIN.

THE BUG. `exon_lengths` comes from `_fetch_exon_lengths`, which reads Ensembl
`/lookup/id?expand=1` and returns `e["end"] - e["start"] + 1` per exon — GENOMIC
exon spans. Their sum is therefore the mRNA length **including both UTRs**, not
the coding sequence. Two places divided by it:

  * `_pvs1_splice_strength` used `exon_len / sum(exon_lengths)`. Denominator too
    big -> fraction too small -> a skip removing 15% of the protein scored ~9%
    and was graded PVS1_Moderate (+2) instead of PVS1_Strong (+4).
  * `_pvs1_nmd_escape_strength` PREFERRED `sum(el) / 3.0` over `uniprot.length`,
    under a comment that called it "CDS nt -> residues". There the inflation
    runs the other way: `(total - pos) / total` with `total` inflated and `pos`
    fixed OVERSTATES the fraction removed, promoting Moderate (+2) to Strong
    (+4). That is the pathogenic direction.

Inflation of sum/3 over true protein length, MANE mRNA length vs UniProt residue
count: KCNQ1 ~1.59x, LMNA ~1.60x, PKP2 ~1.66x.

`exon_len % 3` is a valid frame test only for an INTERNAL exon — exon 1 carries
the 5'UTR and the last exon the 3'UTR, so their genomic length says nothing
about the coding length skipped. Frame is now treated as UNKNOWN there, which
keeps the variant out of the out-of-frame branch (the one returning full +8).
"""
from __future__ import annotations

from backend.acmg.hard_coded import (
    _pvs1_nmd_escape_strength,
    _pvs1_splice_strength,
)

_EXONS = [600, 150, 150, 150, 150, 150, 150, 150, 150, 600]
_UNIPROT_400AA = {"uniprot": {"length": 400}}


def _donor(intron_n, exons=_EXONS, ev=None):
    """splice_donor at intron N hits exon N."""
    return _pvs1_splice_strength(
        {"most_severe_consequence": "splice_donor_variant",
         "intron": f"{intron_n}/9", "exon_lengths": exons}, ev)


def test_internal_exon_skip_graded_against_the_protein_not_the_transcript():
    """Exon 5 is 150 nt. Against the 1200 nt CDS that is 12.5% -> Strong.
    Against the 2400 nt mRNA it was 6.3% -> Moderate. This is the test that
    would have caught the bug."""
    tier, why = _donor(5, ev=_UNIPROT_400AA)
    assert tier == "PVS1_Strong", f"got {tier}: {why}"
    assert "of the protein" in why
    assert "13%" in why or "12%" in why


def test_without_the_protein_length_the_denominator_is_named_as_inflated():
    """Fallback path. Still grades (returning None would leave the full +8),
    but must not silently present an mRNA-based percentage as a protein one."""
    tier, why = _donor(5, ev=None)
    assert tier == "PVS1_Moderate"
    assert "incl. UTRs" in why
    assert "understates" in why


def test_a_genuinely_small_internal_skip_stays_moderate():
    small = [600, 150, 150, 30, 150, 150, 150, 150, 150, 600]
    tier, why = _pvs1_splice_strength(
        {"most_severe_consequence": "splice_donor_variant",
         "intron": "4/9", "exon_lengths": small}, _UNIPROT_400AA)
    assert tier == "PVS1_Moderate", why


def test_terminal_exon_frame_is_not_asserted():
    """Exon 1 is 600 nt (600 % 3 == 0 here by construction), but that number is
    mostly UTR, so neither "in frame" nor "out of frame" is knowable. The
    out-of-frame branch returns full PVS1 (+8) and must not be reachable."""
    tier, why = _donor(1, ev=_UNIPROT_400AA)
    assert tier != "PVS1", f"terminal exon reached the full +8 branch: {why}"
    assert "not determinable" in why


def test_terminal_exon_out_of_frame_length_never_reaches_full_pvs1():
    odd = [601, 150, 150, 150, 150, 150, 150, 150, 150, 599]
    tier, why = _pvs1_splice_strength(
        {"most_severe_consequence": "splice_donor_variant",
         "intron": "1/9", "exon_lengths": odd}, _UNIPROT_400AA)
    assert tier != "PVS1", (
        "601 % 3 != 0 on a UTR-containing exon is not evidence of a frameshift"
    )


def test_terminal_exon_percentage_is_treated_as_an_upper_bound():
    """Exon 1's 600 nt is 50% of the 1200 nt CDS on paper, but part of it is
    5'UTR, so >10% of the PROTEIN is not established. Conservative Moderate."""
    tier, why = _donor(1, ev=_UNIPROT_400AA)
    assert tier == "PVS1_Moderate", why
    assert "not established" in why


def test_internal_out_of_frame_skip_still_reaches_full_pvs1():
    """Guard against over-suppression: the +8 branch must still work where the
    frame test is actually valid."""
    exons = [600, 150, 151, 150, 150, 150, 150, 150, 150, 600]
    tier, why = _pvs1_splice_strength(
        {"most_severe_consequence": "splice_donor_variant",
         "intron": "3/9", "exon_lengths": exons}, _UNIPROT_400AA)
    assert tier == "PVS1", why


def test_nmd_escape_prefers_uniprot_length_over_the_exon_sum():
    """KCNQ1-shaped case. Real protein 676 aa; a PTC at residue 650 removes
    3.8% -> Moderate. With the mRNA-derived 1020 "residues" it scored
    (1020-650)/1020 = 36% -> Strong. A 2-point over-call in the pathogenic
    direction, and the exon sum was PREFERRED over the true length."""
    tier, why = _pvs1_nmd_escape_strength(
        {"protein_start": 650, "exon_lengths": [3060]},
        {"uniprot": {"length": 676}},
    )
    assert tier == "PVS1_Moderate", f"got {tier}: {why}"
    assert "/676" in why


def test_nmd_escape_still_falls_back_to_the_exon_sum():
    tier, _ = _pvs1_nmd_escape_strength(
        {"protein_start": 1271, "exon_lengths": [3822]}, {})
    assert tier == "PVS1_Moderate"


def test_nmd_escape_large_truncation_still_strong():
    tier, _ = _pvs1_nmd_escape_strength(
        {"protein_start": 300}, {"uniprot": {"length": 676}})
    assert tier == "PVS1_Strong"
