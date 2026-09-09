"""Cross-transcript residue numbering in the ClinVar evidence surfaces.

Regression tests for the evidence-retrieval gap reported on GATA4
``NM_002052.5:c.907G>T`` (p.Gly303Trp). These began as xfail(strict=True)
tests documenting the bug; the fix (MANE-anchored matching) landed and the
markers were removed, so they are now live regression guards.

The facts these tests are built on (all verified against live ClinVar
VCV003897161 / VCV001519604 and Ensembl VEP, 2026-08-24):

  * GATA4's MANE Select transcript is ``NM_001308093.3`` / ``NP_001295022.1``
    (443 aa). ``NM_002052.5`` / ``NP_002043.2`` is 442 aa — the MANE isoform
    carries one extra residue (Val) inserted at position 206, so every residue
    downstream of 205 is numbered **+1** relative to NM_002052.5.
  * ClinVar's ``Name`` field — the only place ``clinvar.py`` reads a residue
    number from — is written on ClinVar's preferred transcript, which for
    GATA4 is the MANE Select ``NM_001308093.3``. 1074 of 1086 GATA4 rows in
    the local mirror are named on it.
  * Therefore ``NM_002052.5:p.Gly303`` and ``NM_001308093.3:p.Gly304`` are the
    SAME residue (chr8:11750234-11750236, GRCh38), and
    ``NM_002052.5:p.Gly296Ser`` (c.886G>A) is the SAME variant as
    ``NM_001308093.3:p.Gly297Ser`` (c.889G>A), chr8:11750213 G>A.

``_pm5_evidence_sync`` compares ``protein_position`` (derived from VEP's
hgvsp on the transcript the CURATOR supplied) against the integer parsed out
of ClinVar's Name (always ClinVar's preferred transcript). No transcript
reconciliation happens anywhere on that path, so for any gene whose isoforms
differ in numbering the comparison is between two different residues.

THE FIX. ``_mane_hgvsp_from_vep`` / ``_mane_hgvsc_from_vep`` read the MANE
Select row out of VEP's ``transcript_consequences_all``; the ClinVar client
MATCHES on those and DISPLAYS the curator's own numbering. Where no MANE
annotation exists, everything falls back to the curator's transcript — i.e.
unchanged pre-fix behaviour — rather than guessing an offset.
"""
from __future__ import annotations

import contextlib
import os
import sqlite3
import tempfile
from pathlib import Path

import backend.clients.clinvar as clinvar

_DDL = """
CREATE TABLE variants (
    variation_id            INTEGER,
    obj_type                TEXT,
    name                    TEXT,
    gene_symbol             TEXT,
    clinical_significance   TEXT,
    review_status           TEXT,
    number_submitters       INTEGER,
    phenotype_list          TEXT,
    chromosome_accession    TEXT,
    start                   INTEGER,
    stop                    INTEGER,
    reference_allele        TEXT,
    alternate_allele        TEXT,
    assembly                TEXT,
    last_evaluated          TEXT
)
"""
_INSERT_SQL = "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"

_MULTI = "criteria provided, multiple submitters, no conflicts"
_SINGLE = "criteria provided, single submitter"
_NOCRIT = "no assertion criteria provided"

_GATA4_ROWS = [
    (2820525, "NM_001308093.3(GATA4):c.887G>C (p.Cys296Ser)",
     "Pathogenic", _SINGLE, 1, 11750211),
    (9030, "NM_001308093.3(GATA4):c.889G>A (p.Gly297Ser)",
     "Pathogenic/Likely pathogenic", _MULTI, 4, 11750213),
    (30098, "NM_001308093.3(GATA4):c.889G>T (p.Gly297Cys)",
     "Pathogenic", _NOCRIT, 1, 11750213),
    (30106, "NM_001308093.3(GATA4):c.889G>C (p.Gly297Arg)",
     "Pathogenic", _NOCRIT, 1, 11750213),
    (700648, "NM_001308093.3(GATA4):c.909C>T (p.His303=)",
     "Likely benign", _MULTI, 4, 11750233),
    (3897161, "NM_001308093.3(GATA4):c.910G>A (p.Gly304Arg)",
     "Uncertain significance", _MULTI, 2, 11750234),
    (1519604, "NM_001308093.3(GATA4):c.912G>A (p.Gly304=)",
     "Uncertain significance", _SINGLE, 1, 11750236),
]


def _rows():
    out = []
    for vid, name, sig, rev, nsub, gpos in _GATA4_ROWS:
        out.append((
            vid, "single nucleotide variant", name, "GATA4", sig, rev, nsub,
            "Atrioventricular septal defect 4", "NC_000008.11", gpos, gpos,
            "na", "na", "GRCh38", "2026-01-01",
        ))
    return out


@contextlib.contextmanager
def _temp_db():
    fd, path = tempfile.mkstemp(suffix=".db", prefix="clinvar_txnum_")
    os.close(fd)
    conn = sqlite3.connect(path)
    try:
        conn.execute(_DDL)
        conn.executemany(_INSERT_SQL, _rows())
        conn.commit()
    finally:
        conn.close()
    original = clinvar.DB_PATH
    clinvar.DB_PATH = Path(path)
    try:
        yield
    finally:
        clinvar.DB_PATH = original
        with contextlib.suppress(OSError):
            os.unlink(path)


def test_exact_variant_lookup_finds_record_named_on_another_transcript():
    """``NM_002052.5:c.886G>A`` and ``NM_001308093.3:c.889G>A`` are the SAME
    variant (chr8:11750213 G>A). ClinVar holds it as Pathogenic/Likely
    pathogenic at 2★. Matching only the curator's bare ``c.`` token against
    ClinVar's Name reported it absent, silently removing PP5. Passing the
    MANE coding change finds it, and says which layer matched."""
    with _temp_db():
        res = clinvar._query_sync("GATA4", "c.886G>A", "c.889G>A")
        assert res["found"] is True, (
            "GATA4 NM_002052.5:c.886G>A (= NM_001308093.3:c.889G>A, "
            "VCV000009030, Pathogenic/Likely pathogenic, 2★) reported as "
            "absent from ClinVar"
        )
        assert res["matched_on_mane"] is True
        assert res["records"][0]["variation_id"] == 9030
        assert res["records"][0]["matched_via"] == "hgvs_c_mane"


def test_curator_token_alone_still_misses_it_documented_limitation():
    """The pre-fix path, pinned deliberately. With no MANE token supplied
    there is nothing in the mirror to match, so the answer is still "not
    found". This is why the lookup is chained on VEP — it is not a fallback
    that quietly works, it is a fallback that quietly cannot."""
    with _temp_db():
        res = clinvar._query_sync("GATA4", "c.886G>A")
        assert res["found"] is False
        assert res["matched_on_mane"] is False


def test_pm5_finds_same_residue_candidates_under_supplied_transcript_numbering():
    """Proband ``NM_002052.5:c.886G>A`` (p.Gly296Ser) → VEP hgvsp on the
    supplied transcript gives protein_position 296. The same residue is
    numbered 297 in ClinVar's Names, where two other Pathogenic missense
    changes sit (p.Gly297Cys, p.Gly297Arg). Matching on the MANE number
    retrieves them; the payload still reports 296 for display."""
    with _temp_db():
        ev = clinvar._pm5_evidence_sync(
            "GATA4", 296, proband_alt_aa="Ser", proband_hgvs_c="c.886G>A",
            mane_protein_position=297,
        )
        names = {c["name"] for c in ev["candidates"]}
        assert any("Gly297Cys" in n for n in names), (
            f"same-residue P/LP missense not retrieved; got {sorted(names)}"
        )
        assert any("Gly297Arg" in n for n in names)
        assert ev["protein_position"] == 296
        assert ev["matched_protein_position"] == 297
        assert ev["numbering_differs"] is True


def test_ps1_pm5_do_not_match_a_different_residue_with_the_same_number():
    """The converse failure, and the more dangerous one. At protein_position
    296 the only ClinVar Name carrying "296" is ``p.Cys296Ser`` — Cys295 in
    NM_002052.5 numbering, a DIFFERENT residue from the proband's Gly296.
    Its alt AA (Ser) equals the proband's, so it WAS emitted as a PS1
    comparison variant ("same amino-acid change, different nucleotide") and
    only its 1★ review status stopped PS1_Strong (+4) from firing off a
    variant at another residue. Anchoring the match to MANE removes it."""
    with _temp_db():
        ev = clinvar._pm5_evidence_sync(
            "GATA4", 296, proband_alt_aa="Ser", proband_hgvs_c="c.886G>A",
            mane_protein_position=297,
        )
        surfaced = {c["name"] for c in ev["candidates"]} | {
            c["name"] for c in ev["ps1_candidates"]
        }
        assert not any("Cys296Ser" in n for n in surfaced), (
            "p.Cys296Ser (NM_001308093.3 numbering = residue 295 of "
            f"NM_002052.5) surfaced as evidence for residue 296: {sorted(surfaced)}"
        )


def test_mirror_holds_the_same_codon_records_at_correct_genomic_coordinates():
    """Control. The two records missing from the upstream view ARE in the mirror, at the
    right genomic coordinates and inside the gene landscape payload. The gap
    is retrieval/numbering, not a stale or incomplete mirror."""
    with _temp_db():
        ls = clinvar._gene_variant_landscape_sync("GATA4", None)
        by_pos = {p["gpos"]: p for p in ls["positions"]}
        assert 11750234 in by_pos
        assert by_pos[11750234]["csq"] == "missense"
        assert by_pos[11750234]["tier"] == "VUS"
        assert 11750236 in by_pos
        assert by_pos[11750236]["csq"] == "synonymous"
        assert by_pos[11750236]["tier"] == "VUS"


_VEP_GATA4 = {
    "ok": True,
    "most_severe_consequence": "missense_variant",
    "seq_region_name": "8",
    "start": 11750234,
    "end": 11750234,
    "transcript_strand": 1,
    "transcript_id": "NM_002052.5",
    "hgvsc": "NM_002052.5:c.907G>T",
    "hgvsp": "NP_002043.2:p.Gly303Trp",
    "cds_start": 907,
    "protein_start": 303,
    "is_mane_select": False,
    "is_mane_clinical": False,
    "transcript_consequences_all": [
        {"transcript_id": "NM_001308093.3", "hgvsc": "NM_001308093.3:c.910G>T",
         "hgvsp": "NP_001295022.1:p.Gly304Trp", "is_mane_select": True,
         "is_mane_plus_clinical": False, "is_picked": False,
         "mane_select_accession": "ENST00000532059.6"},
        {"transcript_id": "NM_002052.5", "hgvsc": "NM_002052.5:c.907G>T",
         "hgvsp": "NP_002043.2:p.Gly303Trp", "is_mane_select": False,
         "is_mane_plus_clinical": False, "is_picked": True},
        {"transcript_id": "NM_001374274.1", "hgvsc": "NM_001374274.1:c.165+1149G>T",
         "hgvsp": None, "is_mane_select": False,
         "is_mane_plus_clinical": False, "is_picked": False},
    ],
}


def test_mane_extractors_read_the_mane_row_not_the_picked_row():
    from backend.acmg.hard_coded import (
        _any_mane_protein_position_from_vep,
        _mane_hgvsc_from_vep,
        _mane_hgvsp_from_vep,
    )
    from backend.evidence import (
        _mane_protein_position_from_vep,
        _protein_position_from_vep,
    )
    assert _mane_hgvsp_from_vep(_VEP_GATA4) == "NP_001295022.1:p.Gly304Trp"
    assert _mane_hgvsc_from_vep(_VEP_GATA4) == "c.910G>T"
    assert _protein_position_from_vep(_VEP_GATA4) == 303
    assert _mane_protein_position_from_vep(_VEP_GATA4) == 304
    assert _any_mane_protein_position_from_vep(_VEP_GATA4) == 304


def test_mane_extractors_return_none_when_no_mane_annotation():
    """No MANE row and the picked transcript isn't MANE → None, so callers
    fall back to the curator's numbering instead of guessing an offset."""
    from backend.acmg.hard_coded import _mane_hgvsc_from_vep, _mane_hgvsp_from_vep
    from backend.evidence import _mane_protein_position_from_vep
    vep = dict(_VEP_GATA4)
    vep["transcript_consequences_all"] = [
        {"transcript_id": "NM_002052.5", "hgvsc": "NM_002052.5:c.907G>T",
         "hgvsp": "NP_002043.2:p.Gly303Trp", "is_mane_select": False,
         "is_mane_plus_clinical": False, "is_picked": True},
    ]
    assert _mane_hgvsp_from_vep(vep) is None
    assert _mane_hgvsc_from_vep(vep) is None
    assert _mane_protein_position_from_vep(vep) is None


def test_mane_extractors_use_picked_transcript_when_it_is_itself_mane():
    """Curator supplied the MANE transcript (or a trimmed response omitted
    the table): the picked row IS the MANE row, offset zero."""
    from backend.acmg.hard_coded import _mane_hgvsc_from_vep
    from backend.evidence import _mane_protein_position_from_vep
    vep = {
        "ok": True, "most_severe_consequence": "missense_variant",
        "hgvsc": "NM_001308093.3:c.910G>T",
        "hgvsp": "NP_001295022.1:p.Gly304Trp",
        "is_mane_select": True,
    }
    assert _mane_hgvsc_from_vep(vep) == "c.910G>T"
    assert _mane_protein_position_from_vep(vep) == 304


def test_same_position_records_are_not_counted_as_finding_the_variant():
    """The variant under test, GATA4 NM_002052.5:c.907G>T (= MANE c.910G>T),
    is genuinely absent from ClinVar (verified live 2026-08-24: esearch
    "GATA4[gene] AND c.910G>T" → Count 0). The G>A allele at the SAME base
    IS in ClinVar. Position matching must surface it as context and must NOT
    flip ``found`` to True — the mirror's allele columns are 'na', so the
    allele cannot be confirmed."""
    with _temp_db():
        res = clinvar._query_sync(
            "GATA4", "c.907G>T", "c.910G>T",
            chrom="8", pos=11750234, end=11750234,
        )
        assert res["found"] is False
        assert res["position_matching_used"] is True
        titles = [r["title"] for r in res["same_position_records"]]
        assert any("c.910G>A (p.Gly304Arg)" in t for t in titles), titles
        assert all(
            r["allele_confirmed"] is False for r in res["same_position_records"]
        )
        assert res["same_position_allele_confirmed"] is False


def test_position_matching_is_skipped_for_multi_base_loci():
    """Indel guard. The mirror stores ClinVar's ``Start``, which equals
    ``PositionVCF`` for 499,997/500,000 sampled GRCh38 SNV rows but disagrees
    for 65.1% of non-SNV rows, so position matching is single-base only."""
    with _temp_db():
        res = clinvar._query_sync(
            "GATA4", "c.907_909del", None,
            chrom="8", pos=11750234, end=11750236,
        )
        assert res["position_matching_used"] is False
        assert res["same_position_records"] == []


def test_same_residue_surface_labels_synonymous_and_vus_as_pm5_ineligible():
    """The two records from the original report. Both must be VISIBLE (that
    was the gap) and both must be explicitly PM5-INELIGIBLE, for different
    reasons that must not be conflated:
      * p.Gly304= is Gly→Gly — synonymous, so it cannot support PM5 under
        any numbering, regardless of how it is classified.
      * p.Gly304Arg IS a different missense at the residue, but it is a VUS,
        and PM5 needs an established P/LP comparator.
    """
    from backend.clients.clinvar import _same_residue_sync
    with _temp_db():
        res = _same_residue_sync(
            "GATA4", 304, 303, proband_alt_aa="Trp", proband_hgvs_c="c.907G>T",
        )
        assert res["numbering_differs"] is True
        assert res["display_position"] == 303 and res["matched_position"] == 304
        by_name = {r["name"]: r for r in res["records"]}
        syn = next(r for n, r in by_name.items() if "p.Gly304=" in n)
        mis = next(r for n, r in by_name.items() if "p.Gly304Arg" in n)
        assert syn["consequence"] == "synonymous"
        assert syn["pm5_eligible"] is False
        assert "MISSENSE" in syn["pm5_ineligible_reason"]
        assert mis["consequence"] == "missense"
        assert mis["pm5_eligible"] is False
        assert "Uncertain significance" in mis["pm5_ineligible_reason"]
        assert res["pm5_eligible_count"] == 0


def test_same_residue_surface_marks_a_genuine_pm5_comparator_eligible():
    """Control in the other direction. Proband p.Gly297Trp (MANE numbering):
      * p.Gly297Ser is a DIFFERENT missense at the residue, Pathogenic/Likely
        pathogenic at 2★ — a genuine PM5 comparator, must be eligible.
      * p.Gly297Cys/Arg are Pathogenic but 0★ ("no assertion criteria
        provided"), so the reason must be the star bar specifically."""
    from backend.clients.clinvar import _same_residue_sync
    with _temp_db():
        res = _same_residue_sync(
            "GATA4", 297, 296, proband_alt_aa="Trp", proband_hgvs_c="c.886G>T",
        )
        by = {r["name"]: r for r in res["records"]}
        ser = next(r for n, r in by.items() if "Gly297Ser" in n)
        assert ser["pm5_eligible"] is True, ser["pm5_ineligible_reason"]
        assert res["pm5_eligible_count"] == 1
        cys = next(r for n, r in by.items() if "Gly297Cys" in n)
        assert cys["pm5_eligible"] is False
        assert "2★" in cys["pm5_ineligible_reason"], cys["pm5_ineligible_reason"]


def test_same_residue_surface_calls_the_same_substitution_ps1_territory():
    """Proband p.Gly297Ser: the ClinVar p.Gly297Ser record encodes the SAME
    substitution, so it is PS1 territory (or the proband's own record) and
    must never be offered as PM5 support."""
    from backend.clients.clinvar import _same_residue_sync
    with _temp_db():
        res = _same_residue_sync(
            "GATA4", 297, 296, proband_alt_aa="Ser", proband_hgvs_c="c.886G>A",
        )
        ser = next(r for r in res["records"] if "Gly297Ser" in r["name"])
        assert ser["pm5_eligible"] is False
        assert "PS1 territory" in ser["pm5_ineligible_reason"]


def test_codon_positions_from_cds_start_and_strand():
    """Verified against the live VEP response for NM_002052.5:c.907G>T:
    cds_start 907 → offset 0, strand +1, start 11750234 → the codon spans
    11750234-11750236, which is MANE c.910-912 (codon 304)."""
    from backend.clients.ensembl_vep import codon_genomic_positions
    exons = [{"rank": 6, "start": 11750100, "end": 11750236},
             {"rank": 7, "start": 11751000, "end": 11751200}]
    vep = {"ok": True, "cds_start": 907, "start": 11750234,
           "end": 11750234, "transcript_strand": 1}
    assert codon_genomic_positions(vep, exons) == [11750234, 11750235, 11750236]
    vep3 = dict(vep, cds_start=909, start=11750236, end=11750236)
    assert codon_genomic_positions(vep3, exons) == [11750234, 11750235, 11750236]
    rev = {"ok": True, "cds_start": 907, "start": 11750234,
           "end": 11750234, "transcript_strand": -1}
    rev_exons = [{"rank": 1, "start": 11750100, "end": 11750300}]
    assert codon_genomic_positions(rev, rev_exons) == [11750232, 11750233, 11750234]


def test_codon_positions_use_the_transcript_strand_not_the_allele_strand():
    """⚠ THE BUG THIS FIELD EXISTS FOR, and why GATA4 could never catch it.

    ``strand`` on the VEP result is the ALLELE-REPRESENTATION strand and means
    different things per endpoint — measured for MYH7 R403Q, a minus-strand
    gene: /vep/human/region gives +1 with C/T (the VCF row), /vep/human/hgvs
    gives -1 with G/A (the transcript). So on the COORDINATE path this function
    used to receive +1 for every gene regardless of orientation.

    MYH7 at 14:23429278, reproduced exactly:
        c.1207 (offset 0): correct 277,278,279   allele-strand gave 279,280,281
        c.1208 (offset 1): identical  <- the only one that coincides
        c.1209 (offset 2): correct 277,278,279   allele-strand gave 275,276,277

    The docstring's original verification used GATA4, a PLUS-strand gene, where
    the two conventions agree — which is exactly why it survived."""
    from backend.clients.ensembl_vep import codon_genomic_positions
    exons = [{"rank": 1, "start": 23429000, "end": 23429500}]
    for cds_start, pos in ((1207, 23429279), (1209, 23429277)):
        vep = {"ok": True, "cds_start": cds_start, "start": pos, "end": pos,
               "transcript_strand": -1}
        assert codon_genomic_positions(vep, exons) == [23429277, 23429278, 23429279]

    vep = {"ok": True, "cds_start": 1207, "start": 23429279, "end": 23429279,
           "transcript_strand": -1, "strand": 1}
    assert codon_genomic_positions(vep, exons) == [23429277, 23429278, 23429279]


def test_codon_positions_decline_without_the_transcript_strand():
    """Required, not defaulted. A field that means two different things cannot
    be trusted as the CDS direction, and declining costs only the codon view —
    the caller falls back to the variant's own position — whereas guessing
    produces a plausible, wrong codon."""
    from backend.clients.ensembl_vep import codon_genomic_positions
    exons = [{"rank": 1, "start": 11750100, "end": 11750300}]
    vep = {"ok": True, "cds_start": 907, "start": 11750234, "end": 11750234,
           "strand": 1}
    assert codon_genomic_positions(vep, exons) is None


def test_codon_positions_handles_intron_straddling_codon():
    """Codon whose first base is the last base of an exon: the remaining two
    bases come from the start of the next exon, not from the intron."""
    from backend.clients.ensembl_vep import codon_genomic_positions
    exons = [{"rank": 1, "start": 1000, "end": 1100},
             {"rank": 2, "start": 2000, "end": 2100}]
    vep = {"ok": True, "cds_start": 100, "start": 1100, "end": 1100,
           "transcript_strand": 1}
    assert codon_genomic_positions(vep, exons) == [1100, 2000, 2001]


def test_codon_positions_declines_without_an_exon_ladder_or_for_indels():
    """Two explicit declines. Without the exon ladder we cannot PROVE the
    codon is contiguous, and an indel has no single anchoring codon. Both
    return None, which callers must read as 'no codon view available'."""
    from backend.clients.ensembl_vep import codon_genomic_positions
    vep = {"ok": True, "cds_start": 907, "start": 11750234,
           "end": 11750234, "transcript_strand": 1}
    assert codon_genomic_positions(vep, None) is None
    indel = dict(vep, end=11750236)
    exons = [{"rank": 1, "start": 11750100, "end": 11750300}]
    assert codon_genomic_positions(indel, exons) is None


def test_same_site_frequencies_degrade_gracefully_off_panel(monkeypatch):
    """The local gnomAD frequency DB covers the cardiac panel intervals only.
    Off-panel must say so rather than firing a live per-base query."""
    import asyncio
    import backend.clients.gnomad as gnomad
    monkeypatch.setattr(gnomad, "_freq_db_exists", lambda: True)
    monkeypatch.setattr(gnomad, "_in_panel", lambda vid: False)
    res = asyncio.run(gnomad.fetch_same_site_frequencies("8", [11750234]))
    assert res["ok"] is True
    assert res["available"] is False
    assert res["reason"] == "off-panel"
    assert "off-panel" in res["note"]


def test_ps1_self_exclusion_recognises_the_proband_in_mane_coordinates():
    from backend.acmg.hard_coded import _clinvar_ps1_criterion
    with _temp_db():
        ev = clinvar._pm5_evidence_sync(
            "GATA4", 296, proband_alt_aa="Ser", proband_hgvs_c="c.886G>A",
            mane_protein_position=297, proband_hgvs_c_mane="c.889G>A",
        )
        names = [c["name"] for c in ev["ps1_candidates"]]
        assert names == [], f"proband's own record offered as a PS1 comparator: {names}"
        assert _clinvar_ps1_criterion(ev, "GATA4")["status"] == "not_met"


def test_ps1_self_exclusion_without_the_mane_token_is_the_bug():
    """Pinned so the wiring cannot silently regress: drop the MANE token and
    the proband's own 2★ record comes back as a PS1 comparator."""
    from backend.acmg.hard_coded import _clinvar_ps1_criterion
    with _temp_db():
        ev = clinvar._pm5_evidence_sync(
            "GATA4", 296, proband_alt_aa="Ser", proband_hgvs_c="c.886G>A",
            mane_protein_position=297,
        )
        assert [c["name"] for c in ev["ps1_candidates"]] == [
            "NM_001308093.3(GATA4):c.889G>A (p.Gly297Ser)"
        ]
        assert _clinvar_ps1_criterion(ev, "GATA4")["criteria_strength"] == "PS1_Strong"


def test_pm5_fires_for_a_genuinely_novel_missense_at_the_residue():
    """The positive control the whole fix exists for. Proband p.Gly296Cys
    (NM_002052.5 c.886G>T = MANE c.889G>T): the 2★ Pathogenic/Likely
    pathogenic p.Gly297Ser at the same residue is a different missense and
    is NOT the proband, so PM5 must fire — and the proband's own c.889G>T
    record must be excluded."""
    from backend.acmg.hard_coded import _clinvar_pm5_criterion
    with _temp_db():
        ev = clinvar._pm5_evidence_sync(
            "GATA4", 296, proband_alt_aa="Cys", proband_hgvs_c="c.886G>T",
            mane_protein_position=297, proband_hgvs_c_mane="c.889G>T",
        )
        crit = _clinvar_pm5_criterion(ev, "GATA4")
        assert crit["status"] == "met"
        assert crit["criteria_strength"] == "PM5_Moderate"
        assert "Gly297Ser" in crit["evidence"]
        assert "residue 296 (= residue 297 on ClinVar's MANE transcript)" in crit["evidence"]
        assert "Gly297Cys" not in [c["name"] for c in ev["candidates"]][0]


def test_same_residue_surface_labels_the_proband_own_record_as_such():
    from backend.clients.clinvar import _same_residue_sync
    with _temp_db():
        res = _same_residue_sync(
            "GATA4", 297, 296, proband_alt_aa="Ser", proband_hgvs_c="c.886G>A",
            proband_hgvs_c_mane="c.889G>A",
        )
        own = next(r for r in res["records"] if "Gly297Ser" in r["name"])
        assert own["is_proband_own_record"] is True
        assert own["pm5_eligible"] is False
        assert "own record" in own["pm5_ineligible_reason"]


def test_same_site_excludes_the_probands_own_allele(monkeypatch):
    """VEP REST returns ``vcf_string`` as a LIST (``["8-11750234-G-T"]``), so
    comparing allele ids against it directly never matched and the proband's
    own allele was listed under "other alleles at this base" with AC 0 (D7).
    Exercised through build_same_site_evidence with the gnomAD and ClinVar
    lookups stubbed, so this tests the exclusion and nothing else."""
    import asyncio

    import backend.evidence as evidence

    async def _fake_freqs(chrom, positions):
        return {
            "ok": True, "available": True, "chrom": chrom,
            "positions": positions,
            "alleles": [
                {"variant_id": "8-11750234-G-T", "position": 11750234,
                 "ref": "G", "alt": "T", "rsid": None, "ac": 0, "an": 1460238,
                 "af": 0.0, "faf95_popmax": None},
                {"variant_id": "8-11750234-G-A", "position": 11750234,
                 "ref": "G", "alt": "A", "rsid": "rs1205549216", "ac": 15,
                 "an": 1612518, "af": 9.3e-06, "faf95_popmax": 5.75e-06},
                {"variant_id": "8-11750236-G-A", "position": 11750236,
                 "ref": "G", "alt": "A", "rsid": "rs773684507", "ac": 2,
                 "an": 1460216, "af": 1.37e-06, "faf95_popmax": 7.41e-06},
            ],
            "count": 3,
        }

    async def _no_records(*a, **k):
        return {"ok": True, "records": [], "count": 0, "pm5_eligible_count": 0,
                "display_position": 303, "matched_position": 304,
                "numbering_differs": True}

    async def _ladder(_tx):
        return {"ok": True, "transcript_strand": 1,
                "exons": [{"rank": 6, "start": 11750100, "end": 11750236}]}

    monkeypatch.setattr(evidence, "fetch_same_site_frequencies", _fake_freqs)
    monkeypatch.setattr(evidence, "get_same_residue_records", _no_records)
    monkeypatch.setattr(evidence, "fetch_transcript_exons", _ladder)

    vep = dict(_VEP_GATA4)
    vep["vcf_string"] = ["8-11750234-G-T"]
    ss = asyncio.run(evidence.build_same_site_evidence({"vep": vep}, "GATA4", None))

    same_nt = {a["variant_id"] for a in ss["same_nucleotide"]}
    assert "8-11750234-G-T" not in same_nt, (
        f"the proband's own allele is listed as an 'other' allele: {same_nt}"
    )
    assert same_nt == {"8-11750234-G-A"}
    assert [a["variant_id"] for a in ss["same_codon"]] == ["8-11750236-G-A"]
    assert ss["codon_positions"] == [11750234, 11750235, 11750236]
    assert ss["codon_span_confirmed_on_picked_transcript"] is False
    assert ss["codon_span_source"] == "MANE ENST00000532059"
    assert "(" not in ss["codon_span_source"], (
        "the source token is a UI label — nested parentheses in it blew the "
        "row out to five lines"
    )
    assert ss["codon_span_note"] and "unversioned ENST" in ss["codon_span_note"]


def test_ps1_row_states_that_the_probands_own_record_was_excluded():
    """A residue where the only same-AA record IS the proband read identically
    to a residue with nothing submitted, so the guard that stops PS1 firing off
    the variant's own ClinVar record left no trace on the criterion a curator
    audits."""
    from backend.acmg.hard_coded import _clinvar_ps1_criterion
    with _temp_db():
        ev = clinvar._pm5_evidence_sync(
            "GATA4", 296, proband_alt_aa="Ser", proband_hgvs_c="c.886G>A",
            mane_protein_position=297, proband_hgvs_c_mane="c.889G>A",
        )
        assert ev["ps1_self_excluded_count"] == 1
        crit = _clinvar_ps1_criterion(ev, "GATA4")
        assert crit["status"] == "not_met"
        assert "proband's OWN ClinVar record" in crit["evidence"]
        assert "c.889G>A" in crit["evidence"]
        assert "PS1 not_met." not in crit["evidence"]
