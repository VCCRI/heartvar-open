"""UniProt surfaces must be read in UniProt's own numbering, not the curator's.

Same class of bug as the ClinVar cross-transcript residue mismatch. UniProt
annotates the CANONICAL isoform — for the panel genes, the MANE Select protein.
Three places compared a residue number from the transcript the CURATOR supplied
against UniProt-frame annotation:

  1. `_build_domain_plp_evidence`: `_smallest_containing_domain(domains, aa_pos)`
     — the decision that gates PM1 entirely. A wrong-frame residue landing
     outside every domain returns `in_domain: False` and PM1 is never assessed.
     The MANE residue was already being computed 16 lines further down, purely
     to offset the ClinVar query, so the containment test stayed in the wrong
     frame while the ClinVar window was corrected.
  2. `protein_position_for_uniprot` in the evidence gather: windows UniProt's
     natural-variant list around the wrong residue.
  3. `_bp3_functionless_repeat`: locates the UniProt repeat feature. BP3 is a
     BENIGN criterion, so a wrong-frame hit is a false-benign.

The worked example throughout is GATA4, whose real isoform offset the codebase
already documents: MANE `NM_001308093.3` / `NP_001295022.1` is 443 aa,
`NM_002052.5` / `NP_002043.2` is 442 aa (one extra Val at position 206 in the
MANE isoform), so every residue past 205 is numbered +1 in MANE. Curator residue
303 and MANE residue 304 are the SAME residue.
"""
from __future__ import annotations

import asyncio

import backend.evidence as evidence
from backend.acmg.hard_coded import _bp3_functionless_repeat


def _gata4_vep(curator_pos=303, mane_pos=304, csq="missense_variant"):
    """VEP result for GATA4 entered on NM_002052.5, with the MANE row present."""
    return {
        "ok": True,
        "most_severe_consequence": csq,
        "hgvsp": f"NP_002043.2:p.Gly{curator_pos}Trp",
        "transcript_consequences_all": [
            {"transcript_id": "NM_002052.5",
             "hgvsp": f"NP_002043.2:p.Gly{curator_pos}Trp"},
            {"transcript_id": "NM_001308093.3", "is_mane_select": True,
             "hgvsp": f"NP_001295022.1:p.Gly{mane_pos}Trp",
             "hgvsc": "c.910G>T"},
        ],
    }


_UNIPROT = {
    "ok": True,
    "found": True,
    "length": 443,
    "domains": [{"name": "GATA-type 2", "start": 304, "end": 349}],
    "features": [],
    "natural_variants": [],
}


def test_domain_containment_is_tested_in_uniprot_numbering(monkeypatch):
    """THE test that would have caught the bug. Curator residue 303 is outside
    the domain; MANE residue 304 is inside it. Pre-fix the test used 303, so
    `in_domain` was False and PM1 was never assessed for this variant."""
    captured = {}

    async def _fake_domain_plp(gene, start, end, exclude_position=None,
                              numbering_offset=0):
        captured.update(gene=gene, start=start, end=end,
                        exclude_position=exclude_position,
                        numbering_offset=numbering_offset)
        return {"ok": True, "plp_count": 4, "max_stars": 3, "all_hits": []}

    monkeypatch.setattr(evidence, "get_domain_plp_evidence", _fake_domain_plp)
    out = asyncio.run(evidence._build_domain_plp_evidence(
        "GATA4", _gata4_vep(), _UNIPROT))

    assert out.get("in_domain") is True, (
        "curator residue 303 vs a UniProt domain starting at 304: the domain "
        "test must run in UniProt's frame, not the curator's"
    )
    assert out["domain_name"] == "GATA-type 2"
    assert out["variant_position"] == 303
    assert out["matched_protein_position"] == 304
    assert out["numbering_differs"] is True
    assert captured["exclude_position"] == 304
    assert captured["numbering_offset"] == 0, (
        "double-shifting the query would move the window off the domain again"
    )


def test_curator_numbering_is_kept_when_there_is_no_mane_row(monkeypatch):
    """Fail-back, not fail-guess. With no MANE annotation the behaviour is the
    unchanged pre-fix one."""
    async def _fake(gene, start, end, exclude_position=None, numbering_offset=0):
        return {"ok": True, "plp_count": 1, "all_hits": [],
                "_excl": exclude_position}

    monkeypatch.setattr(evidence, "get_domain_plp_evidence", _fake)
    vep = {"ok": True, "most_severe_consequence": "missense_variant",
           "hgvsp": "NP_002043.2:p.Gly310Trp"}
    up = dict(_UNIPROT, domains=[{"name": "d", "start": 300, "end": 350}])
    out = asyncio.run(evidence._build_domain_plp_evidence("GATA4", vep, up))
    assert out["variant_position"] == 310
    assert out["matched_protein_position"] == 310
    assert out["numbering_differs"] is False


def test_mane_position_off_the_end_of_the_uniprot_sequence_is_rejected(monkeypatch):
    """Guard on the guard. If the MANE residue falls past the end of the UniProt
    sequence but the curator's does not, UniProt's canonical is NOT the MANE
    isoform for this gene, so the MANE number is the wrong frame. Keep the
    curator's rather than querying a residue that cannot exist."""
    async def _fake(gene, start, end, exclude_position=None, numbering_offset=0):
        return {"ok": True, "plp_count": 1, "all_hits": []}

    monkeypatch.setattr(evidence, "get_domain_plp_evidence", _fake)
    up = dict(_UNIPROT, length=310,
              domains=[{"name": "d", "start": 300, "end": 310}])
    out = asyncio.run(evidence._build_domain_plp_evidence(
        "GATA4", _gata4_vep(curator_pos=305, mane_pos=999), up))
    assert out["matched_protein_position"] == 305
    assert out["numbering_differs"] is False


def test_natural_variant_window_helper_prefers_the_mane_residue():
    assert evidence._mane_protein_position_from_vep(_gata4_vep()) == 304
    assert evidence._protein_position_from_vep(_gata4_vep()) == 303


def test_bp3_repeat_lookup_uses_the_mane_residue():
    """The repeat spans MANE 304-320 only. On the curator's 303 the variant
    reads as outside the repeat; in UniProt's frame it is inside. Because BP3 is
    benign evidence, getting this frame wrong in the other direction asserts a
    functionless-repeat variant that is not in the repeat at all."""
    ev = {"uniprot": {"ok": True, "features": [
        {"type": "Repeat", "start": 304, "end": 320, "description": "poly-Gln"},
    ]}}
    qualifies, reason = _bp3_functionless_repeat(ev, _gata4_vep())
    assert "304" in reason or qualifies, (
        f"repeat lookup did not run in UniProt numbering: {reason}"
    )


def test_bp3_still_fails_closed_without_a_position():
    ev = {"uniprot": {"ok": True, "features": []}}
    qualifies, reason = _bp3_functionless_repeat(ev, {"ok": True})
    assert qualifies is False
    assert "fail-closed" in reason
