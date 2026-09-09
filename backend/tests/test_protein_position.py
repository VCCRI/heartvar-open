"""Regression coverage for the VEP protein-position extractors.

These three helpers feed PM1 (domain query, any protein-altering consequence)
and PM5 (missense-only). They were silently broken once during the app.py ->
acmg/ package split when `_hgvsp_int_position` lost its `return` line — a defect
invisible to the rest of the suite and the golden harness because no other
fixture threads a VEP ``hgvsp`` through them. This file closes that blind spot.

Imported via ``backend.app`` to also exercise the re-export surface (the split
re-imports these from ``backend.acmg.hard_coded``).
"""
from backend.app import (
    _hgvsp_int_position,
    _any_protein_position_from_vep,
    _protein_position_from_vep,
)


class TestHgvspIntPosition:
    def test_parses_residue_position(self):
        assert _hgvsp_int_position("NP_004324.2:p.Phe595Leu") == 595

    def test_takes_first_number_after_colon_split(self):
        assert _hgvsp_int_position("ENSP00000123456.1:p.Arg403Gln") == 403

    def test_none_when_no_digits_in_suffix(self):
        assert _hgvsp_int_position("NP_004324.2:p.?") is None

    def test_none_on_empty(self):
        assert _hgvsp_int_position("") is None


def _vep(consequence, hgvsp="NP_004324.2:p.Phe595Leu", ok=True):
    return {"ok": ok, "most_severe_consequence": consequence, "hgvsp": hgvsp}


class TestAnyProteinPositionFromVep:
    """PM1 extractor — accepts any protein-altering consequence."""

    def test_missense(self):
        assert _any_protein_position_from_vep(_vep("missense_variant")) == 595

    def test_stop_gain(self):
        assert _any_protein_position_from_vep(_vep("stop_gained")) == 595

    def test_inframe_deletion(self):
        assert _any_protein_position_from_vep(_vep("inframe_deletion")) == 595

    def test_skips_synonymous(self):
        assert _any_protein_position_from_vep(_vep("synonymous_variant")) is None

    def test_skips_splice_region(self):
        assert _any_protein_position_from_vep(_vep("splice_region_variant")) is None

    def test_skips_intron(self):
        assert _any_protein_position_from_vep(_vep("intron_variant")) is None

    def test_none_when_vep_not_ok(self):
        assert _any_protein_position_from_vep(_vep("missense_variant", ok=False)) is None

    def test_none_when_no_hgvsp(self):
        assert _any_protein_position_from_vep(_vep("missense_variant", hgvsp="")) is None
