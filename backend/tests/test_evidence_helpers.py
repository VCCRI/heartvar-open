"""Behavioural unit tests for the pure helper functions in backend.evidence.

These helpers (the ClinVar-landscape keyword builder, the multi-step GenCC
disease normaliser, the cardiac-root substring filter, the PM1 domain-context
composers, and the missense AA/position extractors) had ZERO direct coverage.
A prior app.py -> acmg/ package split silently broke one such helper
(`_hgvsp_int_position` lost its `return`) and nothing in the suite caught it,
because no other fixture threads these inputs through the helpers. This file
closes that blind spot.

Imports are taken DIRECTLY from ``backend.evidence`` (not the ``backend.app``
re-export shim) so the tests stay decoupled from the legacy re-export surface
and pin the behaviour at the module that actually owns it.
"""
from __future__ import annotations

from backend import evidence
from backend.evidence import (
    build_cardiac_keywords,
    _normalise_gencc_disease,
    _contains_cardiac_root,
    _smallest_containing_domain,
    _describe_outside_domain_context,
    _is_unstructured_region,
    _proband_aa_from_vep,
    _PROBAND_AA_RE,
    _protein_position_from_vep,
)


class TestContainsCardiacRoot:
    def test_matches_anatomical_root(self):
        assert _contains_cardiac_root("Hypertrophic cardiomyopathy") is True

    def test_matches_partial_root_substring(self):
        assert _contains_cardiac_root("Bicuspid aortic valve") is True

    def test_matches_named_chd_syndrome_without_cardiac_word(self):
        assert _contains_cardiac_root("Noonan syndrome") is True
        assert _contains_cardiac_root("CHARGE association") is True

    def test_case_insensitive(self):
        assert _contains_cardiac_root("HYPERTROPHIC CARDIOMYOPATHY") is True

    def test_rejects_non_cardiac(self):
        assert _contains_cardiac_root("Diabetes mellitus") is False

    def test_rejects_empty_string(self):
        assert _contains_cardiac_root("") is False

    def test_rejects_none(self):
        assert _contains_cardiac_root(None) is False


class TestNormaliseGenccDisease:
    def test_blank_input_returns_empty(self):
        assert _normalise_gencc_disease("") == []
        assert _normalise_gencc_disease("   ") == []

    def test_lowercases(self):
        assert _normalise_gencc_disease("Cardiomyopathy") == ["cardiomyopathy"]

    def test_step_b_strips_parenthetical_suffix(self):
        assert _normalise_gencc_disease("cardiomyopathy (hypertrophic)") == [
            "cardiomyopathy"
        ]

    def test_step_c_strips_with_or_without_tail(self):
        assert _normalise_gencc_disease(
            "testicular anomalies with or without congenital heart disease"
        ) == ["testicular anomalies"]

    def test_step_d_strips_type_suffix_singular(self):
        assert _normalise_gencc_disease("Long QT syndrome type 1") == [
            "long qt syndrome"
        ]

    def test_step_d_strips_type_suffix_plural(self):
        assert _normalise_gencc_disease("Long QT syndrome types 1") == [
            "long qt syndrome"
        ]

    def test_step_e_strips_gene_suffix_hyphen(self):
        assert _normalise_gencc_disease("structural heart disease - gata4") == [
            "structural heart disease"
        ]

    def test_step_e_strips_gene_suffix_en_dash(self):
        assert _normalise_gencc_disease("structural heart disease – GATA4") == [
            "structural heart disease"
        ]

    def test_step_f_splits_on_and(self):
        assert _normalise_gencc_disease(
            "atrial septal defect and ventricular septal defect"
        ) == ["atrial septal defect", "ventricular septal defect"]

    def test_step_f_splits_on_comma(self):
        assert _normalise_gencc_disease("noonan syndrome, leopard syndrome") == [
            "noonan syndrome",
            "leopard syndrome",
        ]

    def test_step_g_drops_parts_below_min_length(self):
        assert _normalise_gencc_disease("ASD") == []

    def test_step_g_keeps_long_drops_short_in_compound(self):
        assert _normalise_gencc_disease("cardiomyopathy, cat") == ["cardiomyopathy"]

    def test_dedup_is_case_insensitive_and_order_preserving(self):
        assert _normalise_gencc_disease("Cardiomyopathy and cardiomyopathy") == [
            "cardiomyopathy"
        ]


class TestBuildCardiacKeywords:
    def test_empty_inputs_yield_only_gene_fallbacks(self):
        assert build_cardiac_keywords(None, None, "GATA4") == [
            "gata4-related disorder",
            "gata4 syndrome",
            "cardiovascular phenotype",
        ]

    def test_empty_lists_behave_like_none(self):
        assert build_cardiac_keywords([], [], "GATA4") == [
            "gata4-related disorder",
            "gata4 syndrome",
            "cardiovascular phenotype",
        ]

    def test_no_gene_symbol_yields_only_stub(self):
        assert build_cardiac_keywords([], [], "") == ["cardiovascular phenotype"]

    def test_keeps_cardiac_phenotype_strings_drops_non_cardiac(self):
        out = build_cardiac_keywords(
            ["atrioventricular septal defect 4", "developmental delay"],
            None,
            "GATA4",
        )
        assert "atrioventricular septal defect 4" in out
        assert "developmental delay" not in out

    def test_normalises_and_filters_gencc_disease_names(self):
        out = build_cardiac_keywords(
            None,
            ["congenital heart disease (structural)", "intellectual disability"],
            "GATA4",
        )
        assert "congenital heart disease" in out
        assert "intellectual disability" not in out

    def test_dedup_across_sources_and_fallbacks(self):
        out = build_cardiac_keywords(
            ["hypertrophic cardiomyopathy"],
            ["hypertrophic cardiomyopathy (familial)"],
            "MYH7",
        )
        assert out.count("hypertrophic cardiomyopathy") == 1
        assert out[-3:] == [
            "myh7-related disorder",
            "myh7 syndrome",
            "cardiovascular phenotype",
        ]

    def test_all_keywords_lowercased(self):
        out = build_cardiac_keywords(["Hypertrophic Cardiomyopathy"], None, "MYH7")
        assert all(k == k.lower() for k in out)


class TestSmallestContainingDomain:
    def _domains(self):
        return [
            {"name": "Kinase", "start": 100, "end": 400},
            {"name": "SH3", "start": 150, "end": 200},
            {"name": "PH", "start": 500, "end": 600},
        ]

    def test_picks_tightest_nested_domain(self):
        d = _smallest_containing_domain(self._domains(), 175)
        assert d is not None and d["name"] == "SH3"

    def test_picks_only_containing_domain(self):
        d = _smallest_containing_domain(self._domains(), 350)
        assert d is not None and d["name"] == "Kinase"

    def test_none_when_outside_all_domains(self):
        assert _smallest_containing_domain(self._domains(), 450) is None

    def test_boundaries_are_inclusive(self):
        dom = [{"name": "D", "start": 10, "end": 20}]
        assert _smallest_containing_domain(dom, 10)["name"] == "D"
        assert _smallest_containing_domain(dom, 20)["name"] == "D"

    def test_none_when_aa_pos_is_none(self):
        assert _smallest_containing_domain(self._domains(), None) is None

    def test_none_when_aa_pos_non_positive(self):
        assert _smallest_containing_domain(self._domains(), 0) is None
        assert _smallest_containing_domain(self._domains(), -5) is None

    def test_none_when_no_domains(self):
        assert _smallest_containing_domain(None, 175) is None
        assert _smallest_containing_domain([], 175) is None

    def test_skips_domains_with_non_int_bounds(self):
        assert _smallest_containing_domain([{"name": "X", "start": None, "end": 20}], 15) is None


class TestIsUnstructuredRegion:
    def test_disordered_region_is_unstructured(self):
        assert _is_unstructured_region(
            {"type": "Region", "description": "Disordered"}
        ) is True

    def test_low_complexity_region_is_unstructured(self):
        assert _is_unstructured_region(
            {"type": "Region", "description": "Low complexity region"}
        ) is True

    def test_compositionally_biased_region_is_unstructured(self):
        assert _is_unstructured_region(
            {"type": "Region", "description": "Compositionally biased"}
        ) is True

    def test_region_type_match_is_case_insensitive(self):
        assert _is_unstructured_region(
            {"type": "region", "description": "compositionally biased"}
        ) is True

    def test_non_region_type_is_not_unstructured(self):
        assert _is_unstructured_region(
            {"type": "Domain", "description": "Disordered"}
        ) is False

    def test_region_without_unstructured_hint_is_false(self):
        assert _is_unstructured_region(
            {"type": "Region", "description": "Interaction with partner"}
        ) is False

    def test_missing_fields_default_false(self):
        assert _is_unstructured_region({}) is False


class TestDescribeOutsideDomainContext:
    def _doms(self):
        return [{"name": "A", "start": 10, "end": 50}, {"name": "B", "start": 200, "end": 300}]

    def test_inter_domain_linker(self):
        short, sentence = _describe_outside_domain_context([], self._doms(), 120)
        assert short == "in_inter_domain_linker"
        assert "linker between A (aa 10-50) and B (aa 200-300)" in sentence

    def test_n_terminal_of_first_domain(self):
        short, sentence = _describe_outside_domain_context([], self._doms(), 5)
        assert short == "n_terminal_of_first_domain"
        assert "N-terminal of the first annotated domain" in sentence

    def test_c_terminal_of_last_domain(self):
        short, sentence = _describe_outside_domain_context([], self._doms(), 400)
        assert short == "c_terminal_of_last_domain"
        assert "C-terminal of the last annotated domain" in sentence

    def test_no_domains_annotated(self):
        short, sentence = _describe_outside_domain_context([], [], 120)
        assert short == "no_domains_annotated"
        assert "no Domain-type features" in sentence

    def test_no_domains_with_none_inputs(self):
        short, _sentence = _describe_outside_domain_context(None, None, 120)
        assert short == "no_domains_annotated"

    def test_disordered_region_between_two_domains(self):
        feats = [{"type": "Region", "description": "Disordered", "start": 60, "end": 190}]
        short, sentence = _describe_outside_domain_context(feats, self._doms(), 120)
        assert short == "in_disordered_region"
        assert "disordered region (aa 60-190) between A (aa 10-50) and B (aa 200-300)" in sentence

    def test_disordered_n_terminal_region(self):
        feats = [{"type": "Region", "description": "Disordered", "start": 1, "end": 8}]
        short, sentence = _describe_outside_domain_context(feats, self._doms(), 5)
        assert short == "in_disordered_region"
        assert "disordered N-terminal region" in sentence

    def test_disordered_c_terminal_region(self):
        feats = [{"type": "Region", "description": "Disordered", "start": 310, "end": 420}]
        short, sentence = _describe_outside_domain_context(feats, self._doms(), 400)
        assert short == "in_disordered_region"
        assert "disordered C-terminal region" in sentence

    def test_disordered_with_no_flanking_domains(self):
        feats = [{"type": "Region", "description": "Disordered", "start": 1, "end": 500}]
        short, sentence = _describe_outside_domain_context(feats, [], 120)
        assert short == "in_disordered_region"
        assert "Variant lies in a disordered region (aa 1-500)" in sentence

    def test_outside_all_domains_with_no_flank(self):
        doms = [{"name": "Z", "start": "x", "end": "y"}]
        short, sentence = _describe_outside_domain_context([], doms, 120)
        assert short == "outside_all_domains"
        assert "not within an annotated UniProt domain" in sentence


def _vep(consequence, hgvsp="NP_004324.2:p.Phe595Leu", ok=True):
    return {"ok": ok, "most_severe_consequence": consequence, "hgvsp": hgvsp}


class TestProbandAaRegex:
    def test_captures_ref_pos_alt_groups(self):
        m = _PROBAND_AA_RE.search("p.Arg403Gln")
        assert m is not None
        assert m.group(1) == "Arg"
        assert m.group(2) == "403"
        assert m.group(3) == "Gln"

    def test_matches_inside_full_hgvsp_suffix(self):
        m = _PROBAND_AA_RE.search("NP_004324.2:p.Phe595Leu")
        assert m is not None and m.groups() == ("Phe", "595", "Leu")

    def test_does_not_match_nonsense_token(self):
        assert _PROBAND_AA_RE.search("p.?") is None


class TestProbandAaFromVep:
    def test_missense_returns_ref_and_alt(self):
        assert _proband_aa_from_vep(_vep("missense_variant")) == ("Phe", "Leu")

    def test_missense_with_protein_only_hgvsp(self):
        assert _proband_aa_from_vep(
            _vep("missense_variant", hgvsp="p.Arg403Gln")
        ) == ("Arg", "Gln")

    def test_non_missense_returns_none_pair(self):
        assert _proband_aa_from_vep(
            _vep("stop_gained", hgvsp="NP_004324.2:p.Arg403Ter")
        ) == (None, None)

    def test_synonymous_returns_none_pair(self):
        assert _proband_aa_from_vep(_vep("synonymous_variant")) == (None, None)

    def test_vep_not_ok_returns_none_pair(self):
        assert _proband_aa_from_vep(_vep("missense_variant", ok=False)) == (None, None)

    def test_non_dict_returns_none_pair(self):
        assert _proband_aa_from_vep(None) == (None, None)

    def test_unparseable_missense_hgvsp_returns_none_pair(self):
        assert _proband_aa_from_vep(
            _vep("missense_variant", hgvsp="NP_004324.2:p.?")
        ) == (None, None)

    def test_empty_hgvsp_returns_none_pair(self):
        assert _proband_aa_from_vep(_vep("missense_variant", hgvsp="")) == (None, None)


class TestProteinPositionFromVep:
    def test_missense_returns_position(self):
        assert _protein_position_from_vep(_vep("missense_variant")) == 595

    def test_rejects_stop_gain(self):
        assert _protein_position_from_vep(_vep("stop_gained")) is None

    def test_rejects_synonymous(self):
        assert _protein_position_from_vep(_vep("synonymous_variant")) is None

    def test_rejects_inframe_deletion(self):
        assert _protein_position_from_vep(_vep("inframe_deletion")) is None

    def test_none_when_vep_not_ok(self):
        assert _protein_position_from_vep(_vep("missense_variant", ok=False)) is None

    def test_none_when_non_dict(self):
        assert _protein_position_from_vep(None) is None

    def test_none_when_missense_but_no_position(self):
        assert _protein_position_from_vep(
            _vep("missense_variant", hgvsp="NP_004324.2:p.?")
        ) is None


def test_a_bare_coding_change_is_not_sent_to_the_recoder():
    """The only shape this function ever actually received."""
    for bare in ("c.301G>A", "c.611G>A", "c.5058del", "c.1207_1209del"):
        assert evidence._recoder_can_resolve(bare) is False, bare


def test_an_anchored_hgvs_or_an_rsid_is_still_allowed():
    """The guard is about resolvability, not about banning the endpoint — a
    caller that does have an anchor is not blocked."""
    assert evidence._recoder_can_resolve("NM_000257.4:c.1208G>A") is True
    assert evidence._recoder_can_resolve("MYH7:c.1208G>A") is True
    assert evidence._recoder_can_resolve("rs121913642") is True
    assert evidence._recoder_can_resolve("RS121913642") is True


def test_empty_input_is_not_sent_anywhere():
    assert evidence._recoder_can_resolve("") is False
    assert evidence._recoder_can_resolve(None) is False
    assert evidence._recoder_can_resolve("   ") is False
