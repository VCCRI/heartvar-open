"""Phenotype-input vocabulary: ontology coverage, acronyms, compounds.

Driven by a real curator report — a proband entered as ``DORV-TGA, AS``
came back "phenotype terms not recognised". Investigation found three
independent defects, one per section below.

1. COMPOUND TOKENS. ``DORV`` and ``TGA`` were BOTH already recognised
   individually, but joined by a hyphen they failed: `_normalise_phrase`
   turns the hyphen into a space giving "dorv tga", which is no longer an
   exact key, and the substring fallback requires a table key of >= 6
   characters, so a 4- and a 3-character acronym can never rescue it.

2. ONTOLOGY ROOTS. The cardiac descendant closure was built from five
   roots covering cardiac MORPHOLOGY, cardiomyopathy and arrhythmia. The
   whole ``HP:0011025`` cardiovascular-PHYSIOLOGY subtree was absent, so
   valve function, conduction disease, QT/repolarisation and pulmonary
   circulation matched no category. That inverted the intended design:
   free text "long QT syndrome" matched via the phrase table while
   HP:0001657 (Prolonged QT interval) — the canonical ontology term for
   the same thing — matched nothing, silently suppressing PP4 on exactly
   the channelopathy genes it should support.

3. ACRONYMS. The hand-written phrase table held 81 entries. hp.obo
   supplies 678 EXACT synonyms across the cardiac subtree for free
   (including DORV, TGA, D-TGA, ASD, TAPVR, HCM, SVT, PFO), leaving only
   a small curated top-up for clinical shorthand HPO lacks — notably TOF,
   whose only HPO synonym is the typo "Tetrology of fallot".
"""
from __future__ import annotations

import pytest

from backend.clients import panelapp
from backend.clients.panelapp import (
    _parse_phenotype_input,
    check_hpo_relevance,
    resolve_phenotype_terms,
)


@pytest.fixture(autouse=True)
def _descendants():
    panelapp._load_local_descendants()


def _cats(term, cats):
    return [c for c in cats if check_hpo_relevance(term, c)]


def test_the_reported_input_resolves_end_to_end():
    """The exact string the curator submitted, nothing left unrecognised."""
    hpo, cats, recognised, unrecognised = _parse_phenotype_input("DORV-TGA, AS")
    assert unrecognised == []
    assert len(recognised) >= 2
    assert check_hpo_relevance("DORV-TGA", "Congenital heart disease")
    assert check_hpo_relevance("AS", "Congenital heart disease")


@pytest.mark.parametrize("raw", ["DORV-TGA", "DORV/TGA", "DORV + TGA", "DORV & TGA"])
def test_compound_separators_all_split(raw):
    _hpo, _cats_, recognised, unrecognised = _parse_phenotype_input(raw)
    assert unrecognised == [], raw
    assert recognised, raw
    assert check_hpo_relevance(raw, "Congenital heart disease"), raw


def test_whole_chunk_wins_before_splitting():
    """D-TGA is an HPO synonym in its own right (dextrotransposition) and must
    NOT be shredded into "D" + "TGA"; same for hyphenated eponyms."""
    for raw in ("D-TGA", "L-TGA", "Wolff-Parkinson-White syndrome"):
        _h, _c, recognised, unrecognised = _parse_phenotype_input(raw)
        assert unrecognised == [], raw
        assert recognised == [raw], raw


def test_partial_compound_reports_only_the_unresolved_half():
    _h, _c, recognised, unrecognised = _parse_phenotype_input("TGA-zzzznotaphenotype")
    assert recognised, "the TGA half should still count"
    assert unrecognised == ["zzzznotaphenotype"]


def test_single_letter_fragments_are_not_treated_as_terms():
    _h, _c, recognised, unrecognised = _parse_phenotype_input("A-B")
    assert recognised == []
    assert unrecognised == ["A-B"]


def test_qualifier_prose_after_with_is_swallowed_not_flagged():
    """"with" joins a term to the curator's own description, so the tail must
    not be reported as a missing phenotype the way a "-" compound's half is."""
    _h, _c, recognised, unrecognised = _parse_phenotype_input(
        "HCM with apical involvement")
    assert unrecognised == []
    assert recognised == ["HCM with apical involvement"]
    assert check_hpo_relevance("HCM with apical involvement", "hcm")


ARRHYTHMIA = ["channelopathy", "conduction", "Arrhythmia / channelopathy"]
CHD = ["Congenital heart disease", "Other cardiac"]


def test_prolonged_qt_hpo_id_matches_channelopathy():
    """The inversion this fixes: the canonical LQTS ontology term used to
    match nothing while the free-text acronym matched."""
    assert check_hpo_relevance("HP:0001657", "channelopathy")
    assert check_hpo_relevance("HP:0005184", "channelopathy")
    assert check_hpo_relevance("HP:0012232", "channelopathy")


def test_free_text_and_ontology_id_now_agree():
    for cat in ARRHYTHMIA:
        assert check_hpo_relevance("LQTS", cat) == check_hpo_relevance("HP:0001657", cat), cat


def test_conduction_terms_match_conduction():
    assert check_hpo_relevance("HP:0001678", "conduction")
    assert check_hpo_relevance("HP:0001709", "conduction")


def test_valve_function_terms_match_congenital_heart_disease():
    for hid in ("HP:0001650", "HP:0001642", "HP:0001718", "HP:0001659"):
        assert check_hpo_relevance(hid, "Congenital heart disease"), hid
        assert check_hpo_relevance(hid, "Other cardiac"), hid


def test_pulmonary_arterial_hypertension_matches_its_own_bucket():
    assert check_hpo_relevance("HP:0002092", "Pulmonary hypertension")


def test_pulmonary_root_is_narrow_enough_to_exclude_embolism():
    """HP:0004890 (Elevated pulmonary artery pressure) is used rather than the
    broader HP:0030875 (Abnormality of pulmonary circulation) precisely so
    pulmonary/air embolism and intrapulmonary shunt stay out of the
    Pulmonary-hypertension panel bucket."""
    assert not check_hpo_relevance("HP:0002204", "Pulmonary hypertension")
    assert not check_hpo_relevance("HP:0033426", "Pulmonary hypertension")


def test_new_roots_do_not_broad_match_every_category():
    """CARDIAC_PARENT_TERMS is the "matches EVERY category" shortcut, reserved
    for the five top-level ancestors. The physiology roots are mid-level and
    must stay out of it, or submitting one would light up every panel."""
    for root in ("HP:0031653", "HP:0031546", "HP:0031547", "HP:0004890"):
        assert root in panelapp._PARENT_TO_CATEGORIES, root
        assert root not in panelapp.CARDIAC_PARENT_TERMS, root


def test_non_cardiac_terms_still_match_nothing():
    for hid in ("HP:0001250", "HP:0000407"):
        assert _cats(hid, ARRHYTHMIA + CHD + ["Cardiomyopathy"]) == [], hid


def test_hpo_exact_synonyms_are_recognised():
    """Free from hp.obo — no hand curation."""
    for term in ("ASD", "Aortic coarctation", "TAPVR", "PFO"):
        _h, _c, recognised, unrecognised = _parse_phenotype_input(term)
        assert unrecognised == [], term
        assert recognised == [term], term


def test_curated_acronyms_hpo_lacks():
    assert check_hpo_relevance("TOF", "Congenital heart disease")
    assert check_hpo_relevance("AS", "Congenital heart disease")
    assert check_hpo_relevance("AVB", "conduction")
    assert check_hpo_relevance("PAH", "Pulmonary hypertension")


def test_family_history_acronym_is_deliberately_not_mapped():
    """"FH" is far more often "family history" than familial
    hypercholesterolaemia in clinical shorthand. Mapping it would invent a
    phenotype the curator never claimed, so it stays unrecognised."""
    _h, _c, _r, unrecognised = _parse_phenotype_input("FH")
    assert unrecognised == ["FH"]


def test_short_lowercase_words_cannot_match_via_synonyms():
    """Synonym keys of <= 3 characters are admitted only from the curated map,
    never auto-derived from hp.obo, so ordinary English words can never
    become phenotype terms."""
    for word in ("is", "or", "at", "on", "no"):
        _h, _c, recognised, unrecognised = _parse_phenotype_input(word)
        assert recognised == [], word
        assert unrecognised == [word], word


def test_resolve_reports_the_interpretation_of_each_term():
    """Short acronyms are ambiguous (MS is mitral stenosis here, multiple
    sclerosis elsewhere) and phenotype gates PP4, so the reading has to be
    visible to the curator rather than silent."""
    out = resolve_phenotype_terms("AS, HP:0001644, zzzznotaphenotype")
    by_token = {t["token"]: t for t in out}
    assert by_token["AS"]["hpo_id"] == "HP:0001650"
    assert by_token["AS"]["label"] == "Aortic valve stenosis"
    assert by_token["AS"]["expanded"] is True
    assert by_token["HP:0001644"]["expanded"] is False
    assert by_token["zzzznotaphenotype"]["hpo_id"] is None


def test_compound_interpretation_lists_both_halves():
    out = resolve_phenotype_terms("DORV-TGA")
    assert len(out) == 1
    assert out[0]["token"] == "DORV-TGA"
    labels = " ".join(out[0]["parts"])
    assert "Double outlet right ventricle" in labels
    assert "ransposition" in labels


def test_existing_phrase_table_entries_all_still_resolve():
    for phrase in panelapp.PHENOTYPE_PHRASE_TO_CATEGORIES:
        _h, _c, recognised, unrecognised = _parse_phenotype_input(phrase)
        assert unrecognised == [], phrase
        assert recognised == [phrase], phrase


def test_parse_signature_and_empty_inputs_unchanged():
    assert _parse_phenotype_input(None) == ([], set(), [], [])
    assert _parse_phenotype_input("") == ([], set(), [], [])
    assert _parse_phenotype_input("   ,  ; ") == ([], set(), [], [])


def test_hpo_ids_are_still_uppercased_and_lists_accepted():
    hpo_ids, _c, _r, _u = _parse_phenotype_input(["hp:0001644", "DCM"])
    assert "HP:0001644" in hpo_ids
