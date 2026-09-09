"""Unit tests for the per-term phenotype-match audit that backs the
"these terms match nothing for this gene" line at the bottom of the
Gene-disease validity card.

Collaborator request: the validity dropdown should say so when the
curator submitted a phenotype term that matches NONE of the gene's
curated disease records — today a non-matching term is silently absent
from the per-row ✓ ticks, which reads the same as "no phenotype given".

All offline / pure-function: `audit_phenotype_terms` takes the raw
phenotype string plus the set of disease categories the gene actually
has, and reports each term's match state. Category resolution reuses
`check_hpo_relevance`, the same predicate that drives the per-row ticks,
so this line can never contradict a ✓ shown above it.
"""
from __future__ import annotations

import pytest

from backend.clients import panelapp
from backend.clients.panelapp import audit_phenotype_terms


def _load_descendants() -> None:
    """Load the offline HPO descendant closure the app ships.

    Without it `check_hpo_relevance` silently drops to the 61-ID curated seed
    fallback, which is NOT production behaviour — and that fallback lists
    HP:0001250 (Seizure) under 'Cardiomyopathy', so a test written against it
    would assert the opposite of the truth.
    """
    panelapp._load_local_descendants()


@pytest.fixture(autouse=True)
def _descendants():
    _load_descendants()


def test_unmatched_term_is_reported_with_its_label():
    out = audit_phenotype_terms("HP:0001250", {"Cardiomyopathy", "dcm"})
    assert [t["token"] for t in out["unmatched"]] == ["HP:0001250"]
    assert out["unmatched"][0]["label"] == "Seizure"
    assert out["matched"] == []


def test_matching_term_is_not_reported():
    out = audit_phenotype_terms("HP:0001644", {"Cardiomyopathy", "dcm"})
    assert out["unmatched"] == []
    assert [t["token"] for t in out["matched"]] == ["HP:0001644"]


def test_mixed_input_reports_only_the_unmatched_terms():
    out = audit_phenotype_terms(
        "HP:0001644, HP:0001250", {"Cardiomyopathy", "dcm"}
    )
    assert [t["token"] for t in out["matched"]] == ["HP:0001644"]
    assert [t["token"] for t in out["unmatched"]] == ["HP:0001250"]


def test_free_text_phrase_matches_via_the_phrase_table():
    out = audit_phenotype_terms("dilated cardiomyopathy", {"Cardiomyopathy"})
    assert out["unmatched"] == []
    assert [t["token"] for t in out["matched"]] == ["dilated cardiomyopathy"]


def test_free_text_phrase_that_matches_nothing_is_labelled_with_its_ontology_term():
    out = audit_phenotype_terms("long QT syndrome", {"Cardiomyopathy"})
    assert [t["token"] for t in out["unmatched"]] == ["long QT syndrome"]
    assert out["unmatched"][0]["label"] == "Prolonged QT interval"


def test_free_text_with_no_ontology_term_falls_back_to_the_curators_wording():
    out = audit_phenotype_terms("hypertrophic cardiomyopathy 1", {"Congenital heart disease"})
    assert out["unmatched"][0]["label"] == "hypertrophic cardiomyopathy 1"


def test_no_phenotype_submitted_reports_nothing():
    for raw in ("", "   ", None):
        out = audit_phenotype_terms(raw, {"Cardiomyopathy"})
        assert out["unmatched"] == []
        assert out["matched"] == []


def test_gene_with_no_known_categories_is_suppressed():
    out = audit_phenotype_terms("HP:0001250, HP:0001644", set())
    assert out["unmatched"] == []
    assert out["matched"] == []
    assert out["suppressed"] is True


def test_unrecognised_tokens_are_not_reported_as_unmatched():
    out = audit_phenotype_terms("zzzznotaphenotype", {"Cardiomyopathy"})
    assert out["unmatched"] == []
    assert out["unrecognised"] == ["zzzznotaphenotype"]


def test_malformed_hpo_id_is_unrecognised_not_unmatched():
    out = audit_phenotype_terms("HP:123", {"Cardiomyopathy"})
    assert out["unmatched"] == []
    assert out["unrecognised"] == ["HP:123"]


def test_chd_term_matches_when_the_gene_is_a_chdgene_entry():
    out = audit_phenotype_terms("HP:0001631", {"Congenital heart disease"})
    assert out["unmatched"] == []
    assert [t["token"] for t in out["matched"]] == ["HP:0001631"]


def test_chd_phrase_matches_the_congenital_heart_category():
    out = audit_phenotype_terms(
        "atrial septal defect", {"Congenital heart disease"}
    )
    assert out["unmatched"] == []


def test_non_chd_term_still_unmatched_on_a_chd_only_gene():
    out = audit_phenotype_terms("HP:0001663", {"Congenital heart disease"})
    assert [t["token"] for t in out["unmatched"]] == ["HP:0001663"]


def test_audit_agrees_with_check_hpo_relevance_per_category():
    from backend.clients.panelapp import check_hpo_relevance

    cats = {"Cardiomyopathy", "Congenital heart disease"}
    for term in ("HP:0001644", "HP:0001631", "HP:0001250", "HP:0001279"):
        out = audit_phenotype_terms(term, cats)
        expected = any(check_hpo_relevance(term, c) for c in cats)
        assert bool(out["matched"]) is expected, term
        assert bool(out["unmatched"]) is (not expected), term


def test_duplicate_terms_are_reported_once():
    out = audit_phenotype_terms(
        "HP:0001250, HP:0001250", {"Cardiomyopathy"}
    )
    assert [t["token"] for t in out["unmatched"]] == ["HP:0001250"]


def test_categories_union_across_gencc_clingen_and_panelapp():
    from backend.evidence import _gene_disease_categories

    ev = {
        "gencc": {
            "submissions": [{"category": "dcm"}, {"category": None}],
            "clingen_submissions": [{"category": "Cardiomyopathy"}],
        },
        "panelapp": {"panels_found": [{"category": "Arrhythmia / channelopathy"}]},
    }
    assert _gene_disease_categories(ev) == {
        "dcm", "Cardiomyopathy", "Arrhythmia / channelopathy",
    }


def test_chdgene_listing_contributes_the_congenital_bucket():
    from backend.evidence import _gene_disease_categories

    ev = {"chdgene": {"ok": True, "listed": True}}
    assert _gene_disease_categories(ev) == {"Congenital heart disease"}


def test_chdgene_not_listed_or_failed_contributes_nothing():
    from backend.evidence import _gene_disease_categories

    assert _gene_disease_categories({"chdgene": {"ok": True, "listed": False}}) == set()
    assert _gene_disease_categories({"chdgene": {"ok": False, "listed": True}}) == set()


def test_empty_evidence_yields_no_categories():
    from backend.evidence import _gene_disease_categories

    assert _gene_disease_categories({}) == set()


def test_chd_proband_on_a_chdgene_only_gene_is_not_flagged():
    """End-to-end of the user's call: NKX2-5-style gene whose only record is a
    CHDgene listing, proband submitted an ASD term → no warning line."""
    from backend.evidence import _gene_disease_categories

    ev = {"chdgene": {"ok": True, "listed": True}}
    out = audit_phenotype_terms("HP:0001631", _gene_disease_categories(ev))
    assert out["unmatched"] == []
    assert out["suppressed"] is True


if __name__ == "__main__":  # pragma: no cover - manual run helper
    import sys

    _load_descendants()
    mod = sys.modules[__name__]
    failures = 0
    for name in sorted(n for n in dir(mod) if n.startswith("test_")):
        try:
            getattr(mod, name)()
            print(f"PASS {name}")
        except AssertionError as exc:
            failures += 1
            print(f"FAIL {name}: {exc}")
    sys.exit(1 if failures else 0)
