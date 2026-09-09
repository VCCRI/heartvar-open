"""The ClinGen SVI gene-disease-validity ceiling.

The rule, in the RASopathy VCEP's own words on MRAS NM_001085049.3:c.67G>C:

    "Given the Moderate strength of gene-disease relationship between MRAS and
     autosomal dominant RASopathy, ClinGen's sequence variant interpretation
     working group does not recommend the classification of variants in this
     gene beyond likely pathogenic."

HeartVar has shipped the ClinGen Gene-Disease Validity table since 2026-08 and
reads its `classification` field, but ONLY as an adequacy gate for
mode-of-inheritance (_ADEQUATE_TIERS in gene_inheritance_modes). The tier was
never applied as a ceiling, so an MRAS variant could score Pathogenic where the
VCEP says the gene does not support a call beyond Likely pathogenic.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from backend.acmg.hard_coded import gene_validity_ceiling
from backend.acmg.tiers import (
    apply_gene_validity_ceiling,
    classification_for_criteria,
    compute_points_total,
)


def _c(code, strength=None):
    return {"code": code, "criteria_strength": strength or code, "status": "met"}


def test_mras_is_moderate_so_it_gets_a_ceiling():
    assert gene_validity_ceiling("MRAS") == "Moderate"


def test_a_gene_definitive_for_anything_gets_no_ceiling():
    assert gene_validity_ceiling("MYH7") == "Definitive"
    assert gene_validity_ceiling("PTPN11") == "Definitive"
    assert gene_validity_ceiling("LZTR1") == "Definitive"


def test_unknown_gene_fails_open():
    for bad in ("", "NOTAGENE", None, "  "):
        assert gene_validity_ceiling(bad) is None


def test_lookup_is_case_insensitive():
    assert gene_validity_ceiling("mras") == "Moderate"


def test_moderate_validity_caps_pathogenic_at_likely_pathogenic():
    assert apply_gene_validity_ceiling("Pathogenic", "Moderate") == "Likely pathogenic"


def test_weaker_than_moderate_also_caps():
    for tier in ("Limited", "Disputed", "Refuted", "No Known Disease Relationship"):
        assert apply_gene_validity_ceiling("Pathogenic", tier) == "Likely pathogenic"


def test_definitive_and_strong_do_not_cap():
    for tier in ("Definitive", "Strong"):
        assert apply_gene_validity_ceiling("Pathogenic", tier) == "Pathogenic"


def test_the_ceiling_never_touches_any_other_tier():
    for tier in ("Likely pathogenic", "VUS", "Likely benign", "Benign"):
        assert apply_gene_validity_ceiling(tier, "Moderate") == tier


def test_absent_validity_fails_open():
    assert apply_gene_validity_ceiling("Pathogenic", None) == "Pathogenic"
    assert apply_gene_validity_ceiling("Pathogenic", "") == "Pathogenic"


def test_mras_c67G_C_lands_on_the_panels_tier():
    """MRAS c.67G>C: HeartVar scored +12 = Pathogenic, panel said LP."""
    criteria = [_c("PS2", "PS2_Strong"), _c("PS4"), _c("PM1"),
                _c("PM2", "PM2_Supporting"), _c("PP3", "PP3_Supporting")]
    pts = compute_points_total(criteria)
    assert pts >= 10, pts
    assert classification_for_criteria(pts, criteria) == "Pathogenic"
    assert classification_for_criteria(
        pts, criteria, validity="Moderate") == "Likely pathogenic"


def test_a_definitive_gene_with_the_same_points_is_untouched():
    criteria = [_c("PS2", "PS2_Strong"), _c("PS4"), _c("PM1"),
                _c("PM2", "PM2_Supporting"), _c("PP3", "PP3_Supporting")]
    pts = compute_points_total(criteria)
    assert classification_for_criteria(
        pts, criteria, validity="Definitive") == "Pathogenic"


def test_forced_classification_still_wins():
    assert classification_for_criteria(
        12, [_c("PS4")], "Benign", validity="Moderate") == "Benign"


def test_ceiling_ignores_genes_with_no_vcep_spec():
    """ATP2A2 holds exactly one ClinGen assertion — Refuted, for EPILEPSY. A
    refuted epilepsy link must not hold an ATP2A2 cardiac call at Likely
    pathogenic. Same shape: PIK3CA (Refuted, breast carcinoma), COL11A1
    (Moderate, hearing loss), CACNA1S (Moderate, malignant hyperthermia)."""
    from backend.acmg.hard_coded import (
        gene_validity_ceiling, gene_validity_caps_pathogenic,
    )
    for gene in ("ATP2A2", "PIK3CA", "COL11A1", "CACNA1S"):
        assert gene_validity_ceiling(gene) is None, gene
        assert gene_validity_caps_pathogenic(gene) is False, gene


def test_mras_is_still_capped():
    """The rule's motivating case must survive the scoping: MRAS has a
    RASopathy VCEP spec and a Moderate ClinGen assertion."""
    from backend.acmg.hard_coded import (
        gene_validity_ceiling, gene_validity_caps_pathogenic,
    )
    assert gene_validity_ceiling("MRAS") == "Moderate"
    assert gene_validity_caps_pathogenic("MRAS") is True


def test_exactly_one_gene_in_the_table_is_capped():
    """A count pin. If a future harvest widens the spec'd gene set, this fails
    and the new caps get looked at deliberately rather than shipping silently."""
    import json
    from pathlib import Path
    from backend.acmg.hard_coded import gene_validity_caps_pathogenic
    path = (Path(__file__).resolve().parents[1]
            / "data" / "clingen_gene_validity.json")
    genes = json.loads(path.read_text())["genes"]
    capped = sorted(g for g in genes if gene_validity_caps_pathogenic(g))
    assert capped == ["MRAS"], capped


def test_specd_genes_with_a_strong_relationship_are_never_capped():
    from backend.acmg.hard_coded import gene_validity_caps_pathogenic
    for gene in ("MYH7", "TNNI3", "FBN1", "KCNQ1", "MYBPC3", "PTPN11"):
        assert gene_validity_caps_pathogenic(gene) is False, gene
