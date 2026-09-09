"""Unit tests for backend.app._gene_lof_mechanism.

Guards the mechanism-exclusion gate: a gain-of-function / dominant-negative
gene must not be assigned a LoF disease mechanism for a non-canonical-LoF
(e.g. missense) variant, even when it is strongly LoF-constrained in gnomAD.
This is the root cause of the spurious PVS1 over-calls on BRAF c.1024A>G and
HRAS c.508A>T (eRepo benchmark, severity-4 over-calls #6/#7).

No pytest dependency — runnable with pytest or directly
(``python -m backend.tests.test_gene_lof_mechanism``).
"""
from __future__ import annotations

from backend.app import _gene_lof_mechanism

_CONSTRAINED = {
    "gnomad": {"gene": {"gnomad_constraint": {"pLI": 0.99, "oe_lof_upper": 0.2}}},
}
_CHD_AD = {"chdgene": {"listed": True, "inheritance": ["AD"]}}

_GENCC_AR_DEFINITIVE = {"gencc": {"submissions": [
    {"classification": "Moderate", "moi": "Autosomal recessive"},
    {"classification": "Definitive", "moi": "Autosomal recessive"},
    {"classification": "Limited", "moi": "Autosomal dominant"},
]}}
_CHD_AR = {"chdgene": {"listed": True, "inheritance": ["AR", "CH"]}}
_GENCC_XL_STRONG = {"gencc": {"submissions": [
    {"classification": "Strong", "moi": "X-linked recessive"},
]}}
_GENCC_AR_MODERATE = {"gencc": {"submissions": [
    {"classification": "Moderate", "moi": "Autosomal recessive"},
]}}


def test_gof_gene_missense_constraint_suppressed():
    """BRAF/HRAS missense: LoF-constrained but GoF mechanism → False."""
    assert _gene_lof_mechanism(_CONSTRAINED, "BRAF", "missense_variant") is False
    assert _gene_lof_mechanism(_CONSTRAINED, "HRAS", "missense_variant") is False


def test_gof_gene_canonical_lof_constraint_applies():
    """A genuine LoF event (frameshift/stop) in a dominant-negative-list gene
    that the CSpec does NOT exclude from PVS1 → constraint tier re-enabled, so
    True (PVS1 routing stays correct). SCN5A is the canonical example named in
    the _gene_lof_mechanism docstring; it is in _DOMINANT_NEGATIVE_GENES but its
    VCEP keeps PVS1 applicable, so the consequence-aware carve-out governs it.
    (BRAF cannot be used here — its RASopathy CSpec marks PVS1 Not Applicable, so
    the hard override below returns False for it on every consequence.)"""
    assert _gene_lof_mechanism(_CONSTRAINED, "SCN5A", "frameshift_variant") is True
    assert _gene_lof_mechanism(_CONSTRAINED, "SCN5A", "stop_gained") is True


def test_cspec_pvs1_not_applicable_overrides_canonical_lof():
    """CSpec HARD override: a gene whose VCEP marks PVS1 Not Applicable (e.g.
    BRAF/HRAS, RASopathy GoF mechanism) gets no LoF mechanism for ANY
    consequence — even a frameshift/stop/splice canonical-LoF — because the
    VCEP determined LoF is not the disease mechanism. This is stronger than the
    dominant-negative hand-list carve-out and prevents the spurious PVS1 that
    over-called HRAS/BRAF truncating variants to LP."""
    assert _gene_lof_mechanism(_CONSTRAINED, "BRAF", "frameshift_variant") is False
    assert _gene_lof_mechanism(_CHD_AD, "BRAF", "splice_donor_variant") is False
    assert _gene_lof_mechanism(_CONSTRAINED, "HRAS", "stop_gained") is False


def test_gof_gene_curator_tier_suppressed_for_missense():
    """Curator (CHDgene AD) tier is gated out for a dominant-negative-list gene's
    missense, but re-enabled for a canonical LoF consequence. Uses SCN5A (DN-list
    but CSpec-PVS1-applicable) so the carve-out — not the CSpec override — is
    what's exercised."""
    assert _gene_lof_mechanism(_CHD_AD, "SCN5A", "missense_variant") is False
    assert _gene_lof_mechanism(_CHD_AD, "SCN5A", "splice_donor_variant") is True


def test_non_gof_gene_constraint_unaffected():
    """A gene outside the exclusion list still gets LoF from constraint
    on any consequence (e.g. KCNQ1 haploinsufficiency)."""
    assert _gene_lof_mechanism(_CONSTRAINED, "KCNQ1", "missense_variant") is True


def test_no_evidence_returns_false():
    assert _gene_lof_mechanism({}, "BRAF", "missense_variant") is False
    assert _gene_lof_mechanism({}, "KCNQ1", "missense_variant") is False


def test_recessive_gencc_definitive_canonical_lof_admitted():
    """The DNAH9 fix: a GenCC-Definitive autosomal-recessive gene gets a LoF
    mechanism for a canonical null variant. Per ClinGen SVI, PVS1 fires on the
    null allele regardless of MoI; zygosity (PM3) governs the genotype call."""
    assert _gene_lof_mechanism(_GENCC_AR_DEFINITIVE, "DNAH9", "frameshift_variant") is True
    assert _gene_lof_mechanism(_GENCC_AR_DEFINITIVE, "DNAH9", "stop_gained") is True
    assert _gene_lof_mechanism(_GENCC_AR_DEFINITIVE, "DNAH9", "splice_donor_variant") is True


def test_chdgene_recessive_inheritance_admitted():
    """CHDgene tier now admits AR/CH (not only AD) — DNAH9 is a 5-star CHDgene
    Heterotaxy gene with inheritance ['AR','CH']."""
    assert _gene_lof_mechanism(_CHD_AR, "DNAH9", "frameshift_variant") is True


def test_xlinked_gencc_strong_admitted():
    """X-linked Definitive/Strong gene also admitted for a canonical null."""
    assert _gene_lof_mechanism(_GENCC_XL_STRONG, "SOMEGENE", "stop_gained") is True


def test_recessive_weak_gencc_floor_preserved():
    """The Definitive/Strong floor still holds: a recessive gene with only
    Moderate/Limited GenCC and no constraint / CHDgene listing → False."""
    assert _gene_lof_mechanism(_GENCC_AR_MODERATE, "SOMEGENE", "frameshift_variant") is False


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
