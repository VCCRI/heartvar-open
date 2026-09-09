"""Unit tests for backend.acmg.hard_coded.gene_mechanism.

The deterministic, source-attributed disease-mechanism classifier merges the
curated signals HeartVar holds — ClinGen Gene-Dosage scores plus the live VCEP
PVS1/PP2 applicability, the CHDgene dominant-negative flag, the hand GoF/DN
list, and gnomAD constraint — into a single LoF / GoF / DN / mixed / undetermined
call. It must never invent a mechanism (no curated signal → "undetermined") and
must report genes with evidence on BOTH axes as "mixed".

Informational only: this never changes an ACMG criterion or score.

No pytest dependency — runnable with pytest or directly.
"""
from __future__ import annotations

from backend.app import gene_mechanism


def _ev(mis_z: float | None = None, loeuf: float | None = None) -> dict:
    constraint = {}
    if mis_z is not None:
        constraint["mis_z"] = mis_z
    if loeuf is not None:
        constraint["oe_lof_upper"] = loeuf
    return {"gnomad": {"gene": {"gnomad_constraint": constraint}}} if constraint else {}


def test_haploinsufficiency_genes():
    for gene in ("GATA4", "NKX2-5", "TBX5", "MYBPC3"):
        m = gene_mechanism({}, gene)
        assert m["mechanism"] == "haploinsufficiency", (gene, m)
        assert m["confidence"] == "established"
        assert any(s["name"] == "ClinGen Dosage" for s in m["sources"])


def test_recessive_lof_gene():
    m = gene_mechanism({}, "ACADVL")
    assert m["mechanism"] == "recessive_lof"
    assert m["confidence"] == "established"


def test_gof_or_dn_genes():
    for gene in ("HRAS", "BRAF", "MYH7"):
        m = gene_mechanism({}, gene)
        assert m["mechanism"] in ("gof_or_dn", "dominant_negative"), (gene, m)
        # Attributed to the expert/curated source, not invented.
        names = {s["name"] for s in m["sources"]}
        assert names & {"RASopathy VCEP", "Cardiomyopathy VCEP", "HeartVar curated"}


def test_mixed_mechanism_genes():
    for gene in ("SCN5A", "PTPN11"):
        m = gene_mechanism({}, gene)
        assert m["mechanism"] == "mixed", (gene, m)
        assert "both" in m["summary"].lower()


def test_undetermined_when_no_curated_signal():
    m = gene_mechanism({}, "RYR2")
    assert m["mechanism"] == "undetermined"
    assert m["confidence"] == "undetermined"


def test_constraint_is_only_a_hint_never_the_mechanism():
    m = gene_mechanism(_ev(loeuf=0.10), "RYR2")
    assert m["mechanism"] == "undetermined"
    assert m["constraint_hint"] is not None
    assert "proxy" in m["constraint_hint"].lower()


def test_emerging_haploinsufficiency():
    m = gene_mechanism({}, "TTN")
    assert m["mechanism"] == "haploinsufficiency"
    assert m["confidence"] == "emerging"


def test_empty_gene_is_undetermined():
    m = gene_mechanism({}, "")
    assert m["mechanism"] == "undetermined"
    assert m["sources"] == []


def test_return_shape_is_stable():
    m = gene_mechanism({}, "GATA4")
    for key in ("mechanism", "label", "confidence", "sources", "hi_score", "summary"):
        assert key in m, key
    assert isinstance(m["sources"], list)
    assert m["label"] in (
        "Loss-of-function (haploinsufficiency)",
        "Loss-of-function (recessive / biallelic)",
        "Dominant-negative",
        "Gain-of-function or dominant-negative",
        "Mixed mechanism (LoF and GoF/DN both reported)",
        "Mechanism not established",
    )


def test_case_insensitive_gene_lookup():
    assert gene_mechanism({}, "gata4")["mechanism"] == "haploinsufficiency"


if __name__ == "__main__":  # pragma: no cover
    import sys

    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"ok  {fn.__name__}")
    print(f"\n{len(fns)} tests passed")
    sys.exit(0)
