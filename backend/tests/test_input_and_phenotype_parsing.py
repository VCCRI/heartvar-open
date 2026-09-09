"""Unit tests for the variant-input parser hardening and the ClinVar
phenotype-placeholder filter (deployment-risk audit Part 1B + Part 3 #1).

All offline / pure-function. Covers:
  - parse_variant_input: gene-symbol prefix strip ("MYH7:c…" → "c…"),
    unicode-dash + thousands-comma normalisation, and regression guards that
    transcript prefixes, bare HGVS, intronic hyphens, genomic accessions and
    coordinates are unchanged.
  - clinvar._split_phenotypes: the literal "N conditions" placeholder is
    dropped while real disease titles (incl. ones containing the word
    "condition") are kept.

No pytest dependency — runnable directly
(``python -m backend.tests.test_input_and_phenotype_parsing``).
"""
from __future__ import annotations

from backend.clients.ensembl_vep import parse_variant_input
from backend.clients.clinvar import _split_phenotypes


def test_gene_symbol_prefix_is_stripped():
    assert parse_variant_input("MYH7:c.1208G>A") == {
        "format": "hgvs", "hgvs": "c.1208G>A", "transcript": None,
    }
    assert parse_variant_input("NKX2-5:c.325C>T") == {
        "format": "hgvs", "hgvs": "c.325C>T", "transcript": None,
    }


def test_transcript_prefixes_still_split_not_treated_as_gene():
    assert parse_variant_input("NM_000257.4:c.1208G>A") == {
        "format": "hgvs", "hgvs": "c.1208G>A", "transcript": "NM_000257.4",
    }
    assert parse_variant_input("ENST00000355349:c.1208G>A") == {
        "format": "hgvs", "hgvs": "c.1208G>A", "transcript": "ENST00000355349",
    }


def test_bare_hgvs_unchanged():
    assert parse_variant_input("c.1208G>A") == {
        "format": "hgvs", "hgvs": "c.1208G>A", "transcript": None,
    }


def test_genomic_accession_not_stripped_as_gene():
    out = parse_variant_input("NC_000007.14:g.117548628G>A")
    assert out["format"] == "hgvs"
    assert "117548628G>A" in out["hgvs"]


def test_thousands_comma_coords():
    assert parse_variant_input("chr14:23,429,278:C:T") == {
        "format": "coordinates", "chrom": "14", "pos": 23429278,
        "ref": "C", "alt": "T", "build": None,
    }


def test_unicode_endash_coords():
    assert parse_variant_input("14–23429278–C–T") == {
        "format": "coordinates", "chrom": "14", "pos": 23429278,
        "ref": "C", "alt": "T", "build": None,
    }


def test_intronic_ascii_hyphen_preserved():
    out = parse_variant_input("c.88-2A>C")
    assert out == {"format": "hgvs", "hgvs": "c.88-2A>C", "transcript": None}


def test_plain_coords_and_unknown():
    assert parse_variant_input("7-117548628-C-T")["format"] == "coordinates"
    assert parse_variant_input("")["format"] == "unknown"
    assert parse_variant_input("just some prose")["format"] == "unknown"


def test_rsid_detected():
    assert parse_variant_input("rs727503113") == {"format": "rsid", "rsid": "rs727503113"}
    assert parse_variant_input("RS727503113")["format"] == "rsid"
    assert parse_variant_input("rs727503113")["rsid"] == "rs727503113"


def test_protein_only_flagged_not_unknown():
    assert parse_variant_input("p.Arg403Gln") == {
        "format": "protein", "hgvs_p": "p.Arg403Gln",
    }
    assert parse_variant_input("NP_000248.2:p.Arg403Gln")["format"] == "protein"


def test_coords_with_internal_whitespace():
    assert parse_variant_input("chr14 : 23429278 : C : T") == {
        "format": "coordinates", "chrom": "14", "pos": 23429278,
        "ref": "C", "alt": "T", "build": None,
    }


def test_n_conditions_placeholder_dropped():
    out = _split_phenotypes(
        "Hypertrophic cardiomyopathy 1|Cardiovascular phenotype|6 conditions|"
        "MYH7-related disorder|not provided"
    )
    assert "6 conditions" not in out
    assert "Hypertrophic cardiomyopathy 1" in out
    assert "Cardiovascular phenotype" in out
    assert "MYH7-related disorder" in out
    assert "not provided" not in out


def test_real_condition_titles_retained():
    out = _split_phenotypes("LMNA-associated condition|unspecified heart condition")
    assert out == ["LMNA-associated condition", "unspecified heart condition"]


def test_singular_and_caps_placeholder_dropped():
    assert _split_phenotypes("1 condition|30 CONDITIONS|Long QT syndrome") == [
        "Long QT syndrome"
    ]


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
