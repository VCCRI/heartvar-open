"""Unit tests for the PubTator3 PS4-retrieval client's pure helpers.

Covers the precision-critical bits offline (no network): one-letter AA
derivation, the gene-match guard that prevents resolving to a wrong-gene
entity (the ACTC1 I284F → SARS-CoV-2 ORF1ab misresolution), candidate-term
ordering, and the BioC document iterator.

No pytest dependency — runnable directly
(``python -m backend.tests.test_pubtator3``).
"""
from __future__ import annotations

from backend.clients.pubtator3 import (
    _one_letter_change,
    _gene_matches,
    _candidate_terms,
    _iter_documents,
    _bare_c,
)


def test_one_letter_change_missense():
    assert _one_letter_change("Gly265Arg") == "G265R"
    assert _one_letter_change("Arg495Trp") == "R495W"
    assert _one_letter_change("Ile284Phe") == "I284F"


def test_one_letter_change_rejects_non_missense():
    assert _one_letter_change("Arg495Ter") is None
    assert _one_letter_change("") is None
    assert _one_letter_change("c.26-2A>G") is None
    assert _one_letter_change("Arg495fs") is None


def test_gene_matches_rejects_wrong_gene():
    covid = {"_id": "@VARIANT_p.I284F_ORF1ab_Severe_acute_respiratory_syndrome_coronavirus_2",
             "description": "ORF1ab"}
    assert _gene_matches("ACTC1", covid) is False
    actc1 = {"_id": "@VARIANT_p.I284F_ACTC1_human", "description": "ACTC1 (human)"}
    assert _gene_matches("ACTC1", actc1) is True


def test_gene_matches_word_boundary_not_substring():
    braf = {"_id": "@VARIANT_p.G265R_BRAF_human", "description": "BRAF (human)"}
    assert _gene_matches("BRAF", braf) is True
    assert _gene_matches("RAF1", braf) is False


def test_candidate_terms_are_gene_qualified_no_bare_token():
    terms = _candidate_terms("ACTC1", "NM_005159.5:c.850A>T", "Ile284Phe", None)
    assert terms[0] == "ACTC1 I284F"
    assert "I284F" not in terms
    assert "ACTC1 Ile284Phe" in terms and "ACTC1 c.850A>T" in terms


def test_candidate_terms_rsid_first_when_present():
    terms = _candidate_terms("MYBPC3", "NM_000256.3:c.1483C>T", "Arg495Trp", "rs397515905")
    assert terms[0] == "rs397515905"


def test_bare_c_strips_transcript():
    assert _bare_c("NM_000256.3:c.1483C>T") == "c.1483C>T"
    assert _bare_c("c.850A>T") == "c.850A>T"


def test_iter_documents_handles_shapes():
    doc = {"passages": [{"text": "x"}]}
    assert list(_iter_documents({"PubTator3": [doc]})) == [doc]
    assert list(_iter_documents({"documents": [doc]})) == [doc]
    assert list(_iter_documents([doc])) == [doc]
    assert list(_iter_documents(doc)) == [doc]
    assert list(_iter_documents({})) == []


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
