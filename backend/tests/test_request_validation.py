"""Unit tests for the CurationRequest gene-symbol validator (deployment-risk
audit Part 2 #8 — input → URL-path injection surface).

`gene` is interpolated into REST URL paths; the validator rejects anything that
isn't a plausible HGNC symbol while accepting every real cardiac-gene form (incl.
hyphen/dot/digit symbols and the empty-string coordinate-input case).

No pytest dependency — runnable directly
(``python -m backend.tests.test_request_validation``).
"""
from __future__ import annotations

from pydantic import ValidationError

from backend.app import CurationRequest


def test_real_gene_symbols_accepted():
    for g in ["MYH7", "KCNQ1", "NKX2-5", "MT-TL1", "C1orf127", "RAF1", "PKP2"]:
        req = CurationRequest(gene=g, hgvs_c="c.1A>G")
        assert req.gene == g


def test_empty_gene_accepted_for_coordinate_input():
    req = CurationRequest(gene="", hgvs_c="14-23429278-C-T")
    assert req.gene == ""


def test_gene_is_stripped():
    assert CurationRequest(gene="  MYH7  ", hgvs_c="c.1A>G").gene == "MYH7"


def test_injection_like_genes_rejected():
    for bad in ["MYH7; DROP", "MYH7/../etc", "MYH7 OR 1=1", "../secret", "a b", "x%2f"]:
        try:
            CurationRequest(gene=bad, hgvs_c="c.1A>G")
        except ValidationError:
            continue
        raise AssertionError(f"expected ValidationError for gene={bad!r}")


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
