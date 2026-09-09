"""Guards on the local gene → MANE Select lookup.

WHY IT EXISTS. A curator typing a gene and a bare ``c.`` — which is all the web
form asks for — had no transcript to anchor to, so hgvs_resolver returned None
and the curation went to ``rest.ensembl.org``. On 2026-08-31 Ensembl answered
with ReadTimeout, then HTTP 500 four times, then 503, and one CHD7 c.5058del
cost 89 s. Reading MANE Select out of the cdot set we already ship removes that
host from the path.

THE RISK THIS CARRIES, and what these tests hold the line on: choosing a
transcript chooses the coordinates the whole ACMG evaluation runs on. So the
lookup must return the DESIGNATED transcript or nothing at all — never a "first
row", "longest CDS" or "newest version" consolation prize. Every test below is
really about that one property.

Stdlib only: the fixture builds a cdot-shaped SQLite by hand, so this runs
without the 608 MB real database.
"""
from __future__ import annotations

import gzip
import os
import json
import sqlite3

import pytest

from backend.clients import cdot_sqlite

REFSEQ_TAG = "MANE Select"
ENSEMBL_TAG = "CCDS,gencode_basic,gencode_primary,MANE_Select,Ensembl_canonical"


def _record(gene: str, tag: str | None, assembly: str = "GRCh38") -> bytes:
    build: dict = {"contig": "NC_000008.11", "strand": "+", "exons": [[1, 2, 0, 1, 100]]}
    if tag is not None:
        build["tag"] = tag
    return gzip.compress(json.dumps({
        "gene_name": gene,
        "biotype": ["mRNA"],
        "genome_builds": {assembly: build},
    }).encode())


def _db(tmp_path, rows: list[tuple[str, str, str | None]], name="cdot.db") -> str:
    """rows = [(accession, gene, tag_or_None)] → a cdot-shaped database path."""
    path = tmp_path / name
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE transcripts (accession TEXT PRIMARY KEY, "
                "gene TEXT, data BLOB NOT NULL)")
    con.executemany(
        "INSERT INTO transcripts (accession, gene, data) VALUES (?, ?, ?)",
        [(acc, gene, _record(gene, tag)) for acc, gene, tag in rows],
    )
    con.commit()
    con.close()
    cdot_sqlite.mane_select_for_gene.cache_clear()
    return str(path)


@pytest.mark.parametrize("tag", [REFSEQ_TAG, ENSEMBL_TAG, "mane_select",
                                 "MANE_SELECT", "a,MANE Select,b"])
def test_both_namespace_spellings_are_recognised(tag):
    assert cdot_sqlite._is_mane_select(tag) is True


@pytest.mark.parametrize("tag", [
    None, "", "CCDS,gencode_basic", "Ensembl_canonical", "basic",
    "MANE Plus Clinical",
    "CCDS,gencode_basic,MANE_Plus_Clinical",
])
def test_non_mane_select_tags_are_rejected(tag):
    assert cdot_sqlite._is_mane_select(tag) is False


def test_returns_the_mane_select_in_both_namespaces(tmp_path):
    db = _db(tmp_path, [
        ("NM_017780.4", "CHD7", REFSEQ_TAG),
        ("ENST00000423902.7", "CHD7", ENSEMBL_TAG),
        ("NM_017780.3", "CHD7", None),
        ("ENST00000524564.5", "CHD7", "CCDS,gencode_basic"),
    ])
    mane = cdot_sqlite.mane_select_for_gene(db, "CHD7")
    assert mane.ensembl == "ENST00000423902.7"
    assert mane.refseq == "NM_017780.4"


def test_a_gene_with_no_mane_select_returns_nothing(tmp_path):
    """The REST fallback has to survive. Returning the only transcript present
    would be the guess this lookup exists to avoid."""
    db = _db(tmp_path, [
        ("NM_999999.1", "OBSCURE1", None),
        ("ENST00000000001.1", "OBSCURE1", "gencode_basic"),
    ])
    assert cdot_sqlite.mane_select_for_gene(db, "OBSCURE1") == (None, None)


def test_other_genes_are_not_consulted(tmp_path):
    """The query is keyed on gene. A MANE Select on a DIFFERENT gene must not
    leak in — that would annotate the variant against the wrong locus."""
    db = _db(tmp_path, [
        ("ENST00000423902.7", "CHD7", ENSEMBL_TAG),
        ("NM_999999.1", "OBSCURE1", None),
    ])
    assert cdot_sqlite.mane_select_for_gene(db, "OBSCURE1") == (None, None)
    assert cdot_sqlite.mane_select_for_gene(db, "CHD7").ensembl == "ENST00000423902.7"


def test_the_tag_is_read_from_the_requested_assembly(tmp_path):
    """A GRCh37-only MANE tag must not answer a GRCh38 question."""
    db = _db(tmp_path, [("ENST00000423902.7", "CHD7", ENSEMBL_TAG)])
    assert cdot_sqlite.mane_select_for_gene(db, "CHD7", "GRCh37") == (None, None)
    assert cdot_sqlite.mane_select_for_gene(db, "CHD7", "GRCh38").ensembl


@pytest.mark.parametrize("db_path,gene", [
    (None, "CHD7"), ("", "CHD7"), ("/nonexistent/cdot.db", "CHD7"),
])
def test_a_missing_database_is_not_an_error(db_path, gene):
    """Contract shared with the whole offline path: serve nothing, let the
    caller fall back. Never raise into a curation."""
    cdot_sqlite.mane_select_for_gene.cache_clear()
    assert cdot_sqlite.mane_select_for_gene(db_path, gene) == (None, None)


def test_no_gene_is_not_an_error(tmp_path):
    db = _db(tmp_path, [("ENST00000423902.7", "CHD7", ENSEMBL_TAG)])
    assert cdot_sqlite.mane_select_for_gene(db, None) == (None, None)
    assert cdot_sqlite.mane_select_for_gene(db, "") == (None, None)


def test_a_corrupt_record_is_skipped_not_raised(tmp_path):
    """One unreadable blob must not cost the gene its MANE transcript."""
    path = tmp_path / "corrupt.db"
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE transcripts (accession TEXT PRIMARY KEY, "
                "gene TEXT, data BLOB NOT NULL)")
    con.execute("INSERT INTO transcripts VALUES (?, ?, ?)",
                ("BROKEN.1", "CHD7", b"this is not gzipped json"))
    con.execute("INSERT INTO transcripts VALUES (?, ?, ?)",
                ("ENST00000423902.7", "CHD7", _record("CHD7", ENSEMBL_TAG)))
    con.commit()
    con.close()
    cdot_sqlite.mane_select_for_gene.cache_clear()
    assert cdot_sqlite.mane_select_for_gene(str(path), "CHD7").ensembl == \
        "ENST00000423902.7"


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))


def test_the_provider_refuses_to_fetch_sequence_over_the_network():
    """A provider whose seqfetcher can reach the network is the bug."""
    pytest.importorskip("cdot", reason="needs the real cdot/hgvs stack")
    pytest.importorskip("hgvs")
    from hgvs.exceptions import HGVSDataNotAvailableError

    provider = cdot_sqlite.make_provider("/nonexistent/cdot.db")
    assert provider is None, "a missing database must not yield a provider"

    import sqlite3 as _sq, tempfile, pathlib as _pl
    tmp = _pl.Path(tempfile.mkdtemp()) / "cdot.db"
    con = _sq.connect(tmp)
    con.execute("CREATE TABLE transcripts (accession TEXT PRIMARY KEY, "
                "gene TEXT, data BLOB NOT NULL)")
    con.execute("CREATE TABLE genes (key TEXT PRIMARY KEY, data BLOB NOT NULL)")
    con.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)")
    con.execute("INSERT INTO transcripts VALUES (?,?,?)",
                ("ENST00000355349.4", "MYH7", _record("MYH7", ENSEMBL_TAG)))
    con.commit(); con.close()

    provider = cdot_sqlite.make_provider(str(tmp))
    if provider is None:
        pytest.skip("provider could not be built in this environment")
    fetcher = provider.seqfetcher
    assert "offline" in getattr(fetcher, "source", "").lower(), (
        f"the provider kept a network seqfetcher ({fetcher!r}) — del/dup/inv "
        "will silently fetch a reference base from rest.ensembl.org"
    )
    with pytest.raises(HGVSDataNotAvailableError):
        fetcher.fetch_seq("ENST00000355349.4", 0, 10)


@pytest.mark.skipif(
    not (os.environ.get("HEARTVAR_CDOT_DB") and os.environ.get("HEARTVAR_VEP_FASTA")),
    reason="needs HEARTVAR_CDOT_DB + HEARTVAR_VEP_FASTA pointing at real data")
def test_bare_gene_input_lands_on_the_same_locus_as_the_supplied_transcript(
        monkeypatch):
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")
    from backend.clients import hgvs_resolver as hr
    hr._providers.cache_clear()

    for hgvs_c in ("c.1208G>A", "c.1207_1209del"):
        anchored = hr.resolve(hgvs_c, "NM_000257.4", "MYH7")
        bare = hr.resolve(hgvs_c, None, "MYH7")
        assert anchored is not None, f"transcript-anchored {hgvs_c} did not resolve"
        assert bare is not None, (
            f"bare-gene {hgvs_c} did not resolve — this is the case that used "
            "to go to rest.ensembl.org")
        assert (bare.chrom, bare.vep_pos, bare.vep_ref, bare.vep_alt) == \
               (anchored.chrom, anchored.vep_pos, anchored.vep_ref, anchored.vep_alt), (
            f"bare-gene {hgvs_c} resolved to a DIFFERENT locus than the "
            f"curator's transcript: {bare} vs {anchored}")
        assert bare.gnomad_id() == anchored.gnomad_id()
        assert "MANE Select" in (bare.note or "")
