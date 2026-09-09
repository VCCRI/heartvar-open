"""Shape tests for backend.clients.clinvar.

clinvar.py is a LOCAL SQLite client (reads data/clinvar.db), NOT an HTTP
client — there is nothing to mock with httpx. Each test builds a temporary
on-disk SQLite DB with the REAL ``variants`` schema (from
scripts/build_clinvar_db.py), points ``clinvar.DB_PATH`` at it, drives the
async public functions with ``asyncio.run``, then restores ``DB_PATH`` and
removes the temp file in a try/finally.

Covered surfaces (distinct from test_input_and_phenotype_parsing.py, which
covers parse_variant_input + the pure _split_phenotypes placeholder filter):
  - fetch_clinvar end-to-end result/record shape, star mapping, VCV padding
  - _hgvs_match word-boundary discrimination via fetch_clinvar
  - the "N conditions" placeholder filter END-TO-END through fetch_clinvar
    and get_gene_phenotype_strings
  - db-missing error path

Runnable via pytest (``pytest backend/tests/test_clinvar.py``) or directly
(``python -m backend.tests.test_clinvar``).
"""
from __future__ import annotations

import asyncio
import contextlib
import sqlite3
import tempfile
from pathlib import Path

import backend.clients.clinvar as clinvar

_DDL = """
CREATE TABLE variants (
    variation_id            INTEGER,
    obj_type                TEXT,
    name                    TEXT,
    gene_symbol             TEXT,
    clinical_significance   TEXT,
    review_status           TEXT,
    number_submitters       INTEGER,
    phenotype_list          TEXT,
    chromosome_accession    TEXT,
    start                   INTEGER,
    stop                    INTEGER,
    reference_allele        TEXT,
    alternate_allele        TEXT,
    assembly                TEXT,
    last_evaluated          TEXT
)
"""

_INSERT_SQL = "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"


def _row(
    *,
    variation_id,
    obj_type="single nucleotide variant",
    name,
    gene_symbol,
    clinical_significance="Pathogenic",
    review_status="criteria provided, single submitter",
    number_submitters=1,
    phenotype_list="",
    chromosome_accession="NC_000014.9",
    start=1,
    stop=1,
    reference_allele="G",
    alternate_allele="A",
    assembly="GRCh38",
    last_evaluated="2024-01-01",
):
    """Build a positional row tuple matching the DDL column order."""
    return (
        variation_id, obj_type, name, gene_symbol, clinical_significance,
        review_status, number_submitters, phenotype_list, chromosome_accession,
        start, stop, reference_allele, alternate_allele, assembly, last_evaluated,
    )


@contextlib.contextmanager
def _temp_db(rows):
    """Create a temp .db with the real variants schema + ``rows`` inserted,
    point clinvar.DB_PATH at it, then restore DB_PATH and delete the file.

    DB_PATH is a module-level attribute read at call time inside the sync
    helpers, so rebinding the module attribute is enough to redirect them.
    """
    fd, path = tempfile.mkstemp(suffix=".db", prefix="clinvar_test_")
    import os
    os.close(fd)
    conn = sqlite3.connect(path)
    try:
        conn.execute(_DDL)
        if rows:
            conn.executemany(_INSERT_SQL, rows)
        conn.commit()
    finally:
        conn.close()

    original = clinvar.DB_PATH
    clinvar.DB_PATH = Path(path)
    try:
        yield Path(path)
    finally:
        clinvar.DB_PATH = original
        with contextlib.suppress(OSError):
            os.unlink(path)


def test_fetch_clinvar_shape_and_stars():
    row = _row(
        variation_id=14093,
        name="NM_000257.4(MYH7):c.1208G>A (p.Arg403Gln)",
        gene_symbol="MYH7",
        clinical_significance="Pathogenic",
        review_status="reviewed by expert panel",
        number_submitters=23,
        phenotype_list="Hypertrophic cardiomyopathy|Cardiovascular phenotype",
    )
    with _temp_db([row]):
        result = asyncio.run(clinvar.fetch_clinvar("MYH7", "c.1208G>A"))

    assert result["ok"] is True
    assert result["found"] is True
    for key in (
        "uids", "records", "total_records", "total_submissions",
        "phenotype_matched_submissions", "all_conditions",
    ):
        assert key in result, f"missing top-level key {key!r}"

    assert result["total_records"] == 1
    assert result["total_submissions"] == 23
    assert result["phenotype_matched_submissions"] == 23
    assert result["uids"] == ["14093"]

    rec = result["records"][0]
    assert rec["stars"] == 3
    assert rec["accession"] == "VCV000014093"
    assert rec["uid"] == "14093"
    assert rec["variation_id"] == 14093
    assert rec["clinical_significance"] == "Pathogenic"
    assert rec["number_submitters"] == 23
    assert rec["title"] == "NM_000257.4(MYH7):c.1208G>A (p.Arg403Gln)"
    assert rec["conditions"] == [
        "Hypertrophic cardiomyopathy", "Cardiovascular phenotype",
    ]
    assert result["all_conditions"] == [
        "Hypertrophic cardiomyopathy", "Cardiovascular phenotype",
    ]


def test_fetch_clinvar_hgvs_word_boundary():
    rows = [
        _row(
            variation_id=1001,
            name="NM_000257.4(MYH7):c.123A>T (p.Lys41Asn)",
            gene_symbol="MYH7",
        ),
        _row(
            variation_id=1002,
            name="NM_000257.4(MYH7):c.1234A>T (p.Thr412Ser)",
            gene_symbol="MYH7",
        ),
    ]
    with _temp_db(rows):
        result = asyncio.run(clinvar.fetch_clinvar("MYH7", "c.123A>T"))

    assert result["ok"] is True
    assert result["total_records"] == 1
    assert result["records"][0]["variation_id"] == 1001
    assert result["records"][0]["title"].endswith("c.123A>T (p.Lys41Asn)")


def test_n_conditions_placeholder_filtered_end_to_end():
    row = _row(
        variation_id=2222,
        name="NM_000257.4(MYH7):c.5000G>A (p.Arg1667His)",
        gene_symbol="MYH7",
        phenotype_list="Hypertrophic cardiomyopathy|6 conditions|not provided",
    )
    with _temp_db([row]):
        result = asyncio.run(clinvar.fetch_clinvar("MYH7", "c.5000G>A"))
        phen_strings = asyncio.run(clinvar.get_gene_phenotype_strings("MYH7"))

    assert result["ok"] is True
    rec = result["records"][0]
    assert "6 conditions" not in rec["conditions"]
    assert "not provided" not in rec["conditions"]
    assert "Hypertrophic cardiomyopathy" in rec["conditions"]

    assert "6 conditions" not in result["all_conditions"]
    assert "Hypertrophic cardiomyopathy" in result["all_conditions"]

    assert "6 conditions" not in phen_strings
    assert "hypertrophic cardiomyopathy" in phen_strings
    assert "not provided" not in phen_strings


def test_fetch_clinvar_db_missing():
    missing = Path(tempfile.gettempdir()) / "clinvar_does_not_exist_xyz.db"
    with contextlib.suppress(OSError):
        missing.unlink()
    assert not missing.exists()

    original = clinvar.DB_PATH
    clinvar.DB_PATH = missing
    try:
        result = asyncio.run(clinvar.fetch_clinvar("MYH7", "c.1208G>A"))
    finally:
        clinvar.DB_PATH = original

    assert result["ok"] is False
    assert "error" in result
    assert "not found" in result["error"].lower()


def test_gene_landscape_positions_aggregation():
    rows = [
        _row(
            variation_id=1, name="NM_000257.4(MYH7):c.1208G>A (p.Arg403Gln)",
            gene_symbol="MYH7", clinical_significance="Pathogenic",
            review_status="criteria provided, single submitter",
            start=100,
        ),
        _row(
            variation_id=2, name="NM_000257.4(MYH7):c.1208G>T (p.Arg403Leu)",
            gene_symbol="MYH7", clinical_significance="Pathogenic",
            review_status="reviewed by expert panel",
            start=100,
        ),
        _row(
            variation_id=3, name="NM_000257.4(MYH7):c.2146G>A (p.Gly716Arg)",
            gene_symbol="MYH7", clinical_significance="Likely pathogenic",
            review_status="criteria provided, multiple submitters, no conflicts",
            start=200,
        ),
        _row(
            variation_id=4, name="NM_000257.4(MYH7):c.2998G>A (p.Glu1000Lys)",
            gene_symbol="MYH7", clinical_significance="Benign", start=300,
        ),
        _row(
            variation_id=5, name="NM_000257.4(MYH7):c.599T>C (p.Met200Thr)",
            gene_symbol="MYH7", clinical_significance="Uncertain significance",
            start=400,
        ),
        _row(
            variation_id=6, name="NM_000257.4(MYH7):c.1234-5A>G",
            gene_symbol="MYH7", clinical_significance="Uncertain significance",
            start=500,
        ),
    ]
    with _temp_db(rows):
        r = clinvar._gene_variant_landscape_sync("MYH7", None)

    assert r["ok"] is True
    assert r["total_classified"] == 6
    assert r["tier_counts"] == {"P": 2, "LP": 1, "VUS": 2, "LB": 0, "B": 1}

    assert r["tier_counts_by_csq"]["missense"] == {
        "P": 2, "LP": 1, "VUS": 1, "LB": 0, "B": 1,
    }
    assert r["tier_counts_by_csq"]["splice"] == {
        "P": 0, "LP": 0, "VUS": 1, "LB": 0, "B": 0,
    }

    positions = r["positions"]
    assert r["positions_truncated"] is False
    assert len(positions) == 5
    assert [p["gpos"] for p in positions] == [100, 200, 300, 400, 500]
    assert r["g_min"] == 100 and r["g_max"] == 500

    by_gpos = {p["gpos"]: p for p in positions}
    assert by_gpos[100]["tier"] == "P"
    assert by_gpos[100]["csq"] == "missense"
    assert by_gpos[100]["count"] == 2
    assert by_gpos[100]["stars"] == 3
    assert by_gpos[100]["aa"] == 403
    assert by_gpos[500]["csq"] == "splice"
    assert by_gpos[500]["aa"] is None


def test_classify_consequence_buckets():
    f = clinvar._classify_consequence
    assert f("NM_x(MYH7):c.1208G>A (p.Arg403Gln)") == "missense"
    assert f("NM_x(MYH7):c.1207C>T (p.Arg403Ter)") == "truncating"
    assert f("NM_x(MYH7):c.1000dupA (p.Lys334AsnfsTer5)") == "truncating"
    assert f("NM_x(MYH7):c.1209C>T (p.Arg403=)") == "synonymous"
    assert f("NM_x(MYH7):c.1234_1236del (p.Lys412del)") == "inframe"
    assert f("NM_x(MYH7):c.1234-5A>G") == "splice"
    assert f("NM_x(MYH7):c.1234+1G>A") == "splice"
    assert f("NM_x(MYH7):c.-30C>T") == "utr"
    assert f("NM_x(MYH7):c.*30C>T") == "utr"


def _run_all():
    test_fetch_clinvar_shape_and_stars()
    test_fetch_clinvar_hgvs_word_boundary()
    test_n_conditions_placeholder_filtered_end_to_end()
    test_fetch_clinvar_db_missing()
    test_gene_landscape_positions_aggregation()
    test_classify_consequence_buckets()
    print("all clinvar shape tests passed")


if __name__ == "__main__":
    _run_all()
