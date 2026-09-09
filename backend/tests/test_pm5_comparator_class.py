"""PM5 comparators must be MISSENSE changes.

THE BUG. ``_AA_CHANGE_RE`` is ``p\\.([A-Z][a-z]{2})(\\d+)([A-Z][a-z]{2})`` and
its alt-AA group matches ``Ter``. So a nonsense record — ``p.Tyr79Ter`` — parsed
cleanly, cleared the ``tier in ("P","LP")`` test and the ≥2★ review gate, and
entered ``candidates``: the list ``_clinvar_pm5_criterion`` counts to fire
PM5_Moderate (+2). The evidence string then described it as "3 other
Pathogenic/Likely-pathogenic **missense** change(s)", which is a false
statement about a stop-gain in a document that can reach a clinical report.

PM5 is "a different MISSENSE change at the same amino-acid residue". A
stop-gain at the residue is evidence about truncation (PVS1's argument), and
says nothing about whether that residue tolerates substitution.

SCALE (measured against data/clinvar.db, 2026-08-24, 677-gene cardiac panel):
6141 of 11991 (gene, residue) pairs carrying ≥2★ P/LP protein tokens have
**only** nonsense/synonymous comparators — 51.2%. MYBPC3 alone: 99 of 115
residues. So this was very likely the dominant term in PM5's measured 25%
precision.

The rows below are copied verbatim from data/clinvar.db. MYBPC3 residue 79
carries three ≥2★ Pathogenic stop-gains and zero missense records, so a
proband missense at Tyr79 was a guaranteed false PM5.
"""
from __future__ import annotations

import contextlib
import os
import sqlite3
import tempfile
from pathlib import Path

import backend.clients.clinvar as clinvar
from backend.acmg.hard_coded import _clinvar_pm5_criterion

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

_MULTI = "criteria provided, multiple submitters, no conflicts"
_PANEL = "reviewed by expert panel"
_SINGLE = "criteria provided, single submitter"

_ROWS = [
    ("NM_000256.3(MYBPC3):c.237C>G (p.Tyr79Ter)", "Pathogenic", _MULTI),
    ("NM_000256.3(MYBPC3):c.236dup (p.Tyr79Ter)", "Pathogenic", _MULTI),
    ("NM_000256.3(MYBPC3):c.237C>A (p.Tyr79Ter)", "Pathogenic", _MULTI),
    ("NM_000256.3(MYBPC3):c.235del (p.Tyr79Thrfs)", "Pathogenic", _MULTI),
    ("NM_000256.3(MYBPC3):c.237C>T (p.Tyr79Tyr)", "Pathogenic", _MULTI),
]

_ROWS_MISSENSE = [
    ("NM_000256.3(MYBPC3):c.394G>A (p.Ala132Thr)", "Pathogenic", _PANEL),
    ("NM_000256.3(MYBPC3):c.395C>T (p.Ala132Val)", "Pathogenic", _MULTI),
]


@contextlib.contextmanager
def _temp_db(rows):
    fd, path = tempfile.mkstemp(suffix=".db", prefix="clinvar_pm5csq_")
    os.close(fd)
    conn = sqlite3.connect(path)
    try:
        conn.execute(_DDL)
        conn.executemany(_INSERT_SQL, [
            (i + 1, "single nucleotide variant", name, "MYBPC3", sig, rev, 3,
             "Hypertrophic cardiomyopathy", "NC_000011.10", 47340000, 47340000,
             "na", "na", "GRCh38", "2026-01-01")
            for i, (name, sig, rev) in enumerate(rows)
        ])
        conn.commit()
    finally:
        conn.close()
    original = clinvar.DB_PATH
    clinvar.DB_PATH = Path(path)
    try:
        yield
    finally:
        clinvar.DB_PATH = original
        with contextlib.suppress(OSError):
            os.unlink(path)


def test_nonsense_records_never_enter_the_pm5_candidate_list():
    """This is the test that would have caught the bug. Pre-fix, `candidates`
    held all three p.Tyr79Ter rows."""
    with _temp_db(_ROWS):
        ev = clinvar._pm5_evidence_sync("MYBPC3", 79, proband_alt_aa="Cys",
                                        proband_hgvs_c="c.236A>G")
    names = [c["name"] for c in ev["candidates"]]
    assert names == [], (
        f"non-missense records offered as PM5 comparators: {names}"
    )
    assert ev["count"] == 0
    assert ev["count_two_star"] == 0


def test_frameshift_and_synonymous_are_excluded_too():
    with _temp_db(_ROWS):
        ev = clinvar._pm5_evidence_sync("MYBPC3", 79, proband_alt_aa="Cys",
                                        proband_hgvs_c="c.236A>G")
    excluded = {c["name"] for c in ev["non_missense_excluded"]}
    assert any("Thrfs" in n for n in excluded), (
        "a p.LeuNNNProfs-style frameshift name parses as a substitution and "
        f"must be excluded by consequence, not by regex accident: {excluded}"
    )
    assert any("Tyr79Tyr" in n for n in excluded)
    assert ev["non_missense_excluded_count"] == 5
    assert ev["non_missense_excluded_two_star_count"] == 5


def test_pm5_does_not_fire_on_a_nonsense_only_residue():
    with _temp_db(_ROWS):
        ev = clinvar._pm5_evidence_sync("MYBPC3", 79, proband_alt_aa="Cys",
                                        proband_hgvs_c="c.236A>G")
    crit = _clinvar_pm5_criterion(ev, "MYBPC3")
    assert crit["status"] == "not_met", crit["evidence"]
    assert crit["criteria_strength"] is None


def test_evidence_string_never_calls_a_stop_gain_a_missense_change():
    """The report-honesty half of the fix. The criterion row must not assert
    that a p.Tyr79Ter record is 'a Pathogenic/Likely-pathogenic missense
    change', and must not read as though the residue were simply empty."""
    with _temp_db(_ROWS):
        ev = clinvar._pm5_evidence_sync("MYBPC3", 79, proband_alt_aa="Cys",
                                        proband_hgvs_c="c.236A>G")
    text = _clinvar_pm5_criterion(ev, "MYBPC3")["evidence"]
    assert "Ter" not in text, f"stop-gain quoted as PM5 evidence: {text}"
    assert "non-missense" in text, (
        "a residue whose only P/LP records are stop-gains must not read "
        f"identically to an empty residue: {text}"
    )


def test_pm5_still_fires_on_a_genuine_missense_comparator():
    """Guard against over-filtering: the fix must not be a blanket PM5 kill."""
    with _temp_db(_ROWS_MISSENSE):
        ev = clinvar._pm5_evidence_sync("MYBPC3", 132, proband_alt_aa="Gly",
                                        proband_hgvs_c="c.395C>G")
    crit = _clinvar_pm5_criterion(ev, "MYBPC3")
    assert crit["status"] == "met", crit["evidence"]
    assert crit["criteria_strength"] == "PM5_Moderate"
    assert "Ala132Thr" in crit["evidence"] or "Ala132Val" in crit["evidence"]
