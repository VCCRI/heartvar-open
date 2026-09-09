"""Unit tests for the deterministic (server-owned) PS1 criterion.

PS1 was moved from AI-evaluated to server-owned (2026-07-15): the AI fired it
0× even though ``get_pm5_evidence`` already surfaces the same-amino-acid
candidates (``ps1_candidates`` — a DIFFERENT nucleotide encoding the SAME AA
change, the proband's own record removed) off the same ClinVar retrieval PM5
uses. ``_clinvar_ps1_criterion`` now derives PS1 directly from that surface.

These tests pin every branch of the criterion (≥2★ gate, one-star /
no-candidate / not-applicable paths, integration via
``compute_hard_coded_criteria``) plus the clinvar.py fail-closed
self-exclusion: when the proband's own HGVS-c can't be resolved we cannot prove
a same-AA record isn't the proband itself, so the PS1 surface must stay empty
(the key catastrophic-crossing mitigation).
"""
from __future__ import annotations

import asyncio
import contextlib
import os
import sqlite3
import tempfile
from pathlib import Path

import backend.clients.clinvar as clinvar
from backend.acmg.hard_coded import (
    _clinvar_ps1_criterion,
    compute_hard_coded_criteria,
)


def _cand(name, cs="Pathogenic", stars=2):
    return {"name": name, "clinical_significance": cs, "stars": stars, "tier": "P"}


def _ev(ps1_candidates, *, ok=True, not_applicable=False, pos=403, gene="MYH7"):
    """Build the get_pm5_evidence-shaped blob. PS1 reads only the ps1_* surface;
    the PM5 `candidates` side is left empty so the two criteria don't cross."""
    d = {
        "ok": ok,
        "gene": gene,
        "protein_position": pos,
        "candidates": [],
        "count": 0,
        "count_two_star": 0,
        "ps1_candidates": ps1_candidates,
        "ps1_count": len(ps1_candidates),
        "ps1_count_two_star": sum(
            1 for c in ps1_candidates if (c.get("stars") or 0) >= 2
        ),
    }
    if not_applicable:
        d["not_applicable"] = True
    return d


def test_ps1_met_on_two_star_candidate():
    ev = _ev([_cand("NM_000257.4(MYH7):c.1207C>T (p.Arg403Trp)", stars=3)])
    crit = _clinvar_ps1_criterion(ev, "MYH7")
    assert crit["status"] == "met"
    assert crit["criteria_strength"] == "PS1_Strong"
    assert crit["source"] == "hard_coded"
    assert "residue 403" in crit["evidence"]


def test_ps1_not_met_when_only_one_star():
    ev = _ev([_cand("NM_000257.4(MYH7):c.1207C>T (p.Arg403Trp)", stars=1)])
    crit = _clinvar_ps1_criterion(ev, "MYH7")
    assert crit["status"] == "not_met"
    assert crit["criteria_strength"] is None
    assert "below the PS1 confidence bar" in crit["evidence"]


def test_ps1_not_met_when_no_candidate():
    crit = _clinvar_ps1_criterion(_ev([]), "MYH7")
    assert crit["status"] == "not_met"
    assert "No same-amino-acid P/LP record" in crit["evidence"]


def test_ps1_not_applicable_passthrough():
    for ev in ({"ok": False}, _ev([], not_applicable=True), None):
        crit = _clinvar_ps1_criterion(ev, "MYH7")
        assert crit["status"] == "not_met"
        assert crit["criteria_strength"] is None


def test_ps1_fires_on_single_rasopathy_candidate():
    ev = _ev(
        [_cand("NM_002834.5(PTPN11):c.211T>C (p.Phe71Leu)", stars=2)],
        gene="PTPN11", pos=71,
    )
    crit = _clinvar_ps1_criterion(ev, "PTPN11")
    assert crit["status"] == "met"
    assert crit["criteria_strength"] == "PS1_Strong"


def test_compute_hard_coded_emits_ps1():
    ev = {
        "clinvar_pm5_candidates": _ev(
            [_cand("NM_000257.4(MYH7):c.1207C>T (p.Arg403Trp)", stars=3)]
        )
    }
    out = compute_hard_coded_criteria(ev, {"inheritance_input": "AD"}, "MYH7")
    ps1 = [c for c in out if c["code"] == "PS1"]
    assert len(ps1) == 1
    assert ps1[0]["status"] == "met"
    assert ps1[0]["criteria_strength"] == "PS1_Strong"
    assert ps1[0]["source"] == "hard_coded"


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


def _row(variation_id, name, cs="Pathogenic",
         review_status="reviewed by expert panel"):
    return (
        variation_id, "single nucleotide variant", name, "MYH7", cs,
        review_status, 5, "", "NC_000014.9", 1, 1, "G", "A", "GRCh38",
        "2024-01-01",
    )


@contextlib.contextmanager
def _temp_db(rows):
    fd, path = tempfile.mkstemp(suffix=".db", prefix="clinvar_ps1_test_")
    os.close(fd)
    conn = sqlite3.connect(path)
    try:
        conn.execute(_DDL)
        conn.executemany(_INSERT_SQL, rows)
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


_SAME_AA_ROW = _row(
    variation_id=14099,
    name="NM_000257.4(MYH7):c.1207C>T (p.Arg403Trp)",
)


def test_ps1_surface_fires_for_resolved_proband():
    """Positive control: a DIFFERENT nucleotide (proband c.1209G>A) leaves the
    comparator in the PS1 surface."""
    with _temp_db([_SAME_AA_ROW]):
        res = asyncio.run(
            clinvar.get_pm5_evidence(
                "MYH7", 403, proband_alt_aa="Trp", proband_hgvs_c="c.1209G>A",
            )
        )
    assert res["ok"] is True
    assert res["ps1_count_two_star"] == 1
    assert res["ps1_candidates"][0]["name"].endswith("(p.Arg403Trp)")


def test_ps1_surface_fails_closed_when_proband_hgvs_unresolved():
    """When the proband's own HGVS-c is unresolved (None) we cannot prove the
    same-AA record isn't the proband itself → it is dropped and PS1 cannot
    fire off it."""
    with _temp_db([_SAME_AA_ROW]):
        res = asyncio.run(
            clinvar.get_pm5_evidence(
                "MYH7", 403, proband_alt_aa="Trp", proband_hgvs_c=None,
            )
        )
    assert res["ok"] is True
    assert res["ps1_candidates"] == []
    assert res["ps1_count_two_star"] == 0
    assert _clinvar_ps1_criterion(res, "MYH7")["status"] == "not_met"


if __name__ == "__main__":
    import sys
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn(); print(f"PASS {name}")
            except AssertionError as exc:
                failures += 1; print(f"FAIL {name}: {exc}")
    sys.exit(1 if failures else 0)
