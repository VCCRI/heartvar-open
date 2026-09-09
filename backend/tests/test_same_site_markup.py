"""STRUCTURAL invariants for the same-site surfaces — where markup is emitted.

READ THIS BEFORE ADDING TO THIS FILE. Everything here is a substring search
over the TEXT of heartvar.js. That is useful for pinning *where* a block is
emitted (inside the existing gnomAD dropdown rather than a new section; the
same-residue block after the PM5 callout) and useless for anything about
whether the code RUNS. All eight assertions below passed while the page was
blank, because six unqualified ``esc(...)`` calls threw
``ReferenceError: esc is not defined`` at render time — valid syntax, correct
text, dead page.

Behavioural coverage belongs in ``test_frontend_render_executes.py``, which
evaluates heartvar.js in node and calls the renderers. Add assertions there,
not here, unless what you are pinning is genuinely the LOCATION of markup.
"""
from __future__ import annotations

import pathlib

import pytest

_REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
_JS = _REPO_ROOT / "static" / "heartvar.js"


@pytest.fixture(scope="module")
def js() -> str:
    return _JS.read_text(encoding="utf-8")


def test_same_site_rows_are_rendered_inside_the_existing_gnomad_dropdown(js):
    """No new section: the rows are spread into the gnomAD kvGroup (variant
    present) and appended to the absent branch, using the same kv() grammar
    as every other row in that dropdown."""
    assert "const _sameSiteRows = () =>" in js
    assert "..._sameSiteRows()," in js
    assert "+ kvGroup(_sameSiteRows());" in js
    gnomad_row = js.index("rows.push({ source: 'gnomad'")
    helper = js.index("const _sameSiteRows = () =>")
    assert helper < gnomad_row, "helper must be defined before the gnomAD row"


def test_same_site_panel_distinguishes_base_from_codon(js):
    """The three senses of "same site" carry different weight, so the panel
    must not collapse them into one undifferentiated list."""
    assert "Other alleles at this base" in js
    assert "Elsewhere in this codon" in js
    assert "not rsID" in js


def test_off_panel_degrades_visibly_rather_than_showing_nothing(js):
    """"Not available" and "None reported" are different findings; conflating
    them would let a missing DB read as an empty site."""
    assert "gnomad_available" in js
    assert "Not available —" in js
    assert "None reported in gnomAD v4" in js


def test_exactly_one_residue_callout_is_emitted(js):
    """D1. There used to be TWO identically-styled callouts here and they
    contradicted each other: one headed "PM5 — known pathogenic variant(s) at
    residue N" listing records, and one immediately beneath declaring those
    same records "not PM5 evidence". They are merged into a single block with a
    neutral heading, so the source must emit one and only one.
    """
    assert "let sameResidueBlock = ''" in js
    assert "clinvar_same_residue" in js
    assert '-callout-title">PM5' not in js
    assert "${pm5Callout}" not in js
    body = js[js.index("${section2Body}"):]
    inner = body[:body.index("clinvarBrowseLink")]
    assert inner.count("${sameResidueBlock}") == 1


def test_residue_callout_heading_makes_no_claim(js):
    """The heading must describe a list, not assert evidence — a curator
    glancing at the panel must not be able to read it as a verdict."""
    assert "ClinVar records at ${cap}" in js
    assert "PM5 row of the Criteria tab" in js


def test_residue_callout_has_a_fallback_and_an_empty_state(js):
    """"No records at this residue", "per-record verdicts unavailable" and
    "block missing entirely" are three different findings."""
    assert "eligibility was not available" in js
    assert "No other ClinVar record is recorded at this residue" in js


def test_every_same_residue_row_carries_a_pm5_eligibility_verdict(js):
    """A synonymous or VUS record must never be readable as PM5 support."""
    assert "pm5_eligible" in js
    assert "not PM5 evidence" in js
    assert "pm5_ineligible_reason" in js
    assert "DIFFERENT MISSENSE" in js
    assert 'ev-cls-ben">not PM5 evidence' not in js


def test_pm5_callout_names_both_residue_numbers_when_they_differ(js):
    """The candidates are selected in MANE numbering; the header used to be
    stamped with the picked transcript's model alone."""
    assert "matched_protein_position" in js
    assert "numbering_differs" in js
    assert "on ClinVar's MANE transcript" in js


def test_non_mane_transcript_warning_is_rendered_with_the_transcript_rows(js):
    assert "const nonManeRow =" in js
    assert "Numbering caveat" in js
    assert "Not MANE Select" in js
    assert "vep.ok && !vep.is_mane_select && !vep.is_mane_clinical)" in js
    assert js.count("            nonManeRow,") == 2


def test_landscape_stops_hiding_records_at_the_probands_own_codon(js):
    """showVus defaulted to false and the consequence chip defaulted to the
    proband's own bucket, which is how a VUS missense and a synonymous record
    at the proband's own codon both came to be invisible."""
    assert "const _vusAtCodon =" in js
    assert "showVus: _vusAtCodon" in js
    assert "_buildLandscapeBody(lsData, defaultCsq, _vusAtCodon)" in js
    assert "showVus: false }" not in js, "stale hard-coded default still present"
