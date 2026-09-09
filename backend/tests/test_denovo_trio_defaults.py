"""De novo status must not be silently discarded when Trio status is unset.

Reported from production 2026-09-07: a curator selected "De novo status:
Confirmed", left Trio status alone, and got NEITHER PS2 (+4) nor PM6 (+2) — a
silent zero, explained only by "De novo not confirmed by parental testing",
which never mentions the second dropdown that caused it.

Two causes, both fixed:
  * `trio_status` defaults to "" and the control LABELLED that "Not a trio",
    turning a default into a positive claim the curator never made. It now
    reads "Unknown".
  * `denovo_confirmed` required `trio_status == "trio"`, and PM6 eligibility
    required `trio in ("duo", "trio")`, so the default fell through both.

⚠ THE DUO CARVE-OUT IS LOAD-BEARING AND MUST NOT BE "SIMPLIFIED" AWAY. A duo
tests one parent and cannot establish both maternity and paternity, which ACMG
PS2 requires. A duo asserting "confirmed" previously received PS2_Strong with
the false evidence string "confirmed by parental testing (full trio)". An
EXPLICIT duo still routes to PM6: the specific statement beats the general one.
"""
from __future__ import annotations

import pytest

from backend.acmg.hard_coded import _eval_pm6, _eval_ps2

PTS = {"PS2": 4, "PM6": 2}


def _run(trio: str, denovo: str, inheritance: str = "AD"):
    """Mirror app.py's denovo_confirmed / denovo_trio_stated derivation."""
    confirmed = denovo == "confirmed" and trio != "duo"
    ps2 = _eval_ps2(confirmed, "TNNI3", 0, 0, trio_stated=(trio == "trio"))
    pm6 = _eval_pm6(trio, denovo, inheritance, confirmed,
                    ps2["status"] == "met", "TNNI3", 0)
    pts = sum(PTS[c] for c, e in (("PS2", ps2), ("PM6", pm6))
              if e["status"] == "met")
    return ps2, pm6, pts


def test_confirmed_de_novo_with_trio_unset_now_scores():
    """🔴 The regression test for the report. Was +0."""
    ps2, pm6, pts = _run("", "confirmed")
    assert ps2["status"] == "met"
    assert ps2["criteria_strength"] == "PS2_Strong"
    assert pts == 4


def test_absent_from_both_parents_scores_pm6_even_with_trio_unset():
    """CORRECTED 2026-09-07 (third revision, and the label change is why).

    The option now reads "ABSENT FROM BOTH PARENTS — parentage assumed". That
    is a statement that both parents WERE tested and the variant was not found,
    which is exactly ACMG PM6: "Assumed de novo, but without confirmation of
    paternity and maternity." Trio status adds nothing to it.

    The SVI clause "if the parents have not been tested for parentage or for
    the variant, no points should be awarded" is covered by the separate "Not
    assessed" option, which scores nothing — see the test below.

    Earlier revisions had this at +0 because the old label ("Unconfirmed /
    assumed") could be read as "we never tested the parents"."""
    ps2, pm6, pts = _run("", "unconfirmed")
    assert ps2["status"] == "not_met"
    assert pm6["criteria_strength"] == "PM6_Moderate"
    assert pts == 2


def test_not_assessed_is_the_svi_no_points_case():
    """🔴 The SVI clause lives HERE, not on the assumed-parentage option:
    "If the parents have not been tested for parentage or for the variant, no
    points should be awarded." "Not assessed" is that case."""
    for trio in ("", "duo", "trio"):
        ps2, pm6, pts = _run(trio, "")
        assert pts == 0, trio
        assert ps2["status"] == "not_met" and pm6["status"] == "not_met"


def test_a_duo_contradicts_absent_from_both_parents_and_scores_nothing():
    """A duo says one parent was tested; the de-novo option says the variant is
    absent from BOTH. The contradiction is surfaced, not scored."""
    for dn in ("unconfirmed", "confirmed"):
        ps2, pm6, pts = _run("duo", dn)
        assert pts == 0, dn
        assert "one parent" in pm6["evidence"].lower(), dn


def test_no_de_novo_assertion_still_scores_nothing():
    """Non-vacuous: the fix must not make an UNSET de novo status fire."""
    assert _run("", "")[2] == 0
    assert _run("trio", "")[2] == 0
    assert _run("duo", "")[2] == 0


def test_an_explicit_duo_earns_no_de_novo_criterion_and_says_why():
    """🔴 Do not "fix" this to PS2 OR to PM6. One parent cannot establish both
    maternity and paternity (ACMG PS2: "Confirmation of paternity only is
    insufficient"), and PM6 equally requires the variant to be absent from
    BOTH parents. The curator must be told, not shown a silent zero."""
    ps2, pm6, pts = _run("duo", "confirmed")
    assert ps2["status"] == "not_met"
    assert pm6["status"] == "not_met"
    assert pts == 0
    ev = pm6["evidence"].lower()
    assert "one parent" in ev
    assert "inherited from the untested parent" in ev


def test_a_full_trio_is_unchanged():
    ps2, _, pts = _run("trio", "confirmed")
    assert ps2["criteria_strength"] == "PS2_Strong"
    assert pts == 4


def test_ps2_does_not_claim_a_trio_that_was_never_reported():
    """🔴 The original defect was the STRING, not the points: a duo was told
    "confirmed by parental testing (full trio)"."""
    ps2, _, _ = _run("", "confirmed")
    ev = ps2["evidence"].lower()
    assert "full trio" not in ev
    assert "not specified" in ev or "curator's assertion" in ev


def test_ps2_does_claim_the_trio_when_one_was_reported():
    ps2, _, _ = _run("trio", "confirmed")
    assert "full trio" in ps2["evidence"].lower()


def test_the_trio_control_default_is_labelled_unknown():
    """"" is the absence of an answer. Labelling it "Not a trio" turned a
    default into an assertion, and put a self-contradictory pair on screen
    beside "Confirmed de novo"."""
    from pathlib import Path
    html = (Path(__file__).resolve().parents[2] / "index.html").read_text()
    i = html.index('id="hvl-trio"')
    block = html[i:i + 400]
    assert '<option value="">Unknown</option>' in block, block[:200]
    assert "Not a trio" not in block


@pytest.mark.parametrize("inheritance", ["AD", "AR", "XL", ""])
def test_the_fix_holds_across_inheritance_modes(inheritance):
    assert _run("", "confirmed", inheritance)[2] == 4
