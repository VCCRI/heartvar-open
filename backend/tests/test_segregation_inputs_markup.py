"""Both directions of segregation evidence must have a form input.

WHY THIS EXISTS. The BS4 input was silently lost. `index.html` kept the box for
"Affected relatives who carry the variant" (PP1's evidence, pathogenic
direction) while the box for affected relatives who do NOT carry it (BS4's
evidence, benign direction) was removed and `static/heartvar.js` hardcoded the
value to 0.

That was invisible three ways over:
  * `backend/models.py` defaults the field to 0, so the request stayed valid;
  * `backend/app.py` read it happily, so nothing errored;
  * `prompt.py` named that count as BS4's PRIMARY input and told the model to
    cite it, so the primary path could never fire, and any BS4 the model
    inferred from free text fired at the full -4 with no number behind it.

UPDATED 2026-09-08. The counts are now the ONLY route. PM1, PP1, BS4, PS1 and
PM5 became Python-authoritative and left the prompt entirely, because the model
reproduced the Python verdict in substance without adding information, while
moving between two identical runs. BS4's free-text route really
is gone (prose alone used to score -4 and return "Likely benign"); PP1's was
already dead, blocked by the literature hardening for want of a genotyped-carrier
count. So an empty box is now silently worth zero evidence, and the form has to
say so. That is the last group of tests below.

There is no browser in the test run, so the structural precondition is what
gets tested. A parallel case and its reasoning: test_header_layout_markup.py.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

INDEX = Path(__file__).resolve().parents[2] / "index.html"
HEARTVAR_JS = Path(__file__).resolve().parents[2] / "static" / "heartvar.js"


@pytest.fixture(scope="module")
def html():
    return INDEX.read_text()


@pytest.fixture(scope="module")
def js():
    return HEARTVAR_JS.read_text()


SEG_INPUTS = {
    "hvl-seg-affected-carriers": "PP1",
    "hvl-seg-affected-noncarriers": "BS4",
}


@pytest.mark.parametrize("input_id,criterion", sorted(SEG_INPUTS.items()))
def test_the_segregation_input_exists(html, input_id, criterion):
    assert f'id="{input_id}"' in html, (
        f"{input_id} is gone from index.html, so {criterion} has no structured "
        f"input. This is how BS4 broke: the field keeps working end-to-end "
        f"because models.py defaults it to 0, so nothing fails loudly."
    )


@pytest.mark.parametrize("input_id,criterion", sorted(SEG_INPUTS.items()))
def test_the_segregation_input_has_a_label(html, input_id, criterion):
    assert f'for="{input_id}"' in html, f"{input_id} has no <label>"


@pytest.mark.parametrize("input_id,criterion", sorted(SEG_INPUTS.items()))
def test_the_value_is_read_from_the_form_not_hardcoded(js, input_id, criterion):
    """🔴 The actual defect. The input can exist while the JS ignores it."""
    assert f"_segInt('{input_id}')" in js, (
        f"heartvar.js does not read {input_id}. If it has been replaced with a "
        f"literal (the BS4 bug was `const seg_affected_noncarriers = 0;`), the "
        f"box renders and does nothing."
    )


def test_neither_segregation_count_is_assigned_a_literal(js):
    """Catches the regression in the general form, whatever the value."""
    for field in ("seg_affected_carriers", "seg_affected_noncarriers"):
        bad = re.search(rf"const\s+{field}\s*=\s*\d+\s*;", js)
        assert not bad, (
            f"{field} is assigned the literal {bad.group(0)!r} — it must come "
            f"from the form, or the criterion it feeds cannot fire."
        )


def test_the_two_directions_are_symmetric(html):
    """The point of the fix. Collecting only the pathogenic-direction count
    biases the tool, because a curator whose family does NOT segregate has no
    way to say so."""
    carriers = html.count('id="hvl-seg-affected-carriers"')
    noncarriers = html.count('id="hvl-seg-affected-noncarriers"')
    assert carriers == noncarriers == 1, (
        f"asymmetric segregation inputs: {carriers} carrier box(es), "
        f"{noncarriers} non-carrier box(es)"
    )


DENOVO_COUNTS = {
    "hvl-denovo-confirmed-count": "parentage confirmed",
    "hvl-denovo-unconfirmed-count": "parentage assumed",
}


@pytest.mark.parametrize("input_id,qualifier", sorted(DENOVO_COUNTS.items()))
def test_the_de_novo_counts_stay_split(html, input_id, qualifier):
    """Both boxes must exist, and each must say which kind of occurrence it
    counts. They are NOT interchangeable and must not be merged into one field.

    SVI Table 1 prices a confirmed occurrence at twice an assumed one and
    Table 2 turns the total into a strength, so the split is worth a full rung:
    on PTPN11, one confirmed occurrence is PS2_Strong where one assumed
    occurrence is PS2_Moderate; on KCNQ1 it is Moderate against Supporting.
    Merging the boxes silently downgrades or upgrades every RASopathy and
    KCNQ1 de-novo call. Verified by calling _denovo_points_tier directly.
    """
    assert f'id="{input_id}"' in html, (
        f"{input_id} is gone. Merging the two de-novo occurrence counts loses "
        f"the SVI 2x weighting for confirmed parentage."
    )
    label = re.search(rf'<label for="{input_id}">(.*?)</label>', html, re.S)
    assert label, f"{input_id} has no <label>"
    assert qualifier in label.group(1), (
        f"{input_id}'s label no longer says '{qualifier}', so the curator "
        f"cannot tell the two occurrence boxes apart"
    )


def test_the_de_novo_status_options_name_the_parentage_qualifier(html):
    """The dropdown must keep BOTH absent-from-both-parents options and say
    which is which. Confirmed routes to PS2 (+4), assumed to PM6 (+2); a single
    merged option would price every untested pedigree as though parentage had
    been established."""
    for value, qualifier in (("confirmed", "parentage confirmed"),
                             ("unconfirmed", "parentage assumed")):
        m = re.search(rf'<option value="{value}">(.*?)</option>', html, re.S)
        assert m, f"the de novo status option '{value}' is gone"
        text = m.group(1)
        assert "Absent from both parents" in text, (
            f"the '{value}' option no longer says the variant is absent from "
            f"both parents"
        )
        assert qualifier in text, (
            f"the '{value}' option does not say '{qualifier}', so PS2 (+4) and "
            f"PM6 (+2) are indistinguishable in the form"
        )
