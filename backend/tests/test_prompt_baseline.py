"""Tests that pin the production prompt, byte-for-byte.

``build_prompt`` and ``claude._system_blocks`` between them decide everything
the model is told. A silent change to either moves classifications without
failing anything else, so these tests hash both and assert the enumerations
that tell the model which criteria are its job.

History worth knowing: this file used to test a reversible ``prompt_variant``
A/B switch that served a slimmed prompt alongside the production one. The
variant never produced a reported number — every banked benchmark row carries
the baseline value — and it was removed on 2026-09-09
along with ``_SLIM_SYSTEM_REMOVALS`` and ``_slim_system_text``. The baseline
invariants it guarded are worth keeping, so they live on here.
"""

import hashlib
import re

from backend import claude
from backend.prompt import (
    SYSTEM_PROMPT,
    _OUTPUT_SCHEMA_BLOCK,
    build_prompt,
    format_precomputed_criteria_block,
)


_EVIDENCE = {
    "vep": {
        "ok": True,
        "queried_as": "MYH7:c.1988G>A",
        "selected_transcript_id": "ENST00000355349",
        "selected_transcript_source": "MANE Select",
        "is_mane_select": True,
        "transcript_id": "ENST00000355349",
        "assembly_name": "GRCh38",
        "seq_region_name": "14",
        "start": 23424029,
        "allele_string": "C/T",
        "hgvsc": "ENST00000355349.4:c.1988G>A",
        "hgvsp": "ENSP00000347507.3:p.Arg663His",
        "most_severe_consequence": "missense_variant",
        "impact": "MODERATE",
        "sift_prediction": "deleterious",
        "sift_score": 0.0,
        "polyphen_prediction": "probably_damaging",
        "polyphen_score": 0.99,
        "cadd_phred": 32.0,
        "revel_score": 0.91,
    },
}


_CLINICAL_CONTEXT = {
    "zygosity": "het",
    "inheritance_input": "AD",
    "proband_sex": "female",
}


_HARD_CODED = [
    {"code": "PM2", "status": "met", "criteria_strength": "PM2",
     "evidence": "Absent from gnomAD"},
    {"code": "PP3", "status": "met", "criteria_strength": "PP3",
     "evidence": "REVEL 0.91, CADD 32"},
]


def _call():
    """Build the consolidated prompt from the shared fixture."""
    kwargs = dict(
        clinical_context=_CLINICAL_CONTEXT,
        hard_coded_criteria=_HARD_CODED,
        segregation_context="3 affected relatives carry the variant.",
    )
    return build_prompt(
        "MYH7", "c.1988G>A", "HP:0001639", "Autosomal dominant",
        "Cardiac lesion: HCM. Extra-cardiac anomalies: none.",
        _EVIDENCE,
        **kwargs,
    )


_PYTHON_AUTHORITATIVE_CODES = ("PS1", "PM5", "PM1", "PP1", "BS4")


_AI_LIST = "PS3, PS4, PP2, PP4, BS3, BP1, BP5"


_C5B_HARD_CODED = [
    {"code": "PM2", "status": "met", "criteria_strength": "PM2",
     "evidence": "Absent from gnomAD v4 (allele frequency 0)."},
    {"code": "BA1", "status": "not_met", "criteria_strength": "",
     "evidence": "Allele frequency 0.0001 is below the BA1 threshold of 0.05."},
    {"code": "PVS1", "status": "not_met", "criteria_strength": "",
     "evidence": "Missense variant; PVS1 is reserved for predicted LoF."},
]


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def test_system_blocks_is_one_cached_block_of_exactly_system_prompt():
    """_system_blocks() returns EXACTLY one cached text block holding the
    unchanged SYSTEM_PROMPT with ephemeral cache_control. _CACHE_TTL is patched
    to the production-unset state so a local HEARTVAR_CACHE_TTL cannot mask a
    regression."""
    prev = claude._CACHE_TTL
    try:
        claude._CACHE_TTL = ""
        assert claude._system_blocks() == [
            {
                "type": "text",
                "text": SYSTEM_PROMPT,
                "cache_control": {"type": "ephemeral"},
            }
        ]
    finally:
        claude._CACHE_TTL = prev


def test_system_and_user_turn_are_stable_by_sha256():
    """The load-bearing production invariant, by SHA-256: the cached system text
    is exactly SYSTEM_PROMPT, and the assembled user turn is reproducible across
    calls from the same fixture. A prompt edit is meant to fail this — read the
    diff and re-bank deliberately rather than deleting the assertion."""
    prev = claude._CACHE_TTL
    try:
        claude._CACHE_TTL = ""
        assert _sha(claude._system_blocks()[0]["text"]) == _sha(SYSTEM_PROMPT)
    finally:
        claude._CACHE_TTL = prev

    assert _sha(_call()) == _sha(_call())
    assert _call().endswith(_OUTPUT_SCHEMA_BLOCK)


def test_no_variant_asks_the_model_for_a_server_owned_code():
    """The five codes that stopped being sent to the LLM on 2026-09-08.

    PS1 and PM5 were already DISCARDED at the merge (app.py drops any AI
    criterion whose code is in HARD_CODED_CRITERIA_CODES), so asking for them was
    pure token waste. PM1, PP1 and BS4 joined them because the model reproduced
    Python in substance while oscillating between two identical runs.

    Checked in the two places the prompt states the model's job: the AI-EVALUATED
    enumeration and the ordered `criteria` list. Both must name exactly the
    surviving 7. Mentions of a removed code elsewhere are fine and several are
    deliberate, because surviving gates say what a criterion is NOT (the PS3
    gate's "that is PM5 / PM1 evidence, NEVER PS3"), which is why this pins the
    two enumerations rather than grepping the file.

    The enumeration is in SYSTEM_PROMPT; the ordered list is in
    _OUTPUT_SCHEMA_BLOCK, which the user turn appends.
    """
    text = claude._system_blocks()[0]["text"]
    assert f"\n    {_AI_LIST}\n" in text, (
        "the AI-EVALUATED enumeration is not the 7 survivors"
    )
    assert f"EXACTLY 7 entries, in this order: {_AI_LIST}." in _OUTPUT_SCHEMA_BLOCK, (
        "the ordered `criteria` list is not the 7 survivors"
    )


def test_every_precomputed_enumeration_names_all_five():
    """The other half: the model has to be TOLD the five are precomputed.

    Silence is worse than the old wording. A code that vanishes from the
    AI-EVALUATED list without appearing in the precomputed list reads as an
    oversight and the model may emit it anyway.

    The prompt states the precomputed set THREE times (the master
    HYBRID-EVALUATION SCOPE list, the do-NOT-evaluate list, and the field-rules
    reminder in the migrated schema block), so this asserts every occurrence of
    the old 16-code tail carries the 5-code extension. A half-finished edit
    leaves one enumeration short and nothing else catches it.
    """
    tail = "PM3, BP2, PP5, BP6"
    extended = tail + ", PM5, PS1, PM1, PP1, BS4"
    text = claude._system_blocks()[0]["text"]
    assert text.count(extended) == text.count(tail) > 0, (
        f"{text.count(tail) - text.count(extended)} precomputed "
        "enumeration(s) still stop at BP6"
    )
    assert "16 precomputed" not in text, "a stale 'the 16 precomputed' survived"
    assert "12 interpretive" not in text, "a stale 'the 12 interpretive' survived"


def test_precomputed_block_keeps_notmet_prose():
    """NOT-MET codes carry their full prose rationale, not a bare code. The
    model grounds its summary in these lines."""
    base = format_precomputed_criteria_block(_C5B_HARD_CODED)
    assert "BA1 NOT MET — Allele frequency 0.0001" in base
    assert "PVS1 NOT MET — Missense variant" in base


def test_the_prompt_variant_switch_has_not_come_back():
    """Regression guard for the 2026-09-09 removal.

    The slim A/B switch was removed because it had never produced a reported
    number, and the cost of leaving it in was a second prompt path that nothing
    exercised. If it is genuinely wanted again, restore it from history and
    delete this test in the same commit — do not let it reappear by accident on
    the back of a resync.
    """
    import backend.prompt as prompt_mod

    assert not hasattr(prompt_mod, "_slim_system_text")
    assert not hasattr(prompt_mod, "_SLIM_SYSTEM_REMOVALS")
    for fn in (build_prompt, format_precomputed_criteria_block, claude._system_blocks):
        assert "prompt_variant" not in fn.__code__.co_varnames, (
            f"{fn.__name__} has regained a prompt_variant parameter"
        )
