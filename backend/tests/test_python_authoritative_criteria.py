"""The five criteria the LLM stopped being asked for on 2026-09-08.

PM1, PP1, BS4, PS1 and PM5 are now Python-authoritative. Two of them (PS1, PM5)
were already being discarded at the merge, so sending them was pure cost. The
other three were measured: the model reproduced the Python verdict in
substance, differing mostly by emitting a bare code where Python emits an
equal-valued strength suffix (PM1 against PM1_Moderate, both +2), and PM1 and
PP1 moved between two IDENTICAL AI runs. The model was adding run-to-run
variance without adding information.

Two invariants are pinned here, both through the real endpoint rather than
against the helpers, because the wiring is where this can break:

  1. **Nothing the model says about the five can move the score.** Proven
     differentially: the same request and the same evidence, run twice against
     two different model outputs, must produce IDENTICAL rows for all five and
     an identical points_total. A recount would be a second implementation of
     the thing under test, so this compares one run against the other.

  2. **PP1 and BS4 still fire on the AI arm, from the structured curator
     fields.** This is what stops the change being a silent capability loss:
     dropping the five from the prompt without ALSO merging Python's PM1/PP1/BS4
     into the AI arm leaves them at the merge placeholder, which measured 4
     points to 0 on a variant whose PM1 fires.

The free-text route is deliberately given up and is asserted gone in (2c): prose
in the family-history box no longer produces PP1 or BS4 in any mode. The numeric
fields are the only route, which is why index.html says so next to the inputs.

No network and no Anthropic calls: the evidence gather and ``stream_claude`` are
both replaced with canned equivalents, following test_curate_interpret_endpoint.
"""
from __future__ import annotations

import json

import pytest
from starlette.testclient import TestClient

import backend.app as app_module
from backend import ratelimit

SERVER_OWNED = ("PS1", "PM5", "PM1", "PP1", "BS4")

ACCOUNT = {
    "oid": "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee",
    "preferred_username": "curator@hospital.org.au",
    "name": "Test Curator",
}


def _model_json(criteria: list[dict]) -> str:
    return json.dumps({
        "confidence": "moderate",
        "summary": "Canned interpretation for the test.",
        "criteria": criteria,
        "gene_context": {},
    })


def _c(code: str, strength: str | None, status: str = "met") -> dict:
    return {
        "code": code,
        "status": status,
        "criteria_strength": strength,
        "direction": "benign" if code.startswith("B") else "pathogenic",
        "evidence": f"Model-asserted {code}.",
    }


BASELINE_CRITERIA = [_c("PS4", "PS4_Supporting")]

GREEDY_CRITERIA = BASELINE_CRITERIA + [
    _c("PS1", "PS1_Strong"),
    _c("PM5", "PM5_Moderate"),
    _c("PM1", "PM1_Strong"),
    _c("PP1", "PP1_Strong"),
    _c("BS4", "BS4"),
]


@pytest.fixture(autouse=True)
def _clean_state(monkeypatch):
    original_gather = app_module._gather_evidence_sse
    original_stream = app_module.stream_claude
    monkeypatch.setattr(ratelimit, "_ai_quota", {"day": None, "counts": {}})
    monkeypatch.setattr(app_module, "SERVER_AI_AVAILABLE", True)
    monkeypatch.setattr(app_module, "session_account", lambda request: dict(ACCOUNT))
    yield
    app_module._gather_evidence_sse = original_gather
    app_module.stream_claude = original_stream


@pytest.fixture(autouse=True)
def _no_inbound_rate_limit(monkeypatch):
    """Same reasoning as test_curate_interpret_endpoint: TestClient presents as
    one IP, so the 20/minute ceiling on /api/curate/stream is a budget shared
    across the whole suite and exhausting it fails a LATER file with a 429."""
    monkeypatch.setattr(ratelimit.limiter, "enabled", False)
    yield


def _install_gather(evidence: dict, gene: str = "MYH7"):
    async def _fake_gather(req, state):
        state["evidence"] = dict(evidence)
        state["variant_id"] = "14-23424115-G-A"
        state["gene"] = req.gene or gene
        state["hgvs_c"] = req.hgvs_c
        state["timing"] = {}
        return
        yield
    app_module._gather_evidence_sse = _fake_gather


def _install_claude(criteria: list[dict]):
    text = _model_json(criteria)

    async def _fake_stream(prompt, max_tokens=None):
        yield text
    app_module.stream_claude = _fake_stream


def _parse_sse(body: str) -> list[tuple[str, dict]]:
    events = []
    for block in body.split("\n\n"):
        event = data = None
        for line in block.splitlines():
            if line.startswith("event: "):
                event = line[len("event: "):].strip()
            elif line.startswith("data: "):
                data = line[len("data: "):]
        if event and data is not None:
            try:
                events.append((event, json.loads(data)))
            except json.JSONDecodeError:
                events.append((event, {}))
    return events


def _run_ai(model_criteria: list[dict], evidence: dict | None = None,
            gene: str = "MYH7", **body_overrides) -> dict:
    """POST an AI-mode curation and return the stage2_complete payload."""
    _install_gather(evidence or {}, gene=gene)
    _install_claude(model_criteria)
    body = {"gene": gene, "hgvs_c": "c.1988G>A", "ai_mode": "server"}
    body.update(body_overrides)
    resp = TestClient(app_module.app).post("/api/curate/stream", json=body)
    assert resp.status_code == 200, resp.text
    events = _parse_sse(resp.text)
    payloads = [d for e, d in events if e == "stage2_complete"]
    assert len(payloads) == 1, [e for e, _ in events]
    return payloads[0]


def _by_code(payload: dict) -> dict[str, dict]:
    return {c["code"]: c for c in payload["criteria"]}


def test_the_ai_copies_of_the_five_change_nothing():
    """The differential. Same request, same evidence, two model outputs."""
    quiet = _run_ai(BASELINE_CRITERIA)
    greedy = _run_ai(GREEDY_CRITERIA)

    assert greedy["points_total"] == quiet["points_total"], (
        "the model moved the score with a criterion it no longer owns: "
        f"{quiet['points_total']} -> {greedy['points_total']}"
    )
    assert greedy["classification"] == quiet["classification"]

    q, g = _by_code(quiet), _by_code(greedy)
    for code in SERVER_OWNED:
        assert g[code] == q[code], (
            f"{code} differs between the two runs, so the model's copy reached "
            f"the merge: {q[code]!r} -> {g[code]!r}"
        )


def test_the_five_never_arrive_with_an_ai_source():
    """Provenance, not just points. A criterion the curator reads as "ai" was
    written by the model; these five are Python's, and the report has to say so
    even when the Python rule stayed silent."""
    by = _by_code(_run_ai(GREEDY_CRITERIA))
    for code in SERVER_OWNED:
        assert by[code]["source"] != "ai", (
            f"{code} is attributed to the AI on the AI arm"
        )
        assert "Model-asserted" not in (by[code]["evidence"] or ""), (
            f"{code} carries the model's evidence string"
        )


def test_a_surviving_code_still_gets_through():
    """Control. If the filter were dropping everything, the tests above would
    pass for the wrong reason. Strength is deliberately not asserted, because
    the PS4 spec cap rewrites it and that belongs to PS4's own tests."""
    by = _by_code(_run_ai(GREEDY_CRITERIA))
    assert by["PS4"]["status"] == "met"
    assert by["PS4"]["source"] == "ai"


def test_pp1_fires_on_the_ai_arm_from_the_carrier_count():
    """Two affected relatives tested and carrying the variant, model silent.

    Without app._python_authoritative_supplementary this is where PP1 vanished:
    the AI's copy is dropped by the HARD_CODED_CRITERIA_CODES filter and nothing
    supplies a replacement, so the merge placeholder lands instead.
    """
    silent = _run_ai(BASELINE_CRITERIA)
    with_pp1 = _run_ai(BASELINE_CRITERIA, seg_affected_carriers=2)

    pp1 = _by_code(with_pp1)["PP1"]
    assert pp1["status"] == "met", "PP1 did not fire from the carrier count"
    assert pp1["criteria_strength"] == "PP1_Supporting"
    assert pp1["source"] == "inferred", (
        "source must stay 'inferred'. hard_coded.apply_cross_criterion_"
        "exclusions exempts the curator-derived PP1 from the literature "
        "hardening on exactly that value"
    )
    assert with_pp1["points_total"] == silent["points_total"] + 1


def test_bs4_fires_on_the_ai_arm_from_the_noncarrier_count():
    """Two affected relatives tested and NOT carrying it. Biesecker 2024 prices
    a non-segregation at -4.0, so this is the single biggest thing the curator
    can enter, and it is worth pinning that the AI arm still collects it."""
    silent = _run_ai(BASELINE_CRITERIA)
    with_bs4 = _run_ai(BASELINE_CRITERIA, seg_affected_noncarriers=2)

    bs4 = _by_code(with_bs4)["BS4"]
    assert bs4["status"] == "met", "BS4 did not fire from the non-carrier count"
    assert bs4["criteria_strength"] == "BS4"
    assert bs4["source"] == "inferred"
    assert with_bs4["points_total"] == silent["points_total"] - 4


def test_family_history_prose_alone_produces_neither():
    """The capability deliberately given up.

    An LLM BS4 read out of the family-history free text used to score -4 and
    could return "Likely benign" on prose alone. It cannot now: BS4 reads
    seg_affected_noncarriers and nothing else. PP1's prose route was already
    dead, killed by the literature hardening for want of a genotyped-carrier
    fact, so only BS4 loses anything measurable.

    This is the assertion behind the UI hint on the two numeric inputs in
    index.html. If it ever fails, the hint has become a lie.
    """
    prose = ("Two affected sisters were tested for the variant and neither "
             "carries it; three affected cousins do carry it.")
    silent = _run_ai(BASELINE_CRITERIA)
    typed = _run_ai(
        BASELINE_CRITERIA + [_c("PP1", "PP1_Strong"), _c("BS4", "BS4")],
        family=prose,
    )
    by = _by_code(typed)
    assert by["PP1"]["status"] != "met", "prose produced PP1"
    assert by["BS4"]["status"] != "met", "prose produced BS4"
    assert typed["points_total"] == silent["points_total"]


def test_pp1_and_bs4_still_cancel_when_both_counts_are_entered():
    """The mutual exclusion has to keep working now that both arms reach it by
    the same route. One pedigree cannot both support and refute segregation."""
    by = _by_code(_run_ai(
        BASELINE_CRITERIA,
        seg_affected_carriers=2, seg_affected_noncarriers=2,
    ))
    assert by["PP1"]["status"] == "not_met"
    assert by["BS4"]["status"] == "not_met"


def test_the_two_arms_agree_on_the_five():
    """The whole justification for the change: PM1/PP1/BS4 are one derivation
    now, so the AI arm and the evidence-only arm must return the same rows for
    them off the same inputs. Before, they were two implementations and the
    LLM's oscillated."""
    _install_gather({})
    body = {"gene": "MYH7", "hgvs_c": "c.1988G>A", "ai_mode": "none",
            "seg_affected_carriers": 2}
    resp = TestClient(app_module.app).post("/api/curate/stream", json=body)
    assert resp.status_code == 200, resp.text
    no_ai = [d for e, d in _parse_sse(resp.text) if e == "db_only_complete"]
    assert len(no_ai) == 1
    no_ai_by = {c["code"]: c for c in no_ai[0]["criteria"]}

    ai_by = _by_code(_run_ai(BASELINE_CRITERIA, seg_affected_carriers=2))

    for code in SERVER_OWNED:
        for field in ("status", "criteria_strength", "evidence"):
            if no_ai_by[code]["status"] == "not_assessed":
                continue
            assert no_ai_by[code][field] == ai_by[code][field], (
                f"{code}.{field} differs between the arms: "
                f"{no_ai_by[code][field]!r} vs {ai_by[code][field]!r}"
            )


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
