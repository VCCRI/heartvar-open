"""AI-free supplementary criterion inference.

Best-effort, deterministic inference of the judgement criteria derivable from
structured inputs / already-computed evidence (PP1/BS4/PM1/PS1/PM5), plus the
no-key 28-code assembler.

``build_no_ai_criteria`` is still evidence-only (``ai_mode="none"``). Since
2026-09-08 ``infer_supplementary_criteria`` is NOT: the AI path calls it too,
through ``app._python_authoritative_supplementary``, for PM1/PP1/BS4 only.
Those three stopped being sent to the model because it reproduced the Python
verdict in substance without adding information, while moving between identical
runs. So the derivations below are now the single source of
truth for both arms, and a change here changes the AI result too.

Moved VERBATIM from ``backend.app`` (no logic change). Imports only from
``.constants`` / ``.hard_coded`` / stdlib — never from ``backend.app``.
"""
from __future__ import annotations

from .constants import _CRITERION_NAMES, CANONICAL_CRITERIA_ORDER
from .hard_coded import (
    _criterion_applicable,
    _pm1_assessment_from_domain_plp,
    _bs4_strength_from_noncarriers,
    _pp1_strength_from_lod_meioses,
)


_SERVER_OWNED_NOT_ASSESSED: dict[str, str] = {
    "PM1": (
        "Not assessed. HeartVar derives PM1 in Python from the ClinVar "
        "domain P/LP density, which was not available for this variant (or "
        "the gene's VCEP marks PM1 Not Applicable). An API key will not "
        "produce it."
    ),
    "PP1": (
        "Not assessed. HeartVar derives PP1 in Python from the structured "
        "segregation counts. Enter the number of affected relatives tested "
        "for the variant and found to CARRY it. Prose in the family-history "
        "box does not produce PP1, and neither does an API key."
    ),
    "BS4": (
        "Not assessed. HeartVar derives BS4 in Python from the structured "
        "segregation counts. Enter the number of affected relatives tested "
        "for the variant and found NOT to carry it. Prose in the "
        "family-history box does not produce BS4, and neither does an API key."
    ),
}


def server_owned_placeholder(code: str) -> dict:
    """The entry for a Python-authoritative criterion whose derivation did not
    fire on the inputs supplied.

    PM1, PP1 and BS4 became server-owned on 2026-09-08 and are no longer sent to
    the model. ``infer_supplementary_criteria`` emits an entry only when the
    derivation FIRES, so without this the two arms described the same absence
    differently: the no-AI arm said ``not_assessed`` with copy naming the input
    that would produce the criterion, while the AI arm said ``not_met`` with a
    generic "the rule did not fire". ``not_met`` also overclaims — the rule was
    not evaluated and found absent, it had nothing to evaluate.

    Both arms now use this. ``not_assessed`` is ignored by compute_points_total
    and by apply_cross_criterion_exclusions, so it cannot perturb a score.
    """
    return {
        "code": code,
        "name": _CRITERION_NAMES.get(code, code),
        "status": "not_assessed",
        "direction": "benign" if code.startswith("B") else "pathogenic",
        "criteria_strength": None,
        "evidence": _SERVER_OWNED_NOT_ASSESSED.get(code) or (
            "Not assessed — HeartVar derives this criterion in Python and the "
            "inputs it needs were not supplied."
        ),
        "source": "unassessed",
    }


def build_no_ai_criteria(
    hard_coded: list[dict], supplementary: list[dict] | None = None,
) -> list[dict]:
    """Build the full 28-code criteria list for the no-key (AI-free) path.

    Starts from the deterministic hard-coded set, layers any
    ``supplementary`` entries inferred without an LLM (Phase 2 —
    PP1/BS4/PM1/PS1/PM5 best-effort), then fills every remaining code with
    an explicit ``not_assessed`` placeholder. ``not_assessed`` is distinct
    from ``not_met`` on purpose: it tells the curator the criterion was
    NOT evaluated (it needs clinical-geneticist / AI judgement), rather
    than evaluated and found absent. Both compute_points_total (sums only
    ``met``) and apply_cross_criterion_exclusions (acts only on ``met``)
    ignore it, so it never perturbs the deterministic score."""
    by_code: dict[str, dict] = {c["code"]: c for c in hard_coded}
    for c in supplementary or []:
        code = c.get("code") or ""
        if not code or code in by_code:
            continue
        by_code[code] = c
    ordered: list[dict] = []
    for code in CANONICAL_CRITERIA_ORDER:
        entry = by_code.get(code)
        if entry is None:
            entry = {
                "code": code,
                "name": _CRITERION_NAMES.get(code, code),
                "status": "not_assessed",
                "direction": "benign" if code.startswith("B") else "pathogenic",
                "criteria_strength": None,
                "evidence": _SERVER_OWNED_NOT_ASSESSED.get(code) or (
                    "Not assessed — this criterion requires AI / "
                    "clinical-geneticist judgement. Add an API key for a "
                    "complete interpretation."
                ),
                "source": "unassessed",
            }
        ordered.append(entry)
    return ordered


def _inferred_entry(
    code: str, status: str, strength: str | None, evidence_txt: str,
) -> dict:
    """A supplementary (no-key, AI-free) criterion dict. source="inferred"
    so the audit trail distinguishes it from both the deterministic
    hard-coded set and the AI set."""
    return {
        "code": code,
        "name": _CRITERION_NAMES.get(code, code),
        "status": status,
        "direction": "benign" if code.startswith("B") else "pathogenic",
        "criteria_strength": strength,
        "evidence": evidence_txt,
        "source": "inferred",
    }


def infer_supplementary_criteria(
    evidence: dict,
    clinical_context: dict,
    gene: str | None,
    hard_coded: list[dict],
) -> list[dict]:
    """Best-effort, AI-FREE inference of the judgement criteria that can be
    derived from structured inputs / already-computed evidence. Returns
    criterion dicts for whichever of PP1/BS4/PM1/PS1/PM5 can be inferred. Codes
    that need genuine literature reading (PS3/BS3, PS4, PP4, PP2, BP1, BP5) are
    intentionally omitted and surface as ``not_assessed``.

    TWO callers since 2026-09-08, and they take different slices. The no-key
    path takes the whole list (build_no_ai_criteria then refuses to let PS1/PM5
    clobber the hard-coded ones). The AI path takes PM1/PP1/BS4 ONLY, filtered
    by app._python_authoritative_supplementary, because
    compute_hard_coded_criteria already owns PS1/PM5 and
    merge_hard_coded_and_ai has no such guard.

    Two evidence classes:
      - PP1 / BS4 are driven purely by the curator's structured segregation
        fields, so they fire ONLY on explicit input (and are inert in the
        eRepo benchmark, which carries no family counts).
      - PM1 / PS1 / PM5 are driven by already-computed variant/gene evidence
        (domain-P/LP density; ClinVar same-AA / same-codon at ≥2★), so they
        DO appear in the benchmark and are gated there.
    """
    evidence = evidence or {}
    cc = clinical_context or {}
    out: list[dict] = []

    seg_aff_carriers = int(cc.get("seg_affected_carriers") or 0)
    seg_meioses = int(cc.get("seg_meioses") or 0)
    tier = _pp1_strength_from_lod_meioses(
        None, seg_meioses, seg_aff_carriers, 0, gene=gene,
    )
    if tier:
        out.append(_inferred_entry(
            "PP1", "met", f"PP1_{tier}",
            f"Co-segregation: {seg_aff_carriers} affected carrier(s)"
            + (f" across {seg_meioses} informative meioses" if seg_meioses else "")
            + f" — {tier} (curator-entered).",
        ))

    seg_aff_noncarriers = int(cc.get("seg_affected_noncarriers") or 0)
    _bs4_tier = _bs4_strength_from_noncarriers(seg_aff_noncarriers, gene=gene)
    if _bs4_tier is not None:
        out.append(_inferred_entry(
            "BS4", "met", f"BS4_{_bs4_tier}" if _bs4_tier else "BS4",
            f"Non-segregation: {seg_aff_noncarriers} affected relative(s) "
            f"tested and do NOT carry the variant — "
            + (f"{gene} VCEP publishes BS4 at {_bs4_tier} for this count"
               if _bs4_tier else
               "BS4 Strong (-4.0 points per Biesecker 2024 ClinGen PP1/BS4 "
               "guidance)")
            + "; curator-entered.",
        ))

    pm1_payload = evidence.get("domain_plp")
    if isinstance(pm1_payload, dict) and _criterion_applicable(gene, "PM1"):
        pm1_status, pm1_sentence = _pm1_assessment_from_domain_plp(pm1_payload)
        if pm1_status == "met":
            out.append(_inferred_entry("PM1", "met", "PM1_Moderate", pm1_sentence))
        elif pm1_status == "supporting":
            out.append(_inferred_entry("PM1", "met", "PM1_Supporting", pm1_sentence))
        elif pm1_status == "not_met":
            out.append(_inferred_entry("PM1", "not_met", None, pm1_sentence))

    pm5_blob = evidence.get("clinvar_pm5_candidates")
    if isinstance(pm5_blob, dict) and pm5_blob.get("ok"):
        ps1_two_star = int(pm5_blob.get("ps1_count_two_star") or 0)
        pm5_two_star = int(pm5_blob.get("count_two_star") or 0)
        if ps1_two_star >= 1:
            out.append(_inferred_entry(
                "PS1", "met", "PS1_Strong",
                f"{ps1_two_star} ClinVar ≥2★ P/LP record(s) with the same "
                "amino-acid change (different nucleotide) — PS1.",
            ))
        elif pm5_two_star >= 1:
            out.append(_inferred_entry(
                "PM5", "met", "PM5_Moderate",
                f"{pm5_two_star} ClinVar ≥2★ P/LP record(s) with a different "
                "missense at the same codon — PM5.",
            ))

    return out
