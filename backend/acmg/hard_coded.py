"""Deterministic ACMG criterion scorer (the "hard-coded" 21 criteria).

This module holds ``compute_hard_coded_criteria`` and its full helper closure,
moved VERBATIM from ``backend.app`` (no logic change). It owns:

  - input normalisation (``_normalize_criteria`` / ``_hard_coded_entry``);
  - the confidence-tiered ClinVar PP5/BP6 cluster;
  - the gnomAD / SpliceAI / VCEP-frequency extractors + thresholds;
  - the gene LoF-mechanism + PVS1 splice-strength logic;
  - the cross-criterion mutual-exclusion pass + the AI/hard-coded merge;
  - the free-text family-history classifier;
  - the protein-position + PM1-domain helpers shared with the evidence layer.

Imports only from ``.constants`` / clients / stdlib — never from ``backend.app``.
"""
from __future__ import annotations

import json
import logging
import os
import re
from pathlib import Path

from ..clients.clinvar import _classify_tier as _clinvar_classify_tier
from ..clients.panelapp import _categorize_panel, check_hpo_relevance
from .constants import (
    _BARE_CODE_POINTS,
    _CRITERION_NAMES,
    _CLINVAR_STAR_STRENGTH,
    CANONICAL_CRITERIA_ORDER,
    HARD_CODED_CRITERIA_CODES,
)

log = logging.getLogger("heartvar.acmg.hard_coded")


def _normalize_criteria(items, default_source: str | None = None) -> list[dict]:
    """Accept either the legacy dict-per-entry shape or the compact tuple
    shape ``[code, status, criteria_strength, evidence]`` — or its B1
    5-element extension ``[code, status, criteria_strength, evidence, facts]``
    — and return the full dict shape (code, name, status, direction,
    criteria_strength, evidence, facts, source) that downstream consumers +
    the frontend renderer already expect. ``direction`` is derived from the
    code prefix (B → benign, otherwise pathogenic). ``name`` is looked
    up in the server-side _CRITERION_NAMES map. ``source`` is preserved
    from the dict form when present; otherwise defaults to
    ``default_source`` (typically "ai" for AI-returned criteria and
    "hard_coded" for the deterministic Python-evaluated set). The field
    is stored for downstream auditability but is intentionally NOT
    surfaced to the UI.

    ``facts`` (B1) is the optional structured-evidence object the LLM emits
    for quantified criteria (e.g. ``{"proband_count": 12}`` for PS4) — the
    engine scores the raw number deterministically rather than re-parsing it
    from a truncated prose sentence. Only a dict is kept; anything else → None.

    Tolerant by design — passes through the ``facts`` key on dict entries and
    the 5th tuple element (forward compat), and handles short tuples
    (defensive against an early truncation of the array)."""
    out: list[dict] = []
    for it in items or []:
        if isinstance(it, dict):
            code = (it.get("code") or "").strip()
            name = it.get("name") or _CRITERION_NAMES.get(code, code)
            direction = it.get("direction") or (
                "benign" if code.startswith("B") else "pathogenic"
            )
            facts = it.get("facts")
            out.append({
                "code": code,
                "name": name,
                "status": it.get("status"),
                "direction": direction,
                "criteria_strength": it.get("criteria_strength"),
                "evidence": it.get("evidence") or "",
                "facts": facts if isinstance(facts, dict) else None,
                "source": it.get("source") or default_source,
            })
        elif isinstance(it, (list, tuple)):
            code = (it[0] if len(it) > 0 else "") or ""
            status = it[1] if len(it) > 1 else None
            strength = it[2] if len(it) > 2 else None
            evidence = it[3] if len(it) > 3 else ""
            facts = it[4] if len(it) > 4 else None
            out.append({
                "code": code,
                "name": _CRITERION_NAMES.get(code, code),
                "status": status,
                "direction": "benign" if code.startswith("B") else "pathogenic",
                "criteria_strength": strength,
                "evidence": evidence or "",
                "facts": facts if isinstance(facts, dict) else None,
                "source": default_source,
            })
    return out


def _hard_coded_entry(
    code: str, status: str, strength: str | None, evidence: str,
) -> dict:
    """Build a hard-coded criterion dict in the normalised shape.
    Always carries source="hard_coded"."""
    return {
        "code": code,
        "name": _CRITERION_NAMES.get(code, code),
        "status": status,
        "direction": "benign" if code.startswith("B") else "pathogenic",
        "criteria_strength": strength,
        "evidence": evidence,
        "source": "hard_coded",
    }


_PP1_TIER_ORDER = ("Supporting", "Moderate", "Strong")


_PP1_DEFAULT_LADDER = (3, 5, 7)

_PP1_PROBAND_INCLUSIVE = ("affected individual",)

_PP1_LADDER_RE = re.compile(
    r"(?:>=|≥|at least|more than|>|\+)?\s*(\d+)(?:\s*-\s*\d+)?\s+"
    r"(?:informative\s+)?"
    r"(segregations?|meios\w*|affected\s+(?:individuals?|family\s+members?|relatives?))",
    re.IGNORECASE,
)


_BS4_DEFAULT_MIN_NONCARRIERS = 2

_BS4_RUNG_RE = re.compile(
    r"absence of the variant in\s+(\d+)\s*(?:or more\s*)?affected\s+family\s+members?",
    re.IGNORECASE,
)


def _bs4_spec_rungs(gene: str | None) -> dict[str, int]:
    """{strength: minimum affected NON-carriers} from this gene's published BS4
    rows, strongest-first when iterated by the caller.

    Empty dict when the gene publishes no parseable count, in which case the
    caller applies ``_BS4_DEFAULT_MIN_NONCARRIERS`` at bare BS4 (-4.0), which is
    the pre-2026-09-07 behaviour for every gene.
    """
    strengths = (_crit_spec(gene, "BS4").get("strengths") or {})
    out: dict[str, int] = {}
    for label in ("Strong", "Moderate", "Supporting"):
        m = _BS4_RUNG_RE.search(strengths.get(label) or "")
        if m:
            out[label] = int(m.group(1))
    return out


def _bs4_strength_from_noncarriers(
    noncarriers, gene=None,
) -> str | None:
    """The BS4 strength `noncarriers` earns on this gene, or None.

    Returns a bare tier name ("Strong"/"Moderate"/"Supporting") for a gene that
    publishes rungs, or the sentinel "" meaning "bare BS4, -4.0" for the 26
    genes that publish only a Strong row (whose existing behaviour is to emit an
    unsuffixed BS4). Never raises.
    """
    try:
        n = int(noncarriers) if noncarriers is not None and not isinstance(
            noncarriers, bool) else 0
    except (TypeError, ValueError):
        n = 0
    rungs = _bs4_spec_rungs(gene)
    if len(rungs) > 1:
        for label in ("Strong", "Moderate", "Supporting"):
            if label in rungs and n >= rungs[label]:
                return label
        return None
    return "" if n >= _BS4_DEFAULT_MIN_NONCARRIERS else None


def _pp1_spec_ladder(gene: str | None) -> tuple[int, int, int]:
    """(supporting, moderate, strong) PP1 thresholds for `gene`, in AFFECTED
    RELATIVES (proband excluded).

    Falls back to the SVI/Kelly 3/5/7 default when the gene publishes no
    parseable numbers — which is the same ladder 26 of 27 genes publish anyway,
    so the fallback is the rule rather than a guess.

    The number must be FOLLOWED by its counted noun for the regex to take it.
    That is load-bearing: KCNQ1's row contains "Schwartz score ... >3" before
    "the proband + 3 affected family members", and a bare-number pattern would
    read the Schwartz threshold as a segregation count.
    """
    strengths = (_crit_spec(gene, "PP1").get("strengths") or {})
    found: dict[str, int] = {}
    for label in ("Supporting", "Moderate", "Strong"):
        m = _PP1_LADDER_RE.search(strengths.get(label) or "")
        if not m:
            continue
        n, unit = int(m.group(1)), re.sub(r"\s+", " ", m.group(2).lower())
        if any(u in unit for u in _PP1_PROBAND_INCLUSIVE):
            n -= 1
        found[label] = n
    if len(found) != 3:
        return _PP1_DEFAULT_LADDER
    out = (found["Supporting"], found["Moderate"], found["Strong"])
    return out if out[0] <= out[1] <= out[2] else _PP1_DEFAULT_LADDER


def _pp1_strength_from_lod_meioses(
    lod=None, meioses=None, carriers=None, unaff=None, gene=None,
) -> str | None:
    """Single source of truth for the PP1 co-segregation strength ladder,
    shared by the no-key path (no_ai.infer_supplementary_criteria) and the
    literature-sourced PP1 hardening in apply_cross_criterion_exclusions — so
    the two never drift.

    Three axes, each graded against the gene's OWN published thresholds
    (``_pp1_spec_ladder``, normalised to affected relatives with the proband
    excluded). The result is the STRONGEST tier any available axis supports:
      - explicit LOD (preferred): ≥2.1 Strong / ≥1.5 Moderate / ≥0.9 Supporting
      - informative meioses:      the gene's ladder, e.g. 7 / 5 / 3
      - affected carriers:        the same ladder, in the same unit
    Requires ≥1 genotyped affected carrier to fire at all (returns None below
    that).

    🔴 WHY THE AXES ARE MAXED AND NOT CHAINED, and it is the whole point of the
    2026-09-07 change. Carriers used to be capped at Supporting no matter how
    many there were, so 11 affected relatives scored the same +1 as 2. The
    Cardiomyopathy VCEP publishes its ladder in SEGREGATIONS ("≥7 segregations
    (LOD 2.1) for STRONG"), and a count of affected relatives carrying the
    variant IS that count, so refusing to grade it discarded the evidence.
    The cost of that: PP1 came out Supporting on every variant that carried it,
    where the specs' own thresholds would have given Strong or Moderate, and
    variants sat at VUS instead of LP as a result.

    An earlier attempt did this as ``n = meioses if meioses else carriers`` and
    was refuted for NON-MONOTONIC SCORING: 7 carriers with the meiosis field
    blank scored Strong, while the same 7 carriers with ``meioses=1`` scored
    zero — so filling in an optional field truthfully LOWERED the score by 4
    points. Taking the max over the axes makes that impossible by construction:
    supplying one more input can never reduce the tier.

    ⚠ The direction of that choice, stated rather than buried: where a curator
    gives a carrier count and a SMALLER meiosis count, this scores the larger.
    The two are readings of one pedigree and the VCEP ladders publish no rule
    for reconciling them, so the count of observations is used. The
    unaffected-carrier temper below still applies and still lowers the tier —
    that is contrary evidence on a different axis, not the same axis refined.

    Unaffected-carrier temper: if ≥1 genotyped UNAFFECTED carrier is reported,
    drop one tier (min Supporting); if that would zero the tier, keep Supporting
    only when ≥2 affected carriers, else None. Returns a bare tier string
    ("Supporting"/"Moderate"/"Strong") or None. Never raises."""
    def _i(v):
        try:
            return int(v) if v is not None and not isinstance(v, bool) else 0
        except (TypeError, ValueError):
            return 0

    def _f(v):
        try:
            return float(v) if v is not None and not isinstance(v, bool) else None
        except (TypeError, ValueError):
            return None

    carriers, meioses, unaff, lod = _i(carriers), _i(meioses), _i(unaff), _f(lod)
    if carriers < 1:
        return None

    sup_n, mod_n, strong_n = _pp1_spec_ladder(gene)

    def _by_count(n: int) -> str | None:
        """Grade one count against this gene's published thresholds."""
        if n >= strong_n:
            return "Strong"
        if n >= mod_n:
            return "Moderate"
        if n >= sup_n:
            return "Supporting"
        return None

    candidates: list[str] = []
    if lod is not None:
        if lod >= 2.1:
            candidates.append("Strong")
        elif lod >= 1.5:
            candidates.append("Moderate")
        elif lod >= 0.9:
            candidates.append("Supporting")
    for _t in (_by_count(meioses), _by_count(carriers)):
        if _t:
            candidates.append(_t)
    if not candidates and carriers >= 2:
        candidates.append("Supporting")
    if not candidates:
        return None
    tier = max(candidates, key=_PP1_TIER_ORDER.index)
    if unaff >= 1:
        idx = _PP1_TIER_ORDER.index(tier) - 1
        if idx < 0:
            return "Supporting" if carriers >= 2 else None
        tier = _PP1_TIER_ORDER[idx]
    return tier


def _clinvar_assertion_source(rec: dict, stars: int) -> str:
    """Transparency prefix naming the ClinVar source + star level for a
    PP5/BP6 evidence string. Built from the exact-variant record."""
    acc = rec.get("accession") or rec.get("uid") or "ClinVar"
    sig = rec.get("clinical_significance") or "classification"
    rs = rec.get("review_status") or ""
    n_sub = rec.get("number_submitters")
    rs_part = f": '{rs}'" if rs else ""
    sub = f"; {n_sub} submitter(s)" if n_sub else ""
    return (
        f"ClinVar {acc} reports {sig} for this exact variant "
        f"({stars}★{rs_part}){sub}"
    )


def _select_clinvar_assertion_record(
    clinvar_ev: dict | None, directions: set[str],
) -> tuple[dict | None, int, str | None]:
    """From the proband's exact-variant ClinVar records, pick the
    highest-confidence record whose classified tier is in *directions*
    (``{"P","LP"}`` for PP5, ``{"B","LB"}`` for BP6). Selection is by max
    stars, then by submitter count — records[0] is the highest-submitter
    match, NOT the highest-star, so a 3★ expert-panel record must not be
    masked by a higher-submitter 1★ aggregate. Returns (record, stars, tier)
    or (None, 0, None) when no same-direction record exists."""
    cv = clinvar_ev or {}
    recs = cv.get("records") or []
    best: tuple[dict, int, str] | None = None
    best_key: tuple[int, int] = (-1, -1)
    for r in recs:
        tier = _clinvar_classify_tier(r.get("clinical_significance"))
        if tier not in directions:
            continue
        stars = r.get("stars") or 0
        key = (stars, r.get("number_submitters") or 0)
        if key > best_key:
            best_key = key
            best = (r, stars, tier)
    if best is None:
        return None, 0, None
    return best


def _clinvar_no_apply_reason(
    clinvar_ev: dict | None, directions: set[str], code: str,
) -> str:
    """Explanatory evidence string for a not_met PP5/BP6 when an exact-variant
    record exists but does not qualify (wrong direction, conflicting, or
    below the 1★ confidence threshold)."""
    cv = clinvar_ev or {}
    if not cv.get("ok", True):
        return "ClinVar lookup unavailable — PP5/BP6 not evaluated."
    recs = cv.get("records") or []
    if not recs:
        return (
            "No exact-variant ClinVar assertion found for this variant — "
            "PP5/BP6 not applicable."
        )
    rec, stars, _tier = _select_clinvar_assertion_record(cv, directions)
    label = (
        "Pathogenic/Likely pathogenic" if code == "PP5"
        else "Benign/Likely benign"
    )
    if rec is not None:
        return (
            f"{_clinvar_assertion_source(rec, stars)} — below the 1★ "
            f"confidence threshold for {code} (unreviewed / no-assertion "
            "record); not applied."
        )
    top = recs[0]
    top_sig = top.get("clinical_significance") or "—"
    if _clinvar_classify_tier(top_sig) == "VUS":
        return (
            f"ClinVar reports '{top_sig}' for this exact variant "
            f"(uncertain / conflicting) — {code} not applicable."
        )
    return (
        f"No {label} ClinVar assertion for this exact variant "
        f"(top: '{top_sig}') — {code} not applicable."
    )


def _clinvar_pp5_bp6_criteria(
    clinvar_ev: dict | None, *, ba1_met: bool, bs1_met: bool,
) -> list[dict]:
    """Derive confidence-tiered PP5 and BP6 entries from the proband's
    exact-variant ClinVar record. PP5 is withheld (and the conflict surfaced
    in its evidence string) when the engine fired BA1/BS1 — population
    frequency contradicts a pathogenic ClinVar assertion. If the variant has
    BOTH a qualifying pathogenic AND a qualifying benign ClinVar record (rare
    cross-record conflict), neither is applied and the conflict is surfaced.
    The further guards against strong NON-frequency evidence — BP6 vs PVS1/PS*
    and PP5 vs BS2/BS3/BS4 — run post-merge in apply_cross_criterion_exclusions
    (those codes are computed elsewhere / after merge and not visible here).
    Returns [PP5_entry, BP6_entry]."""
    cv = clinvar_ev or {}

    rec, stars, _tier = _select_clinvar_assertion_record(cv, {"P", "LP"})
    if rec is not None and stars in _CLINVAR_STAR_STRENGTH:
        tier_name = _CLINVAR_STAR_STRENGTH[stars]
        src = _clinvar_assertion_source(rec, stars)
        if ba1_met or bs1_met:
            which = "BA1" if ba1_met else "BS1"
            pp5 = _hard_coded_entry(
                "PP5", "not_met", None,
                f"{src}, but HeartVar fired {which} (gnomAD population "
                "frequency above the benign threshold) — population frequency "
                "contradicts the ClinVar assertion; PP5 withheld and conflict "
                "surfaced for curator review.",
            )
        else:
            pp5 = _hard_coded_entry(
                "PP5", "met", f"PP5_{tier_name}",
                f"{src} — ClinVar review status maps to PP5 {tier_name}, but "
                "shown for reference only and excluded from the point total "
                "(reusing ClinVar's own classification would be circular).",
            )
    else:
        pp5 = _hard_coded_entry(
            "PP5", "not_met", None,
            _clinvar_no_apply_reason(cv, {"P", "LP"}, "PP5"),
        )

    rec, stars, _tier = _select_clinvar_assertion_record(cv, {"B", "LB"})
    if rec is not None and stars in _CLINVAR_STAR_STRENGTH:
        tier_name = _CLINVAR_STAR_STRENGTH[stars]
        src = _clinvar_assertion_source(rec, stars)
        bp6 = _hard_coded_entry(
            "BP6", "met", f"BP6_{tier_name}",
            f"{src} — ClinVar review status maps to BP6 {tier_name}, but "
            "shown for reference only and excluded from the point total "
            "(reusing ClinVar's own classification would be circular).",
        )
    else:
        bp6 = _hard_coded_entry(
            "BP6", "not_met", None,
            _clinvar_no_apply_reason(cv, {"B", "LB"}, "BP6"),
        )

    if pp5["status"] == "met" and bp6["status"] == "met":
        conflict = (
            "Conflicting high-confidence ClinVar assertions for this exact "
            "variant — both Pathogenic/Likely pathogenic and Benign/Likely "
            "benign records present; neither PP5 nor BP6 applied. Conflict "
            "surfaced for curator review."
        )
        pp5 = _hard_coded_entry("PP5", "not_met", None, conflict)
        bp6 = _hard_coded_entry("BP6", "not_met", None, conflict)

    return [pp5, bp6]


_PM5_ALT_AA_RE = re.compile(r"p\.[A-Z][a-z]{2}\d+([A-Z][a-z]{2})")

_FUNCTIONAL_ASSAY_RE = re.compile(
    r"functional|assay|in[\s-]*vitro|in[\s-]*vivo|electrophysiolog|"
    r"patch[\s-]*clamp|biochemical|enzymatic|knock[\s-]*(?:in|out)|"
    r"minigene|splic\w*\s+assay|expression\s+(?:assay|study|studies)|"
    r"transfect|western\s+blot|reporter\s+assay"
)

_SPLICING_ASSAY_RE = re.compile(
    r"minigene|splic|rna\s+assay|rt[\s-]*pcr|exon\s+skipping|"
    r"transcript\s+(?:assay|analysis)"
)
_PROTEIN_FUNCTION_RE = re.compile(
    r"functional|in[\s-]*vitro|in[\s-]*vivo|electrophysiolog|"
    r"patch[\s-]*clamp|biochemical|enzymatic|knock[\s-]*(?:in|out)|"
    r"expression\s+(?:assay|study|studies)|"
    r"transfect|western\s+blot|reporter\s+assay"
)


def _is_splicing_only_assay(text: str) -> bool:
    """True when the assay evidence is RNA/splicing and nothing else.

    Per Walker 2023 that is not PS3/BS3 evidence. Caller must pass a
    LOWERCASED string, matching the other gates in this module."""
    return bool(_SPLICING_ASSAY_RE.search(text)) and not _PROTEIN_FUNCTION_RE.search(text)


def _residue_label(ev: dict) -> str:
    """Human-readable residue label for the PM5/PS1 evidence strings.

    The ClinVar retrieval matches on the MANE Select residue number (that is
    the frame ClinVar's ``Name`` field is written in) but the curator is
    working in the numbering of the transcript they supplied. When the two
    differ, name BOTH — otherwise the evidence string cites a residue number
    that does not correspond to the ClinVar records quoted alongside it,
    which is exactly what hid the GATA4 cross-transcript bug.
    """
    pos = ev.get("protein_position")
    matched = ev.get("matched_protein_position")
    if matched is not None and matched != pos:
        return f"residue {pos} (= residue {matched} on ClinVar's MANE transcript)"
    return f"residue {pos}"


def _clinvar_pm5_criterion(pm5_ev: dict | None, gene: str | None) -> dict:
    """Deterministic PM5: a DIFFERENT missense change at the same residue is
    established Pathogenic/Likely-pathogenic in ClinVar at ≥2★ review.

    Previously AI-evaluated. The model honoured the "PM5 MUST-APPLY" prompt
    rule on only a minority of qualifying variants even though
    ``get_pm5_evidence`` had already retrieved the ≥2★
    same-residue candidates — so PM5 is now server-owned and enforced off that
    same retrieval. Strength Moderate (+2). The PS1 › PM5 and (Cardiomyopathy /
    RASopathy) PM1 › PM5 mutual exclusions still run post-merge in
    ``apply_cross_criterion_exclusions``.

    VCEP nuance, corrected 2026-09-01: strength is NOT a constant. The 8
    Cardiomyopathy genes put PM5 at Moderate only when the same-codon comparator
    is classified PATHOGENIC, and at Supporting when it is LIKELY pathogenic.
    And the RASopathy "≥2 distinct changes" requirement belongs to their STRONG
    row, not Moderate — gating Moderate on it was an under-call.
    """
    ev = pm5_ev or {}
    if not ev.get("ok") or ev.get("not_applicable"):
        return _hard_coded_entry(
            "PM5", "not_met", None,
            "PM5 not evaluated: non-missense, no protein residue, or the "
            "ClinVar same-residue lookup was unavailable.",
        )
    res = _residue_label(ev)
    gname = ev.get("gene") or gene
    two_star = [c for c in (ev.get("candidates") or []) if (c.get("stars") or 0) >= 2]
    n_nonmiss = ev.get("non_missense_excluded_two_star_count") or 0
    nonmiss_note = (
        f" ({n_nonmiss} ≥2★ P/LP record(s) at this residue are non-missense "
        "— nonsense/frameshift/synonymous — and are not PM5 comparators.)"
        if n_nonmiss else ""
    )
    if not two_star:
        n_any = ev.get("count", 0)
        note = (
            f"{n_any} P/LP missense record(s) at {res} of {gname} but "
            "none at ≥2★ review — below the PM5 confidence bar."
            if n_any else
            f"No other P/LP missense change at {res} of {gname}."
        )
        return _hard_coded_entry("PM5", "not_met", None, note + nonmiss_note)
    distinct = {
        m.group(1) for c in two_star
        if (m := _PM5_ALT_AA_RE.search(c.get("name") or ""))
    }

    _pm5_st = (_crit_spec(gname, "PM5").get("strengths") or {})
    _ties_to_comparator = (
        "classified as pathogenic" in (_pm5_st.get("Moderate") or "").lower()
        and "classified as likely pathogenic" in (_pm5_st.get("Supporting") or "").lower()
    )
    def _is_full_p(c: dict) -> bool:
        sig = (c.get("clinical_significance") or "").lower()
        return "pathogenic" in sig and not sig.strip().startswith("likely")
    _has_full_p = any(_is_full_p(c) for c in two_star)
    if _ties_to_comparator and not _has_full_p:
        _strength, _why = "PM5_Supporting", (
            f"{gname}'s VCEP allows PM5 at MODERATE only when the same-codon "
            f"comparator is classified PATHOGENIC; the ≥2★ comparator(s) here "
            f"are Likely pathogenic, which its spec puts at SUPPORTING"
        )
    else:
        _strength, _why = "PM5_Moderate", (
            "novel missense at a residue with established pathogenic missense"
        )
    _distinct_note = (
        f" ({len(distinct)} distinct alt-AA change(s) at the codon; the "
        f"RASopathy STRONG row needs ≥2 in ≥5 probands, a proband count "
        f"HeartVar does not hold, so Moderate is the ceiling here.)"
        if bool(_vcep_freq(gname)) and (_vcep_freq(gname) or {}).get("vcep") == "RASopathy"
        else ""
    )
    ex = two_star[0]
    return _hard_coded_entry(
        "PM5", "met", _strength,
        f"{len(two_star)} other Pathogenic/Likely-pathogenic missense change(s) "
        f"at {res} of {gname} in ClinVar at ≥2★ (e.g. {ex.get('name')} — "
        f"{ex.get('clinical_significance')}); {_why} — "
        f"{_strength.replace('_', ' ')}." + _distinct_note + nonmiss_note,
    )


def _clinvar_ps1_criterion(ps1_ev: dict | None, gene: str | None) -> dict:
    """Deterministic PS1: the SAME amino-acid change as the proband is an
    established Pathogenic/Likely-pathogenic variant in ClinVar at ≥2★ review,
    reached via a DIFFERENT nucleotide (ACMG-2015 PS1, Strong +4).

    Previously AI-evaluated, fired 0× — the same under-application that made
    PM5 server-owned — even though ``get_pm5_evidence`` already surfaces the
    same-AA candidates (``ps1_candidates``, the proband's own record removed)
    off the SAME ClinVar retrieval PM5 uses. PS1 is therefore now server-owned
    and enforced off that surface. Missense-only (the ClinVar surface is
    ``name LIKE '%p.%'`` + the p.XxxNNNYyy token); the Walker-2023 splice
    extension of PS1 is deliberately out of scope. The ≥2★ review bar for the
    "established" comparator mirrors the PM5 rule. Unlike PM5 there is NO
    codon-hotspot multiplicity requirement (PS1 is the SAME established
    pathogenic amino-acid change, not a residue-hotspot argument), so the
    RASopathy ≥2-distinct-change guard does not apply. The PS1 › PM5 mutual
    exclusion still runs post-merge in ``apply_cross_criterion_exclusions``.
    """
    ev = ps1_ev or {}
    if not ev.get("ok") or ev.get("not_applicable"):
        _sp = ev.get("splice_ps1") or {}
        _sp_cands = _sp.get("candidates") or []
        if _sp.get("ok") and _sp_cands:
            _names = ", ".join(str(c.get("name")) for c in _sp_cands[:2])
            _self = _sp.get("self_excluded") or []
            _self_note = (
                f" The proband's own ClinVar record "
                f"({', '.join(str(c.get('name')) for c in _self[:2])}) was "
                f"excluded — PS1 requires a different nucleotide change."
                if _self else ""
            )
            return _hard_coded_entry(
                "PS1", "met", "PS1",
                f"{len(_sp_cands)} established Pathogenic/Likely-pathogenic "
                f"ClinVar record(s) at >=2* review sit at the SAME canonical "
                f"splice position c.{_sp.get('position')} via a different "
                f"nucleotide change (e.g. {_names}). Two changes at one "
                f"canonical +/-1,2 position abolish the same splice site, so "
                f"their predicted RNA-splicing effects are equivalent — PS1 per "
                f"ClinGen SVI (Walker 2023).{_self_note}",
            )
        if _sp.get("ok") and not _sp.get("not_applicable"):
            return _hard_coded_entry(
                "PS1", "not_met", None,
                f"No established P/LP ClinVar record at >=2* review at the same "
                f"canonical splice position c.{_sp.get('position')} via a "
                f"different nucleotide change — PS1 not met on the Walker 2023 "
                f"splice route.",
            )
        return _hard_coded_entry(
            "PS1", "not_met", None,
            "PS1 not evaluated: non-missense, no protein residue, or the "
            "ClinVar same-amino-acid lookup was unavailable.",
        )
    res = _residue_label(ev)
    gname = ev.get("gene") or gene
    cands = [c for c in (ev.get("ps1_candidates") or []) if (c.get("stars") or 0) >= 2]
    self_ex = ev.get("ps1_self_excluded") or []
    self_note = ""
    if self_ex:
        _names = ", ".join(str(c.get("name")) for c in self_ex[:2])
        self_note = (
            f" The proband's OWN ClinVar record ({_names}) was excluded from "
            "this comparison — PS1 requires a different nucleotide, and a "
            "record naming the same variant on another transcript is not one."
        )
    if not cands:
        n_any = ev.get("ps1_count", 0)
        note = (
            f"{n_any} same-amino-acid P/LP record(s) at {res} of {gname} "
            "but none at ≥2★ review — below the PS1 confidence bar."
            if n_any else
            f"No same-amino-acid P/LP record (different nucleotide) at "
            f"{res} of {gname}."
        )
        return _hard_coded_entry("PS1", "not_met", None, note + self_note)
    ex = cands[0]
    return _hard_coded_entry(
        "PS1", "met", "PS1_Strong",
        f"ClinVar ≥2★ Pathogenic/Likely-pathogenic record with the SAME "
        f"amino-acid change at {res} of {gname} via a different "
        f"nucleotide (e.g. {ex.get('name')} — {ex.get('clinical_significance')}); "
        "ACMG-2015 PS1 Strong." + self_note,
    )


def _gnomad_popmax_af(ev: dict) -> float | None:
    """Pull the gnomAD popmax FAF95 (preferred) or plain AF (fallback)
    from the nested gnomAD evidence shape. Returns None when gnomAD has
    no record for the variant — the caller distinguishes "absent from
    gnomAD" (AF treated as 0) from "lookup unavailable" via gnomad.ok."""
    gn = ev.get("gnomad") or {}
    if not gn.get("ok"):
        return None
    variant = gn.get("variant") or None
    if not variant:
        return None
    ex = variant.get("exome") or {}
    ge = variant.get("genome") or {}
    ex_faf95 = ex.get("faf95") or {}
    ge_faf95 = ge.get("faf95") or {}
    popmax_vals = [
        v for v in (ex_faf95.get("popmax"), ge_faf95.get("popmax"))
        if v is not None
    ]
    if popmax_vals:
        return float(max(popmax_vals))
    af_vals = [v for v in (ex.get("af"), ge.get("af")) if v is not None]
    return float(max(af_vals)) if af_vals else None


def _int_or_none(v) -> int | None:
    """Coerce a gnomAD allele-number field to int, or None when unusable.

    None is what the BA1 allele-number gate reads as "cannot verify", and it
    fails closed on that, so returning None must stay the only failure mode."""
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _gnomad_point_af(
    ev: dict, *, allow: list[str] | None = None,
    exclude: list[str] | None = None,
    min_alleles: int | None = None,
) -> tuple[float | None, str, int | None]:
    """Highest PER-POPULATION point allele frequency, per a VCEP's own rules.

    Returns ``(af, basis, an)`` where basis names the population the value came
    from, so the report can say what it actually compared rather than claiming
    FAF95 over a different statistic, and ``an`` is that SAME population's
    allele number.

    ``an`` exists because BA1's ">= 2,000 observed alleles" gate has to be
    checked against the dataset the AF actually came from. It used to read the
    block-level AN via ``_gnomad_popmax_an``, which a point_af fixture does not
    carry — so the gate saw None, failed closed, and BA1 could NEVER fire on a
    point_af gene. Found 2026-09-01: a 5% variant in 100,000 alleles was
    withheld on FBN1, KCNQ1 and BMPR2 while firing correctly on MYH7. BS1 defers
    when BA1 is withheld, so those variants got no benign criterion at all —
    a false negative in the toward-pathogenic direction, live since cc29ab6
    (2026-08-26) added the point_af routing without updating the AN gate.

    EXOME-PREFERRED, and that choice matters. Max-of-blocks — taking each of the
    exome/genome blocks' popmax and maxing them — is monotonically upward-biased
    and demonstrably wrong: on MYH7 c.4130C>T it invents a Finnish AF of 9.4e-5
    for a variant with ac=0 in 53,420 Finnish exomes, and that bias is exactly
    what produced the reassuring "every mover is benign" reading in the earlier
    measurement. Summing the blocks instead dilutes exome-ascertained variants.
    Exome-preferred is the most like-for-like with the v2.1.1-era calibration
    these thresholds were derived under.

    ``allow`` is an allow-list (KCNQ1 GN112 names five continental populations
    and nothing else counts); ``exclude`` is a deny-list (FBN1 GN022 bars
    Finnish, Ashkenazi Jewish and "Other" — gnomAD v4 renamed ``oth`` to
    ``remaining``, so both spellings belong in the list); ``min_alleles`` drops a
    population whose AN is below the spec's floor, because an AF computed off a
    handful of alleles is noise and these rules fire in the BENIGN direction.
    """
    gn = ev.get("gnomad") or {}
    if not gn.get("ok"):
        return None, "unavailable", None
    variant = gn.get("variant") or None
    if not variant:
        return None, "absent from gnomAD", None
    block = variant.get("exome") or variant.get("genome") or {}
    which = "exome" if variant.get("exome") else "genome"
    pops = block.get("populations") or []
    allow_set = {p.lower() for p in (allow or [])}
    deny_set = {p.lower() for p in (exclude or [])}
    best_af: float | None = None
    best_pop = ""
    best_an: int | None = None
    skipped_small: list[str] = []
    for pop in pops:
        pid = str(pop.get("id") or "").lower()
        if not pid:
            continue
        if allow_set and pid not in allow_set:
            continue
        if pid in deny_set:
            continue
        try:
            ac, an = int(pop.get("ac") or 0), int(pop.get("an") or 0)
        except (TypeError, ValueError):
            continue
        if an <= 0:
            continue
        if min_alleles is not None and an < min_alleles:
            skipped_small.append(pid)
            continue
        af = ac / an
        if best_af is None or af > best_af:
            best_af, best_pop, best_an = af, pid, an
    if best_af is None:
        if not pops:
            overall = block.get("af")
            if overall is not None:
                return (float(overall),
                        f"{which} overall AF (no per-population data)",
                        _int_or_none(block.get("an")))
            return None, "absent from gnomAD", None
        return 0.0, f"{which} point AF, 0 in every eligible population", None
    note = f"{which} point AF in {best_pop}"
    if skipped_small:
        note += f" (excluded for AN<{min_alleles}: {', '.join(sorted(skipped_small))})"
    return best_af, note, best_an


def _vcep_freq_af(
    ev: dict, gene: str | None,
) -> tuple[float | None, str, int | None]:
    """The allele frequency to compare against this gene's VCEP thresholds.

    Returns ``(af, basis, an)``. The AN is the allele number of whatever dataset
    supplied the AF — the winning POPULATION for a point_af gene, the winning
    exome/genome BLOCK for a FAF95 gene — because BA1's allele-number gate has to
    be checked against the same dataset. This is the single routing point that
    knows which metric is in play, so it is the only place that can answer it.

    PER-VCEP ROUTING, which is the conformant change and not the same as a
    blanket metric switch. The Cardiomyopathy and RASopathy VCEPs specify the
    filtering allele frequency, so FAF95 popmax is ALREADY right for 23 of the 26
    spec'd genes — swapping the metric globally would break the 619 of 852
    records it is right for in order to fix 118. KCNQ1 (GN112), FBN1 (GN022) and
    BMPR2 (GN125) specify a per-population point AF; a scan of every spec'd
    gene's BA1/BS1/PM2 rule text in the registry confirms they are the only
    three. (Read "the only two" before BMPR2's rule set was harvested on
    2026-09-01 — the claim was true of the 25 genes then loaded.)
    """
    spec = _vcep_freq(gene) or {}
    if spec.get("metric") == "point_af":
        return _gnomad_point_af(
            ev,
            allow=spec.get("populations"),
            exclude=spec.get("exclude_populations"),
            min_alleles=spec.get("min_alleles"),
        )
    af = _gnomad_popmax_af(ev)
    return af, "popmax FAF95", _gnomad_popmax_an(ev)


def _gnomad_popmax_an(ev: dict) -> int | None:
    """Observed allele number (AN) for the gnomAD block that supplied the
    frequency. Mirrors ``_gnomad_popmax_af``'s exome/genome selection so BA1's
    ">= 2,000 observed alleles" floor is checked against the dataset the AF
    actually came from. Returns None when AN is absent (caller fails closed)."""
    gn = ev.get("gnomad") or {}
    if not gn.get("ok"):
        return None
    variant = gn.get("variant") or None
    if not variant:
        return None
    ex = variant.get("exome") or {}
    ge = variant.get("genome") or {}
    ex_faf = (ex.get("faf95") or {}).get("popmax")
    ge_faf = (ge.get("faf95") or {}).get("popmax")
    if ex_faf is not None or ge_faf is not None:
        chosen = ex if (ge_faf is None or (ex_faf is not None and ex_faf >= ge_faf)) else ge
    elif ex.get("af") is not None or ge.get("af") is not None:
        ex_af, ge_af = ex.get("af"), ge.get("af")
        chosen = ex if (ge_af is None or (ex_af is not None and ex_af >= ge_af)) else ge
    else:
        chosen = ex or ge
    an = chosen.get("an")
    try:
        return int(an)
    except (TypeError, ValueError):
        return None


_BA1_EXCEPTION_VARIANTS: dict[str, frozenset[str]] = {
    "ACAD9": frozenset({"c.-44_-41duptaag"}),
    "GJB2": frozenset({"c.109g>a"}),
    "HFE": frozenset({"c.187c>g",
                      "c.845g>a"}),
    "MEFV": frozenset({"c.1105c>t",
                       "c.1223g>a"}),
    "PIBF1": frozenset({"c.1214g>a"}),
    "ACADS": frozenset({"c.511c>t"}),
    "BTD": frozenset({"c.1330g>c"}),
}


def _ba1_exception_hit(gene: str | None, vep: dict) -> str | None:
    """The normalised ``c.`` token when this variant is on the SVI BA1 exception
    list, else None. Matching needs an HGVS c.; when VEP supplies none we cannot
    identify the variant and BA1 proceeds normally (documented limitation — the
    list is variant-specific, and suppressing BA1 for the whole gene would be
    far more wrong than missing an unidentifiable variant)."""
    listed = _BA1_EXCEPTION_VARIANTS.get((gene or "").strip().upper())
    if not listed:
        return None
    hgvsc = (vep or {}).get("hgvsc") or ""
    if not isinstance(hgvsc, str) or not hgvsc:
        return None
    token = hgvsc.rpartition(":")[2].strip().lower() or hgvsc.strip().lower()
    return token if token in listed else None


def _gnomad_hom_hemi(ev: dict) -> tuple[int, int]:
    """Return (max hom count, max hemi count) across exome and genome.
    Both default to 0 when gnomAD has no record."""
    variant = (ev.get("gnomad") or {}).get("variant") or {}
    if not variant:
        return 0, 0
    ex = variant.get("exome") or {}
    ge = variant.get("genome") or {}
    hom = max(int(ex.get("ac_hom") or 0), int(ge.get("ac_hom") or 0))
    hemi_vals = [
        int(v) for v in (ex.get("ac_hemi"), ge.get("ac_hemi"))
        if v is not None
    ]
    return hom, max(hemi_vals) if hemi_vals else 0


def _spliceai_per_score_max(ev: dict) -> float:
    """Max SpliceAI delta across DS_AG/DS_AL/DS_DG/DS_DL on the top
    transcript, plus the precomputed top-level ``max_delta`` as a
    fallback. Returns 0.0 when SpliceAI is unavailable / not applicable."""
    sa = ev.get("spliceai") or {}
    if not sa.get("ok"):
        return 0.0
    candidates: list[float] = []
    md = sa.get("max_delta")
    if md is not None:
        try:
            candidates.append(float(md))
        except (TypeError, ValueError):
            pass
    for t in (sa.get("scores_per_transcript") or [])[:1]:
        for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL"):
            v = t.get(k)
            if v is None:
                continue
            try:
                candidates.append(float(v))
            except (TypeError, ValueError):
                pass
    return max(candidates) if candidates else 0.0


_VCEP_FREQ_FILE = Path(__file__).resolve().parent.parent / "data" / "vcep_frequency_thresholds.json"
try:
    _VCEP_FREQ = json.loads(_VCEP_FREQ_FILE.read_text()).get("genes") or {}
except (OSError, ValueError) as _e:  # pragma: no cover - defensive
    log.warning("VCEP frequency table unavailable (%s) — using generic thresholds", _e)
    _VCEP_FREQ = {}


def _vcep_freq(gene: str | None) -> dict | None:
    """VCEP frequency-threshold dict for `gene`, or None when not covered."""
    return _VCEP_FREQ.get((gene or "").upper())


_VCEP_CRIT_FILE = Path(__file__).resolve().parent.parent / "data" / "vcep_criteria_spec.json"
try:
    _VCEP_CRIT = json.loads(_VCEP_CRIT_FILE.read_text()).get("genes") or {}
except (OSError, ValueError) as _e:  # pragma: no cover - defensive
    log.warning("VCEP criteria-specification table unavailable (%s)", _e)
    _VCEP_CRIT = {}


def _crit_spec(gene: str | None, code: str) -> dict:
    """The registry entry for one gene/criterion, or {} when absent."""
    return ((_VCEP_CRIT.get((gene or "").upper()) or {}).get(code) or {})


_PS3_RANK = {
    "PS3": 4, "PS3_Strong": 4, "PS3_Moderate": 2, "PS3_Supporting": 1,
}


_PS3_MODERATE_MIN_CONTROLS = 11


_PS4_OCCURRENCE_RE = re.compile(
    r"PS4(?:_(Supporting|Moderate|Strong))?\s*:\s*"
    r"(?:>=|≥|at least|more than)?\s*(\d+)\s*independent\s+occurrence",
    re.IGNORECASE,
)


_PS4_OR_ROUTE_RE = re.compile(
    r"(?i:odds\s*ratio)|\bOR\b|(?i:95%\s*(?:confidence\s*interval|CI))")

_PS4_ROW_COUNT_RE = re.compile(
    r"(?P<op>>=|>|≥|at least\s|Met by\s)?\s*(?P<n>\d+)"
    r"(?:\s*(?:-|to|–)\s*\d+)?"
    r"\s*(?:or\s+more\s+)?"
    r"(?:(?:unrelated|independent|affected|other|index)\s+)*"
    r"(?P<unit>probands?|patients?|individuals?|observations?"
    r"|occurrences?|families|kindreds?)",
    re.IGNORECASE,
)

_PS4_LADDER_LABELS = ("Supporting", "Moderate", "Strong")


def _ps4_min_count_from_row(text: str) -> int | None:
    """Minimum person-count a single PS4 strength row demands, or None.

    None means "this row publishes no count route" — either it is the
    odds-ratio route, or it is a points ladder (FBN1: "If >= 4 points."), or it
    says nothing countable. Never raises.
    """
    t = (text or "").strip()
    if not t or _PS4_OR_ROUTE_RE.search(t):
        return None
    m = _PS4_ROW_COUNT_RE.search(t)
    if not m:
        return None
    n = int(m.group("n"))
    if (m.group("op") or "").strip() == ">":
        n += 1
    return n if n > 0 else None


def _ps4_spec_occurrence_ladder(gene: str | None) -> dict[str, int]:
    """{strength: minimum independent occurrences} from this gene's own PS4
    record, or {} when it publishes no count route.

    Two sources, in order. The RASopathy panel puts the whole ladder in
    ``comments`` ("PS4: >=5 independent occurrences / PS4_Moderate: >=3 /
    PS4_Supporting: >=1"). KCNQ1 and BMPR2 instead put one threshold in each
    STRENGTH ROW, which the comments-only reader missed entirely — both fell
    through to the generic Kelly ladder and were graded on thresholds their
    VCEP does not publish.

    {} is returned for the eight Cardiomyopathy genes (odds-ratio only, see
    ``_ps4_case_control_only``) and for FBN1 (a points ladder with no published
    points-per-proband mapping). Never raises.
    """
    spec = _crit_spec(gene, "PS4")
    out: dict[str, int] = {}
    for m in _PS4_OCCURRENCE_RE.finditer(spec.get("comments") or ""):
        label = (m.group(1) or "Strong").capitalize()
        n = int(m.group(2))
        if label not in out or n < out[label]:
            out[label] = n
    if out:
        return out
    rows = spec.get("strengths") or {}
    if not isinstance(rows, dict):
        return {}
    for label in _PS4_LADDER_LABELS:
        n = _ps4_min_count_from_row(rows.get(label) or "")
        if n is not None:
            out[label] = n
    if len(out) >= 2:
        return out
    return _ps4_points_ladder(gene, rows)


_PS4_POINTS_RE = re.compile(r"(?:>=|≥|at least)\s*([\d.]+)\s*points?", re.IGNORECASE)


def _ps4_points_ladder(gene: str | None, rows: dict) -> dict[str, int]:
    """{strength: minimum cases} from a POINTS ladder, or {}.

    Only for VCEPs that publish a points-per-case value. Never raises."""
    freq = _vcep_freq(gene) or {}
    if freq.get("vcep") != "RASopathy":
        return {}
    out: dict[str, int] = {}
    for label in _PS4_LADDER_LABELS:
        m = _PS4_POINTS_RE.search(rows.get(label) or "")
        if not m:
            continue
        try:
            pts = float(m.group(1))
        except ValueError:
            continue
        n = int(pts) if pts == int(pts) else int(pts) + 1
        if n > 0:
            out[label] = n
    return out if len(out) >= 2 else {}


def _ps4_case_control_only(gene: str | None) -> bool:
    """True when this gene's PS4 accepts ONLY a case-control odds ratio.

    The eight Cardiomyopathy genes publish an OR bound at every rung and no
    proband route at any of them, so a proband count earns PS4 at no strength
    however large it is. An odds ratio needs a CASE cohort, which HeartVar does
    not have; the honest answer for those genes is that PS4 is unavailable
    rather than gradeable on probands. Never raises.
    """
    spec = _crit_spec(gene, "PS4")
    rows = spec.get("strengths") or {}
    if not isinstance(rows, dict) or not rows:
        return False
    if _ps4_spec_occurrence_ladder(gene):
        return False
    texts = [rows.get(lab) or "" for lab in _PS4_LADDER_LABELS]
    present = [t for t in texts if t.strip()]
    if not present:
        return False
    return all(_PS4_OR_ROUTE_RE.search(t) for t in present)


def _ps3_validated_strength(facts: dict | None) -> tuple[str | None, str]:
    """Strength EARNED by an assay's demonstrated validation, per Brnich 2019.

    Returns ``(strength_label_or_None, why)``. None means "no evidence
    demonstrated", which is the framework's starting point and NOT an error — it
    is the expected answer for a publication that does not describe its controls.
    """
    f = facts if isinstance(facts, dict) else {}
    n = f.get("assay_variant_controls")
    n = n if isinstance(n, int) and not isinstance(n, bool) and n >= 0 else None
    lab = f.get("assay_lab_controls") is True

    if n is not None and n >= _PS3_MODERATE_MIN_CONTROLS:
        return "Moderate", (
            f"the assay was validated against {n} previously classified "
            f"pathogenic/benign variant controls, meeting the "
            f"{_PS3_MODERATE_MIN_CONTROLS}-control bar ClinGen SVI (Brnich 2019) "
            f"sets for moderate-level functional evidence"
        )
    if lab or (n is not None and n > 0):
        detail = (
            f"{n} classified variant control(s), short of the "
            f"{_PS3_MODERATE_MIN_CONTROLS} needed for Moderate"
            if n else "appropriate laboratory controls"
        )
        return "Supporting", (
            f"the assay reports {detail}, which ClinGen SVI (Brnich 2019) allows "
            f"at supporting level for work performed rigorously with appropriate "
            f"laboratory controls"
        )
    return None, (
        "the assay's validation is not described — no classified variant "
        "controls and no laboratory controls were reported. ClinGen SVI "
        "(Brnich 2019) directs that functional evidence 'start from the "
        "assumption of no evidence', with strength earned from demonstrated "
        "validation, so no PS3 weight is applied"
    )


def _ps3_spec_ceiling(gene: str | None) -> str | None:
    """The highest PS3 strength this gene's VCEP publishes, or None.

    "Strong" for 11 of the 26 spec'd genes, "Moderate" for 14 (all RASopathy).
    None when the gene has no harvested PS3 rows — unspec'd, or spec'd but
    omitting PS3 (SHOC2) — meaning there is no VCEP ceiling to impose. It does
    NOT mean "unlimited": for those genes the strength comes from the assay's own
    demonstrated validation via ``_ps3_validated_strength``, which starts at no
    evidence. The two act together, whichever is lower.
    """
    strengths = (_crit_spec(gene, "PS3").get("strengths") or {})
    for label in ("Strong", "Moderate", "Supporting"):
        if label in strengths:
            return label
    return None


_PS3_NO_CONTROLS_FLOOR_RE = re.compile(
    r"if no known variant validation controls"
    r".{0,120}?"
    r"then score at the (supporting|moderate|strong) strength",
    re.IGNORECASE | re.DOTALL,
)


def _ps3_spec_floor(gene: str | None) -> str | None:
    """The strength this gene's VCEP publishes for an assay with NO validation
    controls, or None when it publishes no such floor.

    None is the common case and means "defer to Brnich 2019", i.e. start from
    the assumption of no evidence. It does NOT mean the gene is unspec'd.
    """
    for text in (_crit_spec(gene, "PS3").get("strengths") or {}).values():
        m = _PS3_NO_CONTROLS_FLOOR_RE.search(text or "")
        if m:
            return m.group(1).capitalize()
    return None


def _crit_text(gene: str | None, code: str) -> str:
    """All published rule text for one gene/criterion, strengths joined.

    The registry splits a criterion's text across its strength rows, and rules
    that matter to us are spread over them: MYH7's PM1 codon range sits in the
    Moderate row while FBN1's cbEGF rule is Strong and its EGF/TB rules Moderate.
    Callers that scan for hotspot ranges or a caveat want the union, so joining is
    the correct default rather than a convenience.

    Falls back to the flat ``text`` key so the superseded hand table still reads
    if it is ever pointed at again. Verified equivalent: parsing PM1 hotspots from
    the joined strengths reproduces the old table's ranges and exons for all 25
    genes, with zero differences.
    """
    entry = _crit_spec(gene, code)
    strengths = entry.get("strengths")
    if isinstance(strengths, dict) and strengths:
        return " ".join(str(v) for v in strengths.values() if v)
    return entry.get("text") or ""


def _criterion_applicable(gene: str | None, code: str) -> bool:
    """False ONLY when the gene's VCEP explicitly marks `code` Not Applicable.
    Unknown / uncovered genes default to True (never suppress without evidence)."""
    spec = _VCEP_CRIT.get((gene or "").upper()) or {}
    return (spec.get(code) or {}).get("applicability") != "not_applicable"


_GENE_DOSAGE_FILE = Path(__file__).resolve().parent.parent / "data" / "gene_mechanism.json"
try:
    _GENE_DOSAGE = json.loads(_GENE_DOSAGE_FILE.read_text()).get("genes") or {}
except (OSError, ValueError) as _e:  # pragma: no cover - defensive
    log.warning("ClinGen gene-dosage table unavailable (%s)", _e)
    _GENE_DOSAGE = {}

# CC0, license-clean mode-of-inheritance + evidence-tier layer used by the
_CLINGEN_GV_FILE = Path(
    os.environ.get("CLINGEN_GV_PATH")
    or (Path(__file__).resolve().parent.parent / "data" / "clingen_gene_validity.json")
)
try:
    _CLINGEN_GV_RAW = json.loads(_CLINGEN_GV_FILE.read_text())
    _CLINGEN_GV = _CLINGEN_GV_RAW.get("genes") or {}
    _CLINGEN_GV_META = _CLINGEN_GV_RAW.get("meta") or {}
except (OSError, ValueError) as _e:  # pragma: no cover - defensive
    log.warning("ClinGen Gene-Disease Validity table unavailable (%s)", _e)
    _CLINGEN_GV = {}
    _CLINGEN_GV_META = {}


_GV_TIER_ORDER = (
    "Definitive", "Strong", "Moderate", "Limited", "Disputed", "Refuted",
    "No Known Disease Relationship", "Animal Model Only",
)
_GV_RANK = {t: i for i, t in enumerate(_GV_TIER_ORDER)}
_GV_CEILING_FROM = _GV_RANK["Moderate"]


def gene_validity_ceiling(gene: str | None) -> str | None:
    """The ClinGen Gene-Disease Validity tier that governs this gene, or None.

    THE RULE, in the RASopathy VCEP's own words on MRAS c.67G>C (verbatim from
    that variant's eRepo record):

        "Given the Moderate strength of gene-disease relationship between MRAS
         and autosomal dominant RASopathy, ClinGen's sequence variant
         interpretation working group does not recommend the classification of
         variants in this gene beyond likely pathogenic."

    HeartVar has shipped this table since 2026-08 and already reads its
    ``classification`` field — but only as an adequacy gate for
    mode-of-inheritance (``_ADEQUATE_TIERS`` in ``gene_inheritance_modes``,
    which admits Moderate). The tier was never applied as a ceiling on the
    variant classification, so an MRAS variant could score Pathogenic where the
    VCEP says the gene does not support a call beyond Likely pathogenic.

    ⚠ THE HIGHEST assertion governs, deliberately, and the disease is NOT
    matched. A gene can hold several assertions at different tiers — MYH7 is
    Limited for ARVC and for CHD but Definitive for HCM and DCM. Selecting by
    fuzzy phenotype text would let a mis-parsed indication cap a gene ClinGen
    calls Definitive, which demotes true Pathogenic calls: the dangerous
    direction. Taking the maximum means the cap fires only when ClinGen asserts
    NO strong relationship for the gene at all, which is the case the SVI
    sentence above actually describes.

    ⚠ SCOPED TO THE GENES THAT HAVE A VCEP CRITERIA SPEC. Unscoped, the rule
    capped 85 of the 476 genes in clingen_gene_validity.json, 84 of them outside
    the 27 spec'd genes, and it did so on whatever assertion ClinGen happened to
    have curated — which for some genes is a disease no cardiac curator is
    asking about. ATP2A2 holds exactly one assertion, Refuted for EPILEPSY, so
    an ATP2A2 cardiac call was being held at Likely pathogenic on the strength
    of a refuted epilepsy link. PIK3CA (Refuted, hereditary breast carcinoma)
    and COL11A1 (Moderate, hearing loss) are the same shape. ClinGen having
    never curated a gene for a cardiac indication is not evidence that its
    cardiac relationship is weak.

    Most of the other 84 ARE cardiac assertions (Congenital Heart Disease 25,
    Dilated Cardiomyopathy 9, Hereditary Cardiovascular Disease 7), so whether
    a Moderate CHD assertion ought to cap a cardiac call is a live question —
    but it is a clinical-scope judgement, and answering it by demoting 84 genes
    is not a change to make without review. Inside the spec'd set the rule is
    unambiguous: those genes have a published VCEP specification, MRAS is the
    only one it caps, and the RASopathy VCEP states the cap in the MRAS record
    itself. Everything else fails open, which is the safe direction.

    Returns None for an unknown gene, a gene with no VCEP spec, or a missing
    table, so the ceiling fails open and can never invent a demotion. Never
    raises.
    """
    g = (gene or "").strip().upper()
    if not g or g not in _VCEP_CRIT:
        return None
    best, best_rank = None, len(_GV_TIER_ORDER)
    for a in (_CLINGEN_GV.get(g) or []):
        tier = (a.get("classification") or "").strip()
        rank = _GV_RANK.get(tier)
        if rank is not None and rank < best_rank:
            best, best_rank = tier, rank
    return best


def gene_validity_caps_pathogenic(gene: str | None) -> bool:
    """Whether this gene's validity tier caps its variants at Likely pathogenic."""
    tier = gene_validity_ceiling(gene)
    return tier is not None and _GV_RANK.get(tier, -1) >= _GV_CEILING_FROM


def _as_float(v) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


_APPLICABILITY_GATED_CODES = (
    "PP2", "BP1", "PP4", "PM1", "PVS1",
    "PM3", "PM6", "PS3",
    "BP3", "BP5", "BP6", "BS2", "BS3", "PP5",
    "BP2",
)

_PP4_EVIDENCE_ALLOW_LIST = frozenset({"FBN1"})

_PM1_LOF_CONSEQUENCES = frozenset({
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "transcript_ablation",
})


def _parse_pm1_hotspots(pm1_text: str) -> tuple[list[tuple[int, int]], set[int]]:
    """Parse a VCEP PM1.text hotspot list into (codon_ranges, exon_numbers).

    The text is a curated hotspot enumeration, so any ``N-N`` span is a codon/AA
    range — captures EVERY range, fixing the old ``(?:AA|Codons?)``-anchored
    regex that dropped 2nd+ ranges (e.g. MYBPC3 "Codons 485-502 and 1248-1266"
    kept only 485-502). ``exon N`` tokens are exon-defined hotspots (e.g. BRAF
    "exon 6, exon 11") that the residue-range matcher can't resolve.
    """
    text = pm1_text or ""
    ranges = [
        (int(a), int(b))
        for a, b in re.findall(r"(\d+)\s*[-–]\s*(\d+)", text)
    ]
    ranges += [
        (int(a), int(b))
        for a, b in re.findall(
            r"(?:amino\s*acids?|residues?|codons?|positions?)\s*"
            r"(\d+)\s*(?:to|through)\s*(\d+)", text, flags=re.I)
    ]
    # spec-licensed PM1. That is why PTPN11 NM_002834.4:c.781C>T (p.Leu261Phe),
    _consumed: set[int] = set()
    for _pat in (r"(\d+)\s*[-–]\s*(\d+)",
                 r"(?:amino\s*acids?|residues?|codons?|positions?)\s*"
                 r"(\d+)\s*(?:to|through)\s*(\d+)"):
        for _m in re.finditer(_pat, text, flags=re.I):
            _consumed.update(range(_m.start(), _m.end()))
    singles = sorted({
        int(m.group(1))
        for m in re.finditer(
            r"(?:AA|amino\s*acids?|residues?|codons?|positions?)\s*(\d+)"
            r"(?!\s*(?:[-–]|to\b|through\b)\s*\d)", text, flags=re.I)
        if not (set(range(m.start(), m.end())) & _consumed)
    })
    ranges += [(p, p) for p in singles]
    exons = {int(m) for m in re.findall(r"exon\s*(\d+)", text, flags=re.I)}
    return list(dict.fromkeys(ranges)), exons


_FBN1_CB_EGF = "cbEGF"


def _fbn1_domain_class(name: str | None) -> str | None:
    """Classify one UniProt FBN1 domain description into the spec's vocabulary."""
    d = (name or "").strip().lower()
    if not d:
        return None
    if "egf-like" in d:
        return _FBN1_CB_EGF if "calcium-binding" in d else "EGF-like"
    if d.startswith("tb"):
        return "TB"
    if "hybrid" in d:
        return "hybrid"
    return None


def _cys_index_within(sequence: str | None, start: int, end: int,
                      pos: int) -> int | None:
    """Which inter-cysteine interval of a domain ``pos`` falls in.

    Returns n when ``pos`` lies strictly between the nth and (n+1)th cysteine of
    the domain, so n == 2 means "between Cys2 and Cys3" — the interval FBN1's
    critical-Gly rule names. None when the sequence is unavailable, the domain
    holds too few cysteines, or the residue is not between two of them.
    """
    if not sequence or not (1 <= start <= end <= len(sequence)):
        return None
    cys = [i for i in range(start, end + 1) if sequence[i - 1] == "C"]
    for n in range(1, len(cys)):
        if cys[n - 1] < pos < cys[n]:
            return n
    return None


def _pm1_fbn1(evidence: dict | None, hgvsp: str | None) -> tuple[str | None, str]:
    """Evaluate FBN1 GN022's PM1 enumeration.

    Returns ``(verdict, reason)`` where verdict is "Strong", "Moderate",
    "barred" (a caveat forbids PM1), "no_rule" (enumeration does not cover this
    variant) or None (cannot decide — leave the incoming PM1 untouched).
    """
    ref, alt = _hgvsp_ref_alt_aa(hgvsp)
    pos = _hgvsp_int_position(hgvsp or "")
    if pos is None or not ref:
        return None, "residue or reference amino acid unresolved"

    up = (evidence or {}).get("uniprot") or {}
    domains = up.get("domains") or []
    sequence = up.get("sequence")
    containing = [
        d for d in domains
        if isinstance(d.get("start"), int) and isinstance(d.get("end"), int)
        and d["start"] <= pos <= d["end"]
    ]
    classes = {_fbn1_domain_class(d.get("name")) for d in containing}
    classes.discard(None)
    in_cb = _FBN1_CB_EGF in classes

    if ref == "G" and alt == "A":
        return "barred", ("GN022 records G->A as possibly tolerated and bars PM1 "
                          "for it")
    if ref == "N" and alt == "S" and in_cb:
        return "barred", ("GN022 records N->S in the consensus sequence as "
                          "possibly tolerated and bars PM1 for it (applied to any "
                          "N->S in a cbEGF-like domain, since the consensus "
                          "positions are not enumerated in the spec)")

    if ref == "C" and in_cb:
        names = ", ".join(sorted(d.get("name") or "?" for d in containing))
        return "Strong", f"cysteine residue in a cbEGF-like domain ({names})"

    if ref == "C" and classes & {"EGF-like", "TB", "hybrid"}:
        names = ", ".join(sorted(d.get("name") or "?" for d in containing))
        return "Moderate", f"cysteine residue in {names}"

    if alt == "C":
        return "Moderate", "cysteine-creating variant"

    if ref == "G" and in_cb:
        cb = [d for d in containing
              if _fbn1_domain_class(d.get("name")) == _FBN1_CB_EGF]
        for d in cb:
            n = _cys_index_within(sequence, d["start"], d["end"], pos)
            if n == 2:
                return "Moderate", (f"critical glycine between Cys2 and Cys3 of "
                                    f"{d.get('name')}")
            if n == 3 and any(
                isinstance(o.get("start"), int) and o["start"] < d["start"]
                and _fbn1_domain_class(o.get("name")) == _FBN1_CB_EGF
                for o in domains
            ):
                return "Moderate", (f"glycine between Cys3 and Cys4 of "
                                    f"{d.get('name')}, with an upstream cbEGF "
                                    f"domain")

    if in_cb:
        return None, "inside a cbEGF-like domain; no implemented rule decides it"
    if classes:
        return "no_rule", (f"in {', '.join(sorted(classes))} but not a cysteine; "
                           f"GN022 scopes PM1 there to cysteine residues")
    if not domains:
        return None, "no UniProt domain annotation available"
    return "no_rule", "outside every annotated FBN1 domain"


def _vep_exon_number(vep: dict | None) -> int | None:
    """Leading exon number from VEP's ``exon`` "N/total" string ("6/18" -> 6).
    None for intronic / no-exon / boundary-spanning ("35-36/40") fields."""
    exon = (vep or {}).get("exon")
    if not exon or "/" not in str(exon):
        return None
    head = str(exon).split("/", 1)[0].strip()
    return int(head) if head.isdigit() else None


def _load_hpo_json(name: str) -> dict:
    try:
        return json.loads(
            (Path(__file__).resolve().parent.parent / "data" / name).read_text()
        )
    except (OSError, ValueError) as _e:  # pragma: no cover - defensive
        log.warning("PM1 phenotype gate: %s unavailable (%s)", name, _e)
        return {}


_HPO_CATEGORY_SEEDS = _load_hpo_json("panelapp_hpo_map.json")
_HCM_HPO_SET = frozenset(
    str(t).upper() for t in (_HPO_CATEGORY_SEEDS.get("hcm") or ())
)
_DCM_HPO_SET = frozenset(
    str(t).upper() for t in (_HPO_CATEGORY_SEEDS.get("dcm") or ())
)


_HCM_TEXT_RE = re.compile(r"\bhcm\b|hypertroph\w*\s+cardiomyopath", re.I)
_DCM_TEXT_RE = re.compile(r"\bdcm\b|dilated\s+cardiomyopath", re.I)


def _phenotype_tokens(hpo: str | list[str] | None) -> tuple[set[str], str]:
    """Normalise the phenotype field into (upper-cased tokens, raw text).

    ⚠ THE FIELD IS FREE TEXT. ``CurationRequest.hpo`` is a ``str``, and app.py
    passes ``req.hpo`` straight into the gate — so a helper annotated
    ``list[str]`` and iterating its argument gets CHARACTERS, not terms.
    "HP:0001644" became {'H','P',':','0','1','6','4'}, no HPO ID ever matched,
    and the PM1 phenotype caveat was INERT in production while every unit test
    passed, because the tests all handed it a list. Found 2026-08-26 by checking
    what the request actually carries rather than what the signature claimed.

    So both shapes are accepted, and both HPO IDs and free-text disease names are
    matched: curators type "dilated cardiomyopathy" as readily as HP:0001644, and
    the field has always allowed it.
    """
    if hpo is None:
        return set(), ""
    if isinstance(hpo, str):
        text = hpo
        parts = re.split(r"[,;\n|]+", hpo)
    else:
        parts = [str(x) for x in hpo]
        text = ", ".join(parts)
    return {p.strip().upper() for p in parts if p.strip()}, text


def _pm1_non_hcm_phenotype(hpo: str | list[str] | None) -> str | None:
    """Reason string when the proband's phenotype points at a cardiomyopathy
    OTHER than HCM, else None.

    The Cardiomyopathy VCEP derived its PM1 hotspot ranges from **HCM case
    cohorts** and says so explicitly: "this rule should NOT be applied when
    additional evidence for the variant supports that the variant causes a
    phenotype other than HCM (e.g., variant seen in multiple DCM cases)".
    Verbatim in the PM1 text of MYH7 (GN002), MYBPC3 (GN095), TNNI3 (GN098) and
    TNNT2 (GN099) — and in no other gene, which is why the caller keys on the
    spec sentence rather than a hardcoded gene list.

    HONEST SCOPE. The spec's condition is evidence about the VARIANT; the
    strongest such signal HeartVar holds is the PROBAND's own phenotype, which
    is a narrower thing. So this suppresses PM1 only on an unambiguous signal:
    a DCM term present and no HCM term present. A proband carrying both, or
    carrying only a generic unspecified-cardiomyopathy term (HP:0001638, in
    neither list), is not evidence for a non-HCM phenotype and is left alone.
    Suppress-only, so the failure mode is a missed +2, never an over-call.
    """
    terms, text = _phenotype_tokens(hpo)
    if not terms and not text.strip():
        return None
    if (terms & _HCM_HPO_SET) or _HCM_TEXT_RE.search(text):
        return None
    dcm_hits = sorted(terms & _DCM_HPO_SET)
    if not dcm_hits and not _DCM_TEXT_RE.search(text):
        return None
    what = ", ".join(dcm_hits) if dcm_hits else text.strip()[:60]
    return (
        "the proband's phenotype supports dilated cardiomyopathy rather than "
        f"HCM ({what})"
    )


def _gate_criteria_applicability(
    criteria: list[dict], gene: str | None, evidence: dict | None = None,
    hpo: str | list[str] | None = None,
) -> list[dict]:
    """Demote any met _APPLICABILITY_GATED_CODES criterion (PP2/BP1/PP4/PM1) to
    not_met when the gene's VCEP marks it Not Applicable (CSpec). Spec-faithful
    precision guard for LLM-assigned codes — e.g. removes the spurious PP2 the
    model fires on sarcomeric genes, PP4 on RASopathy genes, or PM1 on genes
    whose CSpec disallows it (MYL2 etc.).

    Additionally enforces BP1's missense-only definition: BP1 ("missense in a
    gene where truncating is the mechanism") is by definition inapplicable to a
    non-missense variant, so a met BP1 is demoted when VEP resolved an effect
    that is NOT ``missense_variant`` (e.g. the LLM fired BP1 on a stop_gained).
    Consequence is read from the in-scope ``evidence`` dict's VEP block."""
    most_severe = (
        ((evidence or {}).get("vep") or {}).get("most_severe_consequence") or ""
    ).lower()
    met_codes = {
        str(x.get("code") or "").upper()
        for x in criteria if x.get("status") == "met"
    }
    pm1_strong_bars_pm5_ps1 = False
    for c in criteria:
        code = c.get("code")
        if c.get("status") != "met":
            continue
        if code in _APPLICABILITY_GATED_CODES and not _criterion_applicable(gene, code):
            gn = (_VCEP_CRIT.get((gene or "").upper()) or {}).get("gn", "CSpec")
            c["status"] = "not_met"
            c["criteria_strength"] = None
            c["evidence"] = (
                f"{gene} VCEP marks {code} Not Applicable ({gn}) — suppressed "
                f"(was: {str(c.get('evidence') or '')[:160]})"
            )
            continue
        if code == "PP4":
            _pp4_spec = (_VCEP_CRIT.get((gene or "").upper()) or {})
            _pp4_affirmative = (
                (_pp4_spec.get("PP4") or {}).get("applicability") == "applicable"
            )
            if not _pp4_affirmative and (gene or "").upper() not in _PP4_EVIDENCE_ALLOW_LIST:
                c["status"] = "not_met"
                c["criteria_strength"] = None
                c["evidence"] = (
                    f"PP4 is applied only where a VCEP affirmatively marks it "
                    f"applicable (KCNQ1 GN112) or expert practice establishes it "
                    f"(FBN1); {gene or 'this gene'} does neither, and PP4 has no "
                    "diagnostic-yield asset behind it here — suppressed "
                    f"(was: {str(c.get('evidence') or '')[:160]})"
                )
                continue

        if (code == "BP1"
                and most_severe
                and most_severe != "missense_variant"):
            c["status"] = "not_met"
            c["criteria_strength"] = None
            c["evidence"] = (
                f"BP1 is missense-only by definition but the variant is "
                f"{most_severe} (not missense) — suppressed "
                f"(was: {str(c.get('evidence') or '')[:160]})"
            )
        if (code == "PP2"
                and most_severe
                and most_severe != "missense_variant"):
            c["status"] = "not_met"
            c["criteria_strength"] = None
            c["evidence"] = (
                f"PP2 is missense-only by definition but the variant is "
                f"{most_severe} (not missense) — suppressed "
                f"(was: {str(c.get('evidence') or '')[:160]})"
            )
        if code == "PM1" and most_severe in _PM1_LOF_CONSEQUENCES:
            c["status"] = "not_met"
            c["criteria_strength"] = None
            c["evidence"] = (
                f"PM1 applies to missense / in-frame variants; this variant is "
                f"{most_severe} (loss-of-function — covered by PVS1) — suppressed "
                f"(was: {str(c.get('evidence') or '')[:140]})"
            )
            continue
        if code == "PM1" and c.get("status") == "met":
            pm1_text = _crit_text(gene, "PM1")
            if "phenotype other than hcm" in pm1_text.lower():
                _why = _pm1_non_hcm_phenotype(hpo)
                if _why:
                    c["status"] = "not_met"
                    c["criteria_strength"] = None
                    c["evidence"] = (
                        f"{gene} VCEP derived its PM1 hotspot ranges from HCM "
                        "case cohorts and bars PM1 when the evidence supports a "
                        f"phenotype other than HCM — {_why} — suppressed "
                        f"(was: {str(c.get('evidence') or '')[:140]})"
                    )
                    continue
            if "cbegf" in pm1_text.lower():
                _vep = (evidence or {}).get("vep") or {}
                _hgvsp = _mane_hgvsp_from_vep(_vep) or _vep.get("hgvsp") or ""
                _verdict, _why = _pm1_fbn1(evidence, _hgvsp)
                _gn = _crit_spec(gene, "PM1") and (
                    (_VCEP_CRIT.get((gene or "").upper()) or {}).get("gn", "CSpec")
                )
                if _verdict == "barred":
                    c["status"] = "not_met"
                    c["criteria_strength"] = None
                    c["evidence"] = (
                        f"{gene} VCEP ({_gn}) bars PM1 here — {_why} — suppressed "
                        f"(was: {str(c.get('evidence') or '')[:140]})"
                    )
                    continue
                if _verdict == "no_rule":
                    c["status"] = "not_met"
                    c["criteria_strength"] = None
                    c["evidence"] = (
                        f"{gene} VCEP ({_gn}) enumerates the residues PM1 covers "
                        f"and this variant is not among them — {_why} — suppressed "
                        f"(was: {str(c.get('evidence') or '')[:140]})"
                    )
                    continue
                if _verdict in ("Strong", "Moderate"):
                    c["criteria_strength"] = f"PM1_{_verdict}"
                    c["evidence"] = (
                        f"{gene} VCEP ({_gn}) applies PM1 at {_verdict}: {_why}."
                        f" (was: {str(c.get('evidence') or '')[:120]})"
                    )
                    if _verdict == "Strong":
                        pm1_strong_bars_pm5_ps1 = True
                    continue

            if "in order to be considered for pm1" in pm1_text.lower() \
                    and "rare" in pm1_text.lower() and "PM2" not in met_codes:
                c["status"] = "not_met"
                c["criteria_strength"] = None
                c["evidence"] = (
                    f"{gene} VCEP requires the variant to be rare (PM2) before "
                    "PM1 can be considered, and PM2 is not met — suppressed "
                    f"(was: {str(c.get('evidence') or '')[:140]})"
                )
                continue
            if pm1_text and "exon" not in pm1_text.lower():
                rngs, _ = _parse_pm1_hotspots(pm1_text)
                pos = _any_protein_position_from_vep((evidence or {}).get("vep") or {})
                if rngs and pos is not None and not any(lo <= pos <= hi for lo, hi in rngs):
                    c["status"] = "not_met"
                    c["criteria_strength"] = None
                    c["evidence"] = (
                        f"{gene} VCEP restricts PM1 to hotspot ranges {rngs}; "
                        f"residue {pos} is outside all of them — suppressed "
                        f"(was: {str(c.get('evidence') or '')[:140]})"
                    )

    most_severe_pm1 = (
        ((evidence or {}).get("vep") or {}).get("most_severe_consequence") or ""
    ).lower()
    pm1_entry = next((c for c in criteria if c.get("code") == "PM1"), None)
    pm1_already_met = bool(pm1_entry and pm1_entry.get("status") == "met")
    pm2_met = any(
        c.get("code") == "PM2" and c.get("status") == "met" for c in criteria
    )
    if (
        not pm1_already_met
        and most_severe_pm1 == "missense_variant"
        and pm2_met
        and _criterion_applicable(gene, "PM1")
    ):
        pm1_text = _crit_text(gene, "PM1")
        _assert_block = (
            _pm1_non_hcm_phenotype(hpo)
            if "phenotype other than hcm" in pm1_text.lower() else None
        )
        if pm1_text and not _assert_block:
            rngs, hotspot_exons = _parse_pm1_hotspots(pm1_text)
            vep = (evidence or {}).get("vep") or {}
            pos = _any_protein_position_from_vep(vep)
            var_exon = _vep_exon_number(vep)
            in_range = pos is not None and any(lo <= pos <= hi for lo, hi in rngs)
            in_exon = var_exon is not None and var_exon in hotspot_exons
            if in_range or in_exon:
                if in_range:
                    hit = next((lo, hi) for lo, hi in rngs if lo <= pos <= hi)
                    where = f"hotspot range {hit} (residue {pos})"
                else:
                    where = f"hotspot exon {var_exon}"
                ev_txt = (
                    f"{gene} VCEP scopes PM1 to enumerated hotspots; rare missense "
                    f"in {where} — PM1_Moderate asserted per the VCEP PM1 "
                    "definition (variant in critical mutational hotspot / "
                    "functional domain)"
                )
                if pm1_entry is not None:
                    pm1_entry["status"] = "met"
                    pm1_entry["criteria_strength"] = "PM1_Moderate"
                    pm1_entry["evidence"] = ev_txt
                else:
                    criteria.append({
                        "code": "PM1",
                        "name": _CRITERION_NAMES.get("PM1", "PM1"),
                        "status": "met",
                        "direction": "pathogenic",
                        "criteria_strength": "PM1_Moderate",
                        "evidence": ev_txt,
                        "source": "vcep_range",
                    })
    if pm1_strong_bars_pm5_ps1:
        for c in criteria:
            if c.get("code") in ("PM5", "PS1") and c.get("status") == "met":
                c["status"] = "not_met"
                c["criteria_strength"] = None
                c["evidence"] = (
                    f"{gene} VCEP bars {c.get('code')} when PM1 is applied at "
                    "Strong for a cysteine in a cbEGF-like domain (same evidence "
                    f"counted twice) — suppressed "
                    f"(was: {str(c.get('evidence') or '')[:140]})"
                )
    return criteria


_CONSTRAINT_PLI_MIN = 0.9
_CONSTRAINT_LOEUF_MAX = 0.35
_CONSTRAINED_PM2_CEILING = 0.00004
_CONSTRAINED_BS1_FLOOR = 0.0001


def _is_lof_constrained(ev: dict) -> bool | None:
    """True when gnomAD marks the gene LoF-intolerant (pLI >= 0.9 or
    LOEUF < 0.35) — the proxy the PVS1 tier uses. None when no gnomAD
    constraint data is available (so callers can keep the generic default)."""
    c = ((ev.get("gnomad") or {}).get("gene") or {}).get("gnomad_constraint") or {}
    pli, loeuf = c.get("pLI"), c.get("oe_lof_upper")
    if pli is None and loeuf is None:
        return None
    return (pli is not None and pli >= _CONSTRAINT_PLI_MIN) or \
           (loeuf is not None and loeuf < _CONSTRAINT_LOEUF_MAX)


def _freq_source_label(ev: dict, gene: str | None, inheritance: str) -> str:
    """Human label for which frequency-threshold basis is in effect, surfaced in
    evidence strings + logged so out-of-table (no-VCEP-spec) genes are visible."""
    if _vcep_freq(gene):
        return "VCEP-specific"
    if inheritance != "AR" and _is_lof_constrained(ev):
        return "constraint-aware fallback (gnomAD LoF-intolerant; no VCEP spec)"
    return f"generic fallback for {inheritance or 'AD/unknown'} inheritance (no VCEP spec)"


_GENERIC_BA1 = 0.05


def _ba1_threshold(gene: str | None) -> float:
    """BA1 stand-alone threshold. Gene-specific VCEP value when published, else
    the ACMG/SVI 0.05 default. NEVER returns below 0.05 for an out-of-table gene
    (the VCEP cutoffs are gene-curated; a generic gene gets the conservative
    5% stand-alone bar).

    gnomAD-version safeguard: the published VCEP cutoffs are gnomAD v2.1.1-
    calibrated but the runtime pipeline queries v4 (see vcep_frequency_thresholds
    `_gnomad_version_note`). The cutoff VALUES are disease-derived (max-credible-AF
    = prevalence x heterogeneity / (2 x penetrance), Whiffin/Ware 2017) and are
    therefore version-independent; the comparison metric (_gnomad_popmax_af) is
    the FAF95 filtering allele frequency, the version-robust statistic the VCEPs
    specify. We never lower the bar below 0.05 for a gene ClinGen has not
    published a CSpec for, so v4 metric drift cannot create a spurious BA1."""
    spec = _vcep_freq(gene)
    if spec and spec.get("ba1") is not None:
        return max(float(spec["ba1"]), 0.0)
    return _GENERIC_BA1


def _bs1_threshold(ev: dict, gene: str | None, inheritance: str) -> float:
    """BS1 floor. Priority: gene-specific VCEP value -> constraint-aware default
    for LoF-intolerant out-of-table AD genes -> prevalence-scaled generic
    fallback."""
    spec = _vcep_freq(gene)
    if spec and spec.get("bs1") is not None:
        return float(spec["bs1"])
    if inheritance == "AR":
        return 0.005
    if _is_lof_constrained(ev):
        return _CONSTRAINED_BS1_FLOOR
    return _bs1_ad_threshold(ev)


def _pm2_threshold(gene: str | None, inheritance: str, ev: dict | None = None) -> float:
    """PM2 ceiling (variant earns PM2_Supporting when FAF is at/below this and
    not precluded by BA1/BS1). Priority: gene-specific VCEP value (incl. the 0.0
    absence-only rule and the LZTR1 AR carve-out pm2_max_ar) -> constraint-aware
    default for LoF-intolerant out-of-table AD genes -> generic (AR 0.01 / AD
    1e-4). `ev` is optional so existing 2-arg callers keep the generic path."""
    spec = _vcep_freq(gene)
    if spec:
        if inheritance == "AR" and spec.get("pm2_max_ar") is not None:
            return float(spec["pm2_max_ar"])
        if spec.get("pm2_max") is not None:
            return float(spec["pm2_max"])
    if inheritance == "AR":
        return 0.01
    if _is_lof_constrained(ev or {}):
        return _CONSTRAINED_PM2_CEILING
    return 0.0001


def _bs1_ad_threshold(ev: dict) -> float:
    """Return the BS1 threshold for AD/unknown inheritance, scaled to
    disease prevalence.

    Priority order:
    1. Phenotype-matched ClinGen tier — used when HPO was submitted and
       matched at least one curated disease. Most precise signal.
    2. Gene-wide ClinGen tier — used as fallback when no HPO was submitted
       (e.g. benchmark runs). Requires corroboration from PanelApp green
       presence AND absence of any Disputed/Refuted record to avoid
       mis-applying a high-prevalence threshold to pleiotropic genes where
       the Definitive curation is for a different disease than the one
       being assessed.
    3. Default 0.001 — conservative fallback when ClinGen data is absent
       or ambiguous.

    Rationale: Richards 2015 defines BS1 as 'allele frequency greater than
    expected for the disorder'. For high-penetrance AD diseases with
    Definitive/Strong ClinGen curation (HCM, LQTS, Noonan), 0.001 is too
    permissive; 0.0004 aligns with ClinGen VCEP practice. For genes with
    weaker or disputed curation, 0.001 is the safe fallback.
    """
    STRONG_TIERS = ("Definitive", "Strong")
    gencc = ev.get("gencc") or {}
    pheno_matched_tier = gencc.get("clingen_phenotype_matched_best_classification") or ""
    if pheno_matched_tier in STRONG_TIERS:
        return 0.0004
    if pheno_matched_tier:
        return 0.001
    gene_wide_tier = gencc.get("best_clingen_classification") or ""
    has_disputed = bool(gencc.get("clingen_has_disputed_or_refuted"))
    panelapp = ev.get("panelapp") or {}
    green_panels = int(panelapp.get("green_panels") or 0)
    if (
        gene_wide_tier in STRONG_TIERS
        and green_panels >= 1
        and not has_disputed
    ):
        return 0.0004
    return 0.001


_DOMINANT_NEGATIVE_GENES = frozenset({
    "PTPN11", "RAF1", "KRAS", "HRAS", "NRAS", "BRAF",
    "MAP2K1", "MAP2K2", "RIT1", "MRAS", "LZTR1",
    "SCN5A",
    "MYH7", "MYH6", "ACTC1", "TPM1", "TPM3",
    "MYL2", "MYL3", "TNNC1", "TNNI3", "TNNT2", "GLA",
})

_LOF_CONSEQUENCE_TYPES: frozenset[str] = frozenset({
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "transcript_ablation",
    "exon_loss_variant",
})

_LOF_INHERITANCE_CODES: frozenset[str] = frozenset(
    {"AD", "AR", "CH", "XL", "XLR", "XLD"}
)
_LOF_GENCC_MOI: tuple[str, ...] = (
    "autosomal dominant",
    "autosomal recessive",
    "x-linked",
)


def _gene_lof_mechanism(
    ev: dict,
    gene_symbol: str = "",
    consequence: str = "",
) -> bool:
    """Three-tier proxy for "loss-of-function is an established disease
    mechanism for this gene". Returns True if ANY of the following
    supports a LoF mechanism:

      1. **CHDgene listing with a Mendelian LoF inheritance** (AD/AR/CH/XL;
         see ``_LOF_INHERITANCE_CODES``) — CHDgene curates only
         high-confidence CHD genes; such a listing is a strong LoF proxy.
      2. **GenCC Definitive or Strong with a Mendelian LoF MoI**
         (autosomal dominant / autosomal recessive / X-linked; see
         ``_LOF_GENCC_MOI``) — a submission at that tier carries
         expert-level support that LoF is the disease mechanism. Recessive
         and X-linked are admitted because PVS1 fires on the null *allele*;
         zygosity (PM3) governs the genotype-level call, not PVS1 itself.
      3. **gnomAD constraint** (pLI > 0.9 or LOEUF < 0.35) — the
         population-genetics LoF-intolerance proxy retained as fallback.

    Why gnomAD constraint alone is insufficient for cardiac genes: KCNQ1,
    LZTR1, and similar dominant cardiac genes cause disease via LoF yet
    tolerate enough population-level LoF variation that they don't cross
    the constraint thresholds. The CHDgene and GenCC tiers add the
    curator-level evidence that the constraint proxy misses.

    Mechanism exclusion is consequence-aware. Genes in
    ``_DOMINANT_NEGATIVE_GENES`` short-circuit the CHDgene/GenCC tiers
    ONLY when the variant consequence is not in ``_LOF_CONSEQUENCE_TYPES``.
    The dominant-negative/GoF exclusion is meant to prevent PVS1 on
    missense variants where the GoF/dominant-negative mechanism applies;
    it should not suppress PVS1 universally for these genes. A frameshift
    in MYH7 or a stop_gained in SCN5A is a genuine LoF event for which
    PVS1 routing remains correct — those are canonical LoF, so
    ``tiers_enabled`` is True and the gnomAD constraint applies as normal.

    The gnomAD constraint tier is gated by ``tiers_enabled`` for the same
    reason: population-level LoF intolerance (pLI > 0.9 / LOEUF < 0.35)
    measures selection against LoF alleles, which is **orthogonal** to
    whether LoF is the *disease* mechanism. BRAF/HRAS and other RASopathy
    genes are strongly LoF-constrained yet cause disease through
    gain-of-function — leaving the constraint tier unconditional fired a
    spurious PVS1 (via the SpliceAI arm) on their missense variants. When
    the gene is in the GoF/dominant-negative exclusion AND the consequence
    is non-canonical-LoF, no tier (curator or constraint) may assert a LoF
    mechanism.
    """
    gene_upper = (gene_symbol or "").strip().upper()
    is_canonical_lof = any(
        c in (consequence or "") for c in _LOF_CONSEQUENCE_TYPES
    )
    if not _criterion_applicable(gene_upper, "PVS1"):
        return False
    tiers_enabled = (gene_upper not in _DOMINANT_NEGATIVE_GENES) or is_canonical_lof

    chdgene = ev.get("chdgene") or {}
    if tiers_enabled and chdgene.get("listed") and (
        _LOF_INHERITANCE_CODES & set(chdgene.get("inheritance") or [])
    ):
        return True

    if tiers_enabled:
        gencc = ev.get("gencc") or {}
        for sub in (gencc.get("submissions") or []):
            classification = (sub.get("classification") or "").strip()
            moi = (sub.get("moi") or "").lower()
            if classification in ("Definitive", "Strong") and any(
                m in moi for m in _LOF_GENCC_MOI
            ):
                return True

    if tiers_enabled:
        constraint = (
            ((ev.get("gnomad") or {}).get("gene") or {}).get("gnomad_constraint")
            or {}
        )
        pli = constraint.get("pLI")
        loeuf = constraint.get("oe_lof_upper")
        if pli is not None:
            try:
                if float(pli) > 0.9:
                    return True
            except (TypeError, ValueError):
                pass
        if loeuf is not None:
            try:
                if float(loeuf) < 0.35:
                    return True
            except (TypeError, ValueError):
                pass
    return False


MECHANISM_LABELS = {
    "haploinsufficiency": "Loss-of-function (haploinsufficiency)",
    "recessive_lof": "Loss-of-function (recessive / biallelic)",
    "dominant_negative": "Dominant-negative",
    "gof_or_dn": "Gain-of-function or dominant-negative",
    "mixed": "Mixed mechanism (LoF and GoF/DN both reported)",
    "undetermined": "Mechanism not established",
}

_MIS_Z_CONSTRAINED = 3.09


def _mechanism_signals(ev: dict, gene_upper: str) -> dict:
    """Collect every raw per-gene mechanism signal HeartVar holds, from the
    ClinGen dosage cache plus the live evidence already gathered for the gene.
    The classifier and the flag both read from here so they can never drift.
    """
    rec = _GENE_DOSAGE.get(gene_upper) or {}
    hi_raw = rec.get("hi_score")
    hi = int(hi_raw) if (hi_raw or "").strip().isdigit() else None
    spec = _VCEP_CRIT.get(gene_upper) or {}
    pvs1_app = (spec.get("PVS1") or {}).get("applicability")
    pp2_app = (spec.get("PP2") or {}).get("applicability")
    chdgene = ev.get("chdgene") or {}
    chd_inh = {str(c).strip().upper() for c in (chdgene.get("inheritance") or [])}
    constraint = (
        ((ev.get("gnomad") or {}).get("gene") or {}).get("gnomad_constraint")
        or {}
    )
    return {
        "hi": hi,
        "hi_raw": hi_raw,
        "hi_desc": rec.get("hi_desc"),
        "clingen_url": rec.get("clingen_url"),
        "vcep_name": spec.get("vcep"),
        "pvs1_na": pvs1_app == "not_applicable",
        "pvs1_applicable": pvs1_app == "applicable",
        "pp2_applicable": pp2_app == "applicable",
        "chdgene_dn": "DN" in chd_inh,
        "hand_dn": gene_upper in _DOMINANT_NEGATIVE_GENES,
        "mis_z": _as_float(constraint.get("mis_z")),
        "loeuf": _as_float(constraint.get("oe_lof_upper")),
        "pli": _as_float(constraint.get("pLI")),
    }


def gene_mechanism(ev: dict, gene: str = "") -> dict:
    """Deterministic, source-attributed disease-mechanism classification.

    Answers the curator question "is this gene's disease caused by loss- or
    gain-of-function?" by merging the curated signals HeartVar holds — never
    by guessing. Two axes are assembled:

      * **LoF axis** — ClinGen Gene-Dosage haploinsufficiency score 3
        (sufficient evidence) or 30 (recessive/biallelic), or a VCEP that keeps
        PVS1 *Applicable*.
      * **non-LoF axis** — a VCEP that marks PVS1 *Not Applicable*, ClinGen
        dosage 40 (dosage sensitivity unlikely), the CHDgene dominant-negative
        inheritance flag, or membership of HeartVar's curated GoF/DN list.

    When **both** axes fire the gene is reported ``mixed`` (e.g. SCN5A:
    LoF→Brugada, GoF→LQT3; PTPN11: GoF→Noonan, LoF→metachondromatosis) — the
    honest answer, and the guard that stops the consistency flag mis-firing on
    these genes. gnomAD constraint is reported only as a hedged *population
    proxy*, never as the mechanism itself, and a gene with no curated signal is
    ``undetermined`` rather than mislabelled.

    Returns ``mechanism`` / ``label`` / ``confidence`` / ``sources`` /
    ``hi_score`` / ``summary`` / ``constraint_hint`` for display. Informational
    only — it does not change any ACMG criterion or score.
    """
    gene_upper = (gene or "").strip().upper()
    if not gene_upper:
        return {
            "mechanism": "undetermined",
            "label": MECHANISM_LABELS["undetermined"],
            "confidence": "undetermined",
            "sources": [],
            "hi_score": None,
            "summary": "No gene specified.",
            "constraint_hint": None,
        }
    s = _mechanism_signals(ev, gene_upper)
    hi = s["hi"]

    sources: list[dict] = []
    if hi is not None:
        detail = f"Haploinsufficiency score {hi}"
        if s["hi_desc"]:
            detail += f" — {s['hi_desc']}"
        src = {"name": "ClinGen Dosage", "detail": detail}
        if s["clingen_url"]:
            src["url"] = s["clingen_url"]
        sources.append(src)
    vcep_label = f"{s['vcep_name']} VCEP" if s["vcep_name"] else "ClinGen VCEP"
    if s["pvs1_na"]:
        sources.append({"name": vcep_label, "detail": "PVS1 (null-variant rule) Not Applicable — loss of function is not the established mechanism"})
    elif s["pvs1_applicable"]:
        sources.append({"name": vcep_label, "detail": "PVS1 Applicable — loss of function is an accepted mechanism"})
    if s["chdgene_dn"]:
        sources.append({"name": "CHDgene", "detail": "Curated dominant-negative inheritance"})
    if s["hand_dn"] and not s["chdgene_dn"]:
        sources.append({"name": "HeartVar curated", "detail": "Listed as a gain-of-function / dominant-negative gene"})

    lof_haplo = (hi == 3) or (s["pvs1_applicable"] and hi != 30)
    recessive_lof = hi == 30
    lof_emerging = hi == 2
    nonlof_strong = s["pvs1_na"] or hi == 40
    nonlof_curated = s["chdgene_dn"] or s["hand_dn"]
    nonlof_any = nonlof_strong or nonlof_curated

    constraint_hint = None
    if s["mis_z"] is not None and s["mis_z"] > _MIS_Z_CONSTRAINED:
        constraint_hint = (
            f"gnomAD missense-constrained (Z={s['mis_z']:.1f}); population "
            f"proxy consistent with gain-of-function / dominant-negative"
        )
    elif s["loeuf"] is not None and s["loeuf"] < 0.35:
        constraint_hint = (
            f"gnomAD LoF-intolerant (LOEUF={s['loeuf']:.2f}); population "
            f"proxy consistent with haploinsufficiency"
        )
    elif s["pli"] is not None and s["pli"] > 0.9:
        constraint_hint = (
            f"gnomAD LoF-intolerant (pLI={s['pli']:.2f}); population proxy "
            f"consistent with haploinsufficiency"
        )

    if (lof_haplo or recessive_lof or lof_emerging) and nonlof_any:
        mech, conf = "mixed", ("established" if (lof_haplo or recessive_lof) else "emerging")
    elif nonlof_any:
        mech = "dominant_negative" if s["chdgene_dn"] else "gof_or_dn"
        conf = "established" if nonlof_strong else "curated"
    elif lof_haplo:
        mech, conf = "haploinsufficiency", "established"
    elif recessive_lof:
        mech, conf = "recessive_lof", "established"
    elif lof_emerging:
        mech, conf = "haploinsufficiency", "emerging"
    else:
        mech, conf = "undetermined", "undetermined"

    label = MECHANISM_LABELS[mech]
    if mech == "undetermined":
        summary = "Established disease mechanism not determined from curated sources."
        if constraint_hint:
            summary += f" {constraint_hint}."
    elif mech == "mixed":
        summary = (
            f"Both loss- and gain-of-function (or dominant-negative) disease have been reported "
            f"for {gene_upper}; the operative mechanism depends on the specific variant and phenotype."
        )
    else:
        summary = label + "."
    return {
        "mechanism": mech,
        "label": label,
        "confidence": conf,
        "sources": sources,
        "hi_score": s["hi_raw"],
        "summary": summary,
        "constraint_hint": constraint_hint,
    }


def mechanism_consistency_flag(
    evidence: dict, gene: str = "", mech: dict | None = None
) -> dict | None:
    """Variant-type ↔ disease-mechanism consistency note — asymmetric by design,
    informational only (never changes a score).

    Three directions, in priority order:

      1. **Mixed-mechanism gene** (e.g. SCN5A, PTPN11) → a transparency ``info``
         note that the consistency check is deliberately *not* applied; no
         variant type is inconsistent when both mechanisms operate.
      2. **Null / LoF variant in a non-LoF gene** → the genuinely suspicious
         direction. ``warning`` when expert/curated-strong evidence says LoF is
         not the mechanism (VCEP PVS1 Not Applicable, or ClinGen dosage 40);
         softer ``info`` when only HeartVar's hand-list / CHDgene-DN flag
         (no VCEP/ClinGen confirmation) supports it.
      3. **Missense in a dominant-haploinsufficiency gene** (ClinGen HI=3) → a
         soft ``info`` mechanism-fit note: a missense here causes disease only
         if it yields loss of function. **Suppressed** when missense is plausibly
         the driver — the gene is missense-constrained (mis_z>3.09) or its VCEP
         keeps PP2 (missense common) Applicable — so genuine pathogenic missense
         in genes like KCNQ1 is never mis-flagged.

    Returns a structured dict (``kind``/``severity``/``title``/``detail``/
    ``gene``/``consequence``/``mechanism``) or ``None``. Genes with no curated
    mechanism signal (``undetermined``) always return ``None``.
    """
    gene_upper = (gene or "").strip().upper()
    if not gene_upper:
        return None
    vep = (evidence or {}).get("vep") or {}
    consequence = (vep.get("most_severe_consequence") or "").lower()
    if not consequence:
        return None
    if mech is None:
        mech = gene_mechanism(evidence, gene_upper)
    m = mech.get("mechanism")
    s = _mechanism_signals(evidence, gene_upper)
    is_canonical_lof = any(c in consequence for c in _LOF_CONSEQUENCE_TYPES)
    is_missense = "missense_variant" in consequence
    pretty_csq = consequence.replace("_variant", "").replace("_", " ").strip()

    if m == "mixed" and (is_canonical_lof or is_missense):
        return {
            "kind": "mixed_mechanism",
            "severity": "info",
            "gene": gene_upper,
            "consequence": consequence,
            "mechanism": m,
            "title": "Mixed disease mechanism",
            "detail": (
                f"{gene_upper} has both loss-of-function and gain-of-function (or dominant-negative) "
                f"disease reported, so the variant-type consistency check is not applied — the operative "
                f"mechanism depends on the specific variant and phenotype."
            ),
        }

    if is_canonical_lof and m in ("gof_or_dn", "dominant_negative"):
        if s["pvs1_na"]:
            return {
                "kind": "lof_in_non_lof_gene",
                "severity": "warning",
                "gene": gene_upper,
                "consequence": consequence,
                "mechanism": m,
                "title": "Variant type vs established disease mechanism",
                "detail": (
                    f"This is a predicted loss-of-function variant ({pretty_csq}), but the ClinGen "
                    f"{s['vcep_name'] or ''} Variant Curation Expert Panel for {gene_upper} marks PVS1 "
                    f"(the null-variant rule) Not Applicable — loss of function is not the established "
                    f"disease mechanism for this gene. Re-assess: a truncating variant here may not be "
                    f"pathogenic by the gene's mechanism, and PVS1 is not applied."
                ).replace("  ", " "),
            }
        if s["hi"] == 40:
            return {
                "kind": "lof_in_non_lof_gene",
                "severity": "warning",
                "gene": gene_upper,
                "consequence": consequence,
                "mechanism": m,
                "title": "Variant type vs established disease mechanism",
                "detail": (
                    f"This is a predicted loss-of-function variant ({pretty_csq}), but ClinGen Gene-Dosage "
                    f"scores {gene_upper} haploinsufficiency 40 (dosage sensitivity unlikely) — loss of "
                    f"function is not an established disease mechanism. Re-assess whether a truncating "
                    f"variant is pathogenic here; PVS1 is not applied."
                ),
            }
        # confirmation → softer, attributed to the curated list.
        return {
            "kind": "lof_in_dn_gene",
            "severity": "info",
            "gene": gene_upper,
            "consequence": consequence,
            "mechanism": m,
            "title": "Variant type vs established disease mechanism",
            "detail": (
                f"This is a predicted loss-of-function variant ({pretty_csq}). {gene_upper} is curated as a "
                f"primarily gain-of-function / dominant-negative gene, where loss of function is not the "
                f"established disease mechanism — a truncating variant may not be pathogenic by that mechanism. "
                f"Confirm against gene-specific guidance."
            ),
        }

    if (
        is_missense
        and m == "haploinsufficiency"
        and mech.get("confidence") == "established"
        and s["hi"] == 3
    ):
        if s["mis_z"] is not None and s["mis_z"] > _MIS_Z_CONSTRAINED:
            return None
        if s["pp2_applicable"]:
            return None
        return {
            "kind": "missense_in_lof_gene",
            "severity": "info",
            "gene": gene_upper,
            "consequence": consequence,
            "mechanism": m,
            "title": "Variant type vs established disease mechanism",
            "detail": (
                f"The established disease mechanism for {gene_upper} is haploinsufficiency (ClinGen "
                f"sufficient-evidence dosage call). Pathogenic missense is less common in such genes and "
                f"typically acts by reducing protein function (destabilising / null-like) or via a distinct "
                f"gene-specific mechanism — weigh functional and in-silico evidence rather than inferring "
                f"pathogenicity from variant type alone. This is mechanism context, not a benign call."
            ),
        }

    return None


_REC_MODES = frozenset({"AR", "XLR", "CH"})
_DOM_MODES = frozenset({"AD", "XLD", "SD", "DN"})
_XL_MODES = frozenset({"XL", "XLR", "XLD"})

# (partly OMIM-derived) feed, so dropping it keeps OMIM-licensed MOI out of the
_ADEQUATE_TIERS = frozenset({"Definitive", "Strong", "Moderate"})

_PAR_GENES = frozenset({
    "SHOX", "CRLF2", "CSF2RA", "IL3RA", "SLC25A6", "ASMT", "P2RY8", "AKAP17A",
    "ASMTL", "DHRSX", "ZBED1", "CD99", "XG", "GTPBP6", "PLCXD1", "PPP2R3B",
    "VAMP7", "SPRY3", "IL9R", "WASH6P",
})


def _norm_gencc_moi(moi: str | None) -> str | None:
    """GenCC free-text MoI → canonical mode. None for Unknown / unrecognised
    (fail closed → the caller treats it as no signal, never a wrong one)."""
    m = (moi or "").strip().lower()
    if not m:
        return None
    if "x-linked" in m or "x linked" in m:
        if "recessive" in m:
            return "XLR"
        if "dominant" in m:
            return "XLD"
        return "XL"
    if "recessive" in m:
        return "AR"
    if "semidominant" in m or "semi-dominant" in m:
        return "SD"
    if "dominant" in m:
        return "AD"
    if "mitochond" in m:
        return "MT"
    if "y-linked" in m or "y linked" in m:
        return "YL"
    return None


def _norm_clingen_gv_moi(moi: str | None) -> str | None:
    """ClinGen GV code (AR/AD/XL/SD/MT/UD) → canonical mode (UD → None)."""
    return {
        "AR": "AR", "AD": "AD", "XL": "XL", "SD": "SD", "MT": "MT", "UD": None,
    }.get((moi or "").strip().upper())


def _norm_panelapp_moi(moi: str | None) -> list[str]:
    """PanelApp free-text MoI → canonical modes (may be several, e.g. BOTH…).
    X-linked entries describe male-hemizygous + female-biallelic in prose, so
    we map them to XL only (never let that prose leak an AR/AD signal)."""
    m = (moi or "").strip().lower()
    if not m:
        return []
    if "x-linked" in m or "x linked" in m:
        return ["XL"]
    out: list[str] = []
    if "biallelic" in m:
        out.append("AR")
    if "monoallelic" in m:
        out.append("AD")
    if "mitochond" in m:
        out.append("MT")
    return out


def gene_inheritance_modes(ev: dict, gene: str = "") -> dict:
    """Roll up mode of inheritance across the free, license-clean sources HeartVar
    holds — ClinGen Gene-Disease Validity (authoritative, per gene-disease pair) +
    GenCC + PanelApp + CHDgene + ClinGen Gene-Dosage — behind a symmetric adequate-
    tier evidence gate, for the carrier-status determination.

    The determinative recessive/dominant signal is captured PER DISEASE, not just
    rolled up to the gene, in ``facts`` — a list of ``{disease, mode, category,
    tier, source}`` from ClinGen GV + adequate-tier GenCC (the high-confidence,
    per-gene-disease sources). ``carrier_status`` matches the proband phenotype
    against ``category`` and scopes the determination to the relevant disease(s),
    so a gene's non-cardiac recessive association (e.g. MYH7's myosin-storage
    myopathy) never bleeds into a cardiac proband's carrier interpretation. With no
    phenotype (or no confident match) it falls back to the gene-wide picture and
    discloses WHICH conditions are dominant vs recessive rather than collapsing to
    a single verdict.

    The gene-level booleans remain for that fallback + the safety invariants.
    Conservative toward dominant: ``has_dominant_signal`` is True for ANY adequate-
    tier dominant mode OR a gene-level dominant indicator (known DN/GoF gene,
    ClinGen haploinsufficiency = 3, CHDgene AD/DN/XLD). That floor can only ADD a
    dominant signal, never remove one — so a dual-MOI cardiac gene (CASQ2, KCNQ1,
    TTN, JPH2 …) can never be mis-called a bare unaffected carrier.

    PanelApp is intentionally NOT a per-disease source: its "BOTH monoallelic and
    biallelic" panel boilerplate is too coarse to attach to a specific disease, so
    it only corroborates the gene-level dominant/X-linked floor here.
    """
    gene_upper = (gene or "").strip().upper()
    ev = ev or {}

    rec_adequate = dom_adequate = xl_any = mt_any = False
    rec_disease = dom_disease = None
    sources: list[str] = []
    facts: list[dict] = []

    def _note(name: str) -> None:
        if name not in sources:
            sources.append(name)

    def _absorb(mode, disease, src, per_disease=False):
        nonlocal rec_adequate, dom_adequate, xl_any, mt_any, rec_disease, dom_disease
        if not mode or mode == "YL":
            return
        if mode == "MT":
            mt_any = True
            return
        _note(src)
        if per_disease and disease:
            facts.append({
                "disease": disease,
                "mode": mode,
                "category": _categorize_panel(disease),
                "source": src,
            })
        if mode in _XL_MODES:
            xl_any = True
        if mode in _REC_MODES:
            rec_adequate = True
            if not rec_disease and disease:
                rec_disease = disease
        if mode in _DOM_MODES:
            dom_adequate = True
            if not dom_disease and disease:
                dom_disease = disease

    for a in (_CLINGEN_GV.get(gene_upper) or []):
        if (a.get("classification") or "").strip() not in _ADEQUATE_TIERS:
            continue
        _absorb(_norm_clingen_gv_moi(a.get("moi")), a.get("disease") or "",
                "ClinGen Gene-Disease Validity", per_disease=True)

    gencc = ev.get("gencc") or {}
    for sub in (gencc.get("submissions") or []):
        if (sub.get("submitter") or "").strip().lower() == "orphanet":
            continue
        if (sub.get("classification") or "").strip() not in _ADEQUATE_TIERS:
            continue
        _absorb(_norm_gencc_moi(sub.get("moi")), sub.get("disease") or "",
                "GenCC", per_disease=True)

    pa = ev.get("panelapp") or {}
    for panel in (pa.get("panels_found") or []):
        if (panel.get("confidence") or "").strip().lower() != "green":
            continue
        for mode in _norm_panelapp_moi(panel.get("moi")):
            _absorb(mode, panel.get("panel_name") or "", "PanelApp")

    chd = ev.get("chdgene") or {}
    chd_inh = set(chd.get("inheritance") or []) if chd.get("listed") else set()
    if chd_inh:
        _note("CHDgene")
    chd_rec = bool(chd_inh & {"AR", "CH", "XLR"})
    if chd_inh & {"XL", "XLR", "XLD"}:
        xl_any = True
    chd_dom = bool(chd_inh & {"AD", "DN", "XLD"})

    dosage = _GENE_DOSAGE.get(gene_upper) or {}
    hi = (dosage.get("hi_score") or "").strip()
    dosage_recessive = hi == "30"
    dosage_dominant = hi == "3"
    if dosage_recessive:
        _note("ClinGen Gene-Dosage")

    dominant_gene_level = (
        gene_upper in _DOMINANT_NEGATIVE_GENES or dosage_dominant or chd_dom
    )
    if gene_upper in _DOMINANT_NEGATIVE_GENES:
        _note("HeartVar dominant-negative/GoF gene list")

    has_dominant_signal = dom_adequate or dominant_gene_level
    return {
        "recessive_established": rec_adequate,
        "recessive_capable": rec_adequate or dosage_recessive or chd_rec,
        "has_dominant_signal": has_dominant_signal,
        "dominant_any": has_dominant_signal,
        "xlinked": xl_any,
        "mito_only": mt_any and not (
            rec_adequate or dom_adequate or xl_any or dominant_gene_level
        ),
        "rec_disease": rec_disease,
        "dom_disease": dom_disease,
        "facts": facts,
        "dominant_floor": dominant_gene_level,
        "sources": sources,
        "snapshot": {
            "clingen_gv_built_at": _CLINGEN_GV_META.get("built_at"),
            "clingen_gv_present": bool(_CLINGEN_GV),
            "gencc_found": bool(gencc.get("found")),
            "panelapp_present": bool(pa.get("panels_found")),
        },
    }


def carrier_status(evidence: dict, clinical_context: dict, gene: str = "",
                   classification: str = "", submitted_hpo: str = "") -> dict | None:
    """Deterministic, INFORMATIONAL carrier-status note for the clinical summary.

    Joins the variant's FINAL classification (P/LP gate), the proband genotype
    (zygosity / in-trans / sex / chromosome from ``clinical_context``) and the
    per-DISEASE MOI facts from ``gene_inheritance_modes``. Returns a flag dict
    (``kind``/``state``/``severity``/``gene``/``title``/``detail``/``basis``) or
    ``None`` when no statement applies. It NEVER changes the ACMG tier.

    Phenotype scoping is the core of the accurate call. Each ClinGen-GV / GenCC
    fact carries a cardiovascular category; when the proband's phenotype matches a
    fact's category, the recessive/dominant determination is scoped to that
    disease. So MYH7 + a cardiac phenotype resolves to its DOMINANT cardiomyopathy
    (the non-cardiac recessive myosin-storage myopathy doesn't match), while
    DSP / DSG2 + an ARVC phenotype resolve to genuine DUAL inheritance (the
    recessive IS the matched arrhythmogenic-cardiomyopathy). With no phenotype (or
    no confident match) it falls back to the gene-wide picture and DISCLOSES which
    conditions are dominant vs recessive rather than reassuring the curator.

    Safety posture: the only reassuring "carrier" branch requires a phenotype-
    matched recessive disease with NO dominant signal (incl. the gene-level dominant
    floor) at a heterozygous genotype; every other branch is a biallelic/affected
    statement, an explicit dominant/dual disclosure, or an honest UNCERTAIN.
    """
    gene_upper = (gene or "").strip().upper()
    if not gene_upper:
        return None
    classification = (classification or "").strip()
    if classification not in ("Pathogenic", "Likely pathogenic", "VUS"):
        return None

    cc = clinical_context or {}
    inh = (cc.get("inheritance_input") or "").strip().upper()
    zyg = (cc.get("zygosity") or "").strip().lower()
    in_trans = (cc.get("in_trans_pathogenic") or "").strip().lower()
    denovo = (cc.get("denovo_status") or "").strip().lower()
    sex = (cc.get("proband_sex") or "").strip().lower()
    chrom = (cc.get("chromosome") or "").strip().upper()
    zyg_inferred = bool(cc.get("zygosity_inferred"))
    has_phenotype = bool((submitted_hpo or "").strip())

    if chrom in ("MT", "Y") or inh == "MT":
        return None

    modes = gene_inheritance_modes(evidence, gene_upper)
    if modes["mito_only"]:
        return None

    facts = modes.get("facts") or []
    matched = [
        f for f in facts
        if f.get("category") and check_hpo_relevance(submitted_hpo, f["category"])
    ] if has_phenotype else []
    scoped = bool(matched)
    scope_facts = matched if scoped else facts

    def _names(diseases: list[str], fallback: str) -> str:
        seen: list[str] = []
        for d in diseases:
            if d and d not in seen:
                seen.append(d)
        seen = seen[:3]
        if not seen:
            return fallback
        if len(seen) == 1:
            return seen[0]
        return ", ".join(seen[:-1]) + f" and {seen[-1]}"

    rec_diseases = [f["disease"] for f in scope_facts if f["mode"] in _REC_MODES]
    dom_diseases = [f["disease"] for f in scope_facts if f["mode"] in _DOM_MODES]
    rec_present = bool(rec_diseases)
    dom_present = bool(dom_diseases) or modes["dominant_floor"]

    rec_name = _names(rec_diseases, "the associated recessive disorder")
    dom_name = _names(dom_diseases, "a dominant mechanism")

    vus = classification == "VUS"
    pretty = {
        "Pathogenic": "Pathogenic", "Likely pathogenic": "Likely Pathogenic",
    }.get(classification, "VUS")
    var_desc = "variant of uncertain significance" if vus else f"{pretty} variant"

    def _c(assertive: str, conditional: str) -> str:
        """Assertive phrasing (P/LP) vs conditional phrasing (VUS) for a claim."""
        return conditional if vus else assertive

    disclose = "" if scoped else (
        " No phenotype was provided, so this is the gene-wide inheritance picture "
        "across all associated conditions."
        if not has_phenotype else
        " The entered phenotype did not match a specific curated disease for this "
        "gene, so this is the gene-wide inheritance picture."
    )

    pm3_biallelic = (
        inh == "AR" and zyg == "het"
        and (in_trans == "yes" or denovo == "inherited_affected")
    )
    homozygous = (zyg == "hom")

    _CAVEAT = " Informational only; it does not change the ACMG classification."

    def _band(state: str, severity: str, title: str, detail: str) -> dict:
        return {
            "kind": "carrier_status",
            "state": state,
            "severity": "info" if vus else severity,
            "gene": gene_upper,
            "classification": pretty,
            "title": title,
            "detail": detail.rstrip() + _CAVEAT,
            "basis": {
                "state": state,
                "sources": modes["sources"],
                "snapshot": modes["snapshot"],
                "zygosity": zyg or "(not provided)",
                "zygosity_inferred": zyg_inferred,
                "in_trans_pathogenic": in_trans or "(not provided)",
                "proband_sex": sex or "(not provided)",
                "inheritance_input": inh or "(not provided)",
                "phenotype_provided": has_phenotype,
                "phenotype_scoped": scoped,
                "recessive_conditions": sorted(set(rec_diseases)),
                "dominant_conditions": sorted(set(dom_diseases)),
            },
        }

    xlinked = (
        (modes["xlinked"] or inh in ("XLR", "XLD") or chrom == "X")
        and gene_upper not in _PAR_GENES
    )
    if xlinked:
        if sex == "male" or zyg == "hemi":
            return _band(
                "xlinked_male", "info", "X-linked — hemizygous (male)",
                f"{gene_upper} is X-linked and the proband is male, so this {var_desc} "
                f"is hemizygous (a single allele is the whole genotype)"
                + _c(
                    " — not carrier-level.",
                    "; if reclassified pathogenic it would be affected-level, not carrier.",
                ),
            )
        if homozygous:
            return _band(
                "biallelic_affected", "info", "Biallelic genotype (homozygous)",
                f"Homozygous {var_desc} in {gene_upper} — a biallelic genotype"
                + _c(
                    f" for {rec_name}, not carrier-level.",
                    f"; if reclassified pathogenic this would be consistent with affected "
                    f"for {rec_name}, not carrier-level.",
                )
                + " (Confirm the apparent homozygosity is not a hemizygous deletion or "
                "uniparental-disomy artifact.)",
            )
        if sex == "female":
            if zyg == "het":
                return _band(
                    "xlinked_female_carrier", "warning",
                    "X-linked heterozygous female",
                    f"Heterozygous female for an X-linked {var_desc} in {gene_upper}. "
                    f"X-linked heterozygous females are not necessarily unaffected — "
                    + _c(
                        "skewed X-inactivation can cause manifestation",
                        "if reclassified pathogenic, skewed X-inactivation could cause "
                        "manifestation",
                    )
                    + " (e.g. GLA / Fabry, LAMP2 / Danon).",
                )
            return _band(
                "uncertain", "info", "X-linked — zygosity not provided",
                f"{gene_upper} is X-linked and the proband is female, but zygosity was "
                f"not provided, so heterozygous-carrier vs biallelic cannot be "
                f"determined.",
            )
        return _band(
            "uncertain", "warning", "X-linked — proband sex not provided",
            f"{gene_upper} is X-linked, but proband sex was not provided, so hemizygous "
            f"(male) vs heterozygous-carrier (female) cannot be determined.",
        )

    if homozygous:
        return _band(
            "biallelic_affected", "info", "Biallelic genotype (homozygous)",
            f"Homozygous {var_desc} in {gene_upper} — a biallelic genotype"
            + _c(
                f" for {rec_name}, not carrier-level.",
                f"; if reclassified pathogenic this would be consistent with affected "
                f"for {rec_name}, not carrier-level.",
            )
            + " (Confirm the apparent homozygosity is not a hemizygous deletion or "
            "uniparental-disomy artifact.)",
        )
    if pm3_biallelic:
        return _band(
            "biallelic_affected", "info", "Biallelic genotype (in trans)",
            f"A {var_desc} in trans with a curator-asserted second pathogenic allele in "
            f"{gene_upper} — a biallelic (compound-heterozygous) genotype"
            + _c(
                f" for {rec_name}, not carrier-level.",
                f"; if reclassified pathogenic this would be consistent with affected "
                f"for {rec_name}, not carrier-level.",
            )
            + " The second allele was not re-classified here.",
        )

    if zyg == "":
        if rec_present:
            return _band(
                "uncertain", "info", "Carrier vs affected — zygosity not provided",
                f"A {var_desc} in {gene_upper}, recessive for {rec_name}. Zygosity was "
                f"not provided, so carrier (heterozygous) vs biallelic cannot be "
                f"determined.",
            )
        return None

    if rec_present and dom_present:
        if scoped:
            return _band(
                "dual_risk", "warning", "Dual inheritance — recessive and dominant",
                f"For {rec_name}, {gene_upper} has both recessive and dominant disease. "
                + _c(
                    f"This heterozygous {var_desc} is carrier-level for the recessive "
                    f"form, but a single allele may itself act via a dominant mechanism "
                    f"— not a simple carrier finding.",
                    f"If reclassified pathogenic, this heterozygous {var_desc} could be "
                    f"carrier-level for the recessive form OR act via a dominant "
                    f"mechanism — not a simple carrier scenario.",
                ),
            )
        return _band(
            "dual_risk", "info", "Dual inheritance — recessive and dominant",
            f"{gene_upper} has recessive disease ({rec_name}) and dominant disease "
            f"({dom_name}). "
            + _c(
                f"A heterozygous {var_desc} is carrier-level for the recessive "
                f"condition, but a single allele may act via the dominant mechanism.",
                f"If reclassified pathogenic, a heterozygous {var_desc} could be "
                f"carrier-level for the recessive condition OR act via the dominant "
                f"mechanism.",
            )
            + disclose,
        )
    if rec_present:
        if scoped:
            return _band(
                "carrier", "info", "Heterozygous carrier (recessive gene)",
                f"Heterozygous {var_desc} in {gene_upper}, recessive for {rec_name}"
                + _c(
                    ": a single allele is carrier-level",
                    ". If reclassified pathogenic, a single allele would be carrier-level",
                )
                + " — biallelic variants are required for the recessive disease. A "
                "second allele (e.g. a CNV or deep-intronic variant) cannot be excluded "
                "from this single finding.",
            )
        return _band(
            "uncertain", "info", "Carrier vs affected — inheritance not scoped",
            f"{gene_upper} is recessive for {rec_name}." + disclose,
        )
    if dom_present:
        return None
    return _band(
        "uncertain", "info", "Carrier status — inheritance not established",
        f"A {var_desc} in {gene_upper}, but no recessive or dominant mode is established "
        f"in the curated high-confidence sources, so carrier vs affected cannot be "
        f"determined." + disclose,
    )


_KCNQ1_PVS1_TREE = (
    (581, "PVS1", "codons 1-581, where nonsense-mediated decay is predicted"),
    (620, "PVS1_Moderate",
     "codons 582-620, where NMD is not predicted but the variant removes the "
     "subunit assembly domain (SAD, residues 589-620) that mediates "
     "tetramerization, so the channel cannot assemble (PMID 10654932)"),
    (676, "PVS1_Supporting",
     "codons 621-676, where NMD is not predicted AND the SAD is retained, so "
     "these distal variants may still yield functional channels"),
)


def _pvs1_kcnq1_strength(vep: dict) -> tuple[str, str] | None:
    """KCNQ1 GN112's own PVS1 ladder for TRUNCATING variants, by codon band.

    Returns ``(tier, rationale)`` or None to fall through to the generic tree —
    when the variant is not truncating, when no protein position is available, or
    when the position is past codon 676 (the spec stops there and inventing a
    band beyond it would be our rule, not theirs).

    Verbatim thresholds from the harvested GN112 rows, so a re-harvest that
    changes them shows up in the spec file and in the test that reads it.
    """
    consequence = (vep.get("most_severe_consequence") or "").lower()
    if not any(c in consequence for c in _LOF_CONSEQUENCE_TYPES):
        return None
    pos = vep.get("protein_position")
    try:
        codon = int(str(pos).split("-", 1)[0].split("/", 1)[0])
    except (TypeError, ValueError):
        return None
    if codon < 1:
        return None
    for last, tier, why in _KCNQ1_PVS1_TREE:
        if codon <= last:
            return tier, (
                f"KCNQ1 GN112 applies PVS1 at "
                f"{tier.split('_', 1)[1] if '_' in tier else 'VeryStrong'} for a "
                f"truncating variant at codon {codon} — {why}"
            )
    return None


def _pvs1_splice_strength(vep: dict, ev: dict | None = None) -> tuple[str, str] | None:
    """ClinGen PVS1 splice-arm strength for a canonical ±1,2 splice variant.

    Returns ``(tier, rationale)`` where tier is ``"PVS1"`` (full / Very
    Strong, +8), ``"PVS1_Strong"`` (+4), or ``"PVS1_Moderate"`` (+2); or
    ``None`` to fall back to base PVS1 when transcript structure is
    unavailable. Only the deterministic structural signals are used — frame
    of the skipped exon, fraction of transcript removed, and NMD competence.
    The "critical region" promotion node is deferred (the >10% rule backstops
    large truncations) and the downgrade never depends on SpliceAI being
    present. See EREPO_PVS1_SPLICE_DESIGN.md (decisions Q1=defer, Q2=frame-
    only, Q3=keep +8).

    The mapping rests on VEP's transcript-rank exon/intron numbering: a
    ``splice_acceptor`` at intron N hits exon N+1; a ``splice_donor`` at
    intron N hits exon N. The change can only ever *lower* strength."""
    consequence = (vep.get("most_severe_consequence") or "").lower()
    exon_lengths = vep.get("exon_lengths")
    intron = vep.get("intron")
    if not exon_lengths or not intron or "/" not in str(intron):
        return None
    try:
        intron_n = int(str(intron).split("/", 1)[0])
    except ValueError:
        return None
    if "splice_acceptor" in consequence:
        affected = intron_n + 1
    elif "splice_donor" in consequence:
        affected = intron_n
    else:
        return None
    if affected < 1 or affected > len(exon_lengths):
        return None
    exon_len = exon_lengths[affected - 1]
    n_exons = len(exon_lengths)
    if sum(exon_lengths) <= 0:
        return None
    _aa_len = _as_float(((ev or {}).get("uniprot") or {}).get("length"))
    cds_nt = int(_aa_len) * 3 if _aa_len and _aa_len > 0 else None
    if cds_nt:
        total, basis = cds_nt, "of the protein"
    else:
        total, basis = sum(exon_lengths), (
            "of the transcript incl. UTRs — protein length unavailable, so this "
            "understates the protein fraction"
        )
    is_terminal_exon = affected == 1 or affected == n_exons
    frame_preserved = None if is_terminal_exon else (exon_len % 3 == 0)
    pct_removed = exon_len / total
    last_or_penult = affected >= n_exons - 1
    pct_txt = f"{pct_removed:.0%} {basis}"
    if frame_preserved is False and not last_or_penult:
        return ("PVS1", f"out-of-frame skip of exon {affected} ({pct_txt}) — "
                "frameshift PTC predicted NMD-competent")
    if frame_preserved is None:
        region = (
            f"skip of terminal exon {affected} (frame effect not determinable — "
            "the exon's genomic length includes UTR)"
        )
    else:
        region = ("in-frame skip" if frame_preserved
                  else "out-of-frame skip escaping NMD")
    if is_terminal_exon and pct_removed > 0.10:
        return ("PVS1_Moderate",
                f"{region} nominally removes {pct_txt}, but part of a terminal "
                "exon is UTR so >10% of the protein is not established — "
                "PVS1_Moderate (conservative)")
    if pct_removed > 0.10:
        return ("PVS1_Strong",
                f"{region} of exon {affected} removes {pct_txt} (>10%)")
    return ("PVS1_Moderate",
            f"{region} of exon {affected} removes {pct_txt} (≤10%)")


def _pvs1_splice_cryptic_strength(vep: dict, ev: dict) -> tuple[str, str] | None:
    """SpliceAI fallback for the PVS1 splice arm when the exon-skip model
    (``_pvs1_splice_strength``) has no transcript structure — which is the case
    for EVERY RefSeq (NM_) input, since Ensembl ``/lookup/id`` rejects RefSeq
    accessions so ``vep['exon_lengths']`` is None. Without this fallback every
    canonical ±1,2 splice variant entered on a RefSeq transcript over-fires
    full PVS1 (+8) regardless of its true consequence.

    Downgrades to ``PVS1_Moderate`` ONLY when SpliceAI confidently predicts the
    canonical site is LOST and a same-type cryptic site is GAINED an in-frame
    distance away — the Abou-Tayoun (2018) node where a cryptic-rescue produces
    a small in-frame indel that preserves the reading frame and is not
    NMD-competent. Conservative gates: confident same-type loss+gain (Δ ≥ 0.5),
    frame-preserving shift (divisible by 3) and small (≤ 30 nt). Returns None
    (→ keep base PVS1) otherwise — it can only ever LOWER strength, never raise.
    """
    sa = (ev or {}).get("spliceai") or {}
    if not sa.get("ok"):
        return None
    per_tx = sa.get("scores_per_transcript") or []
    if not per_tx or not isinstance(per_tx[0], dict):
        return None
    s = per_tx[0]
    consequence = (vep.get("most_severe_consequence") or "").lower()
    if "splice_acceptor" in consequence:
        ds_loss, dp_loss, ds_gain, dp_gain, site = (
            s.get("DS_AL"), s.get("DP_AL"), s.get("DS_AG"), s.get("DP_AG"), "acceptor")
    elif "splice_donor" in consequence:
        ds_loss, dp_loss, ds_gain, dp_gain, site = (
            s.get("DS_DL"), s.get("DP_DL"), s.get("DS_DG"), s.get("DP_DG"), "donor")
    else:
        return None
    if None in (ds_loss, dp_loss, ds_gain, dp_gain):
        return None
    if ds_loss < 0.5 or ds_gain < 0.5:
        return None
    try:
        shift = abs(int(dp_gain) - int(dp_loss))
    except (TypeError, ValueError):
        return None
    if shift > 0 and shift % 3 == 0 and shift <= 30:
        return ("PVS1_Moderate",
                f"canonical {site} site lost (SpliceAI Δ={ds_loss:.2f}) with an "
                f"in-frame cryptic {site} {shift} nt away (Δ={ds_gain:.2f}; "
                f"{shift // 3}-codon in-frame change) — reading frame preserved, "
                "PVS1_Moderate per the Abou-Tayoun splice decision tree")
    return None


def _eval_ba1(gnomad_af: float | None, ba1_threshold: float,
              ba1_basis: str, gnomad_an: int | None = None,
              ba1_exception: str | None = None,
              af_basis: str = "popmax FAF95") -> dict:
    """── BA1 — Benign Stand-alone ─────────────────────────────────────
    ACMG/SVI fix BA1 at 5%, but the cardiac/RASopathy VCEPs publish far
    tighter gene-specific BA1 cutoffs (RASopathy 0.0005, HCM/Marfan 0.001,
    KCNQ1 0.004) under the max-credible-AF framework — same per-gene wiring
    already in force for PM2/BS1. Out-of-table genes hold the 0.05 default
    (see _ba1_threshold + its gnomAD-version safeguard).

    Implements the SVI's UPDATED BA1 definition (Ghosh et al., Genet Med 2018,
    PMID 30311383): "Allele frequency is >0.05 in any general continental
    population dataset OF AT LEAST 2,000 OBSERVED ALLELES and found in a gene
    without a gene- OR VARIANT-specific BA1 modification." Two added gates:

      * >= 2,000 observed alleles (``gnomad_an``). A high AF derived from a
        handful of alleles is noise, and gnomAD AN genuinely drops to the
        hundreds at poorly covered sites. BA1 forces Benign AND zeroes every
        pathogenic criterion, so firing it off a low-AN site is the worst
        false-benign available. Fail closed when AN is unknown.
      * the nine-variant EXCEPTION LIST (``ba1_exception``, see
        _BA1_EXCEPTION_VARIANTS) — variants with MAF >5% that ClinGen holds
        Pathogenic/VUS anyway. Found in practice: HFE c.845G>A (C282Y) has a
        real gnomAD v4 FAF95 popmax of 0.071 and was returning Benign."""
    if ba1_exception is not None:
        return _hard_coded_entry(
            "BA1", "not_met", None,
            f"On the ClinGen SVI BA1 exception list ({ba1_exception}) — one of "
            "the nine variants excluded from BA1 despite population MAF >5% "
            "because there is evidence of pathogenicity (Ghosh 2018, SVI "
            "updated BA1 recommendation). Stand-alone benign is NOT applied; "
            "assess on the remaining evidence"
            + (f". gnomAD {af_basis}={gnomad_af:.6f}"
               if gnomad_af is not None else ""),
        )
    if (
        gnomad_af is not None
        and gnomad_af >= ba1_threshold
        and (gnomad_an is None or gnomad_an < 2000)
    ):
        return _hard_coded_entry(
            "BA1", "not_met", None,
            f"gnomAD {af_basis}={gnomad_af:.6f} exceeds the BA1 threshold "
            f"{ba1_threshold} ({ba1_basis}), but only "
            f"{'an unknown number of' if gnomad_an is None else f'{gnomad_an:,}'} "
            "observed alleles — the SVI updated BA1 definition requires a "
            "population dataset of at least 2,000 observed alleles "
            "(Ghosh 2018), so stand-alone benign is withheld (fail-closed)",
        )
    if gnomad_af is not None and gnomad_af >= ba1_threshold:
        return _hard_coded_entry(
            "BA1", "met", "BA1",
            f"gnomAD {af_basis}={gnomad_af:.6f} is at or above the stand-alone "
            f"benign threshold {ba1_threshold} ({ba1_basis})",
        )
    if gnomad_af is None:
        ba1_ev = "gnomAD data unavailable — BA1 cannot be evaluated"
    else:
        ba1_ev = (
            f"gnomAD {af_basis}={gnomad_af:.6f} at/below stand-alone "
            f"benign threshold {ba1_threshold} ({ba1_basis})"
        )
    return _hard_coded_entry("BA1", "not_met", None, ba1_ev)


def _eval_bs1(gnomad_af: float | None, bs1_threshold: float,
              ba1_threshold: float, inheritance: str, freq_src: str,
              ba1_met: bool = False,
              af_basis: str = "popmax FAF95") -> dict:
    """── BS1 — Benign Strong (AF greater than expected) ───────────────
    Richards 2015 (PMID 25741868): BS1 = "allele frequency greater than
    expected for the disorder"; BA1 = >5% stand-alone benign. BS1 has no
    fixed numeric value in the guideline — it is disease-specific (max-
    credible-AF framework, Whiffin/Ware 2017). ClinGen SVI likewise fixes
    only BA1 at 5% and defers BS1 to gene/disease-specific thresholds.
    AR disorders tolerate higher population AF than AD (carriers are common,
    affecteds need biallelic variants), so the AR threshold sits above the
    AD one — but it MUST be below BA1's 0.05, otherwise the gate window
    `(threshold, 0.05]` below is empty and BS1 can never fire for AR. The
    previous AR value of 0.05 was exactly that empty-window bug. 0.005
    (0.5%, an order of magnitude under BA1) is a provisional generic AR
    threshold pending gene-specific VCEP calibration.

    The generic AD/unknown fallback (see _bs1_ad_threshold) is disease-
    prevalence-aware: 0.0004 for Definitive/Strong ClinGen genes, else 0.001.
    Gene-specific VCEP BS1 floors now take precedence when available (see
    vcep_frequency_thresholds.json): RASopathy 0.00025, cardiomyopathy
    0.0001 (MYBPC3 0.0002), KCNQ1 0.0004 — all on FAF95 popmax. These recover
    the eRepo BS1 calls the old 0.0004 floor missed (e.g. HRAS c.510G>A,
    FAF95 ~0.000252 >= 0.00025). RESIDUAL: MAP2K2 c.813C>T (FAF95 ~4.4e-5)
    sits below even the RASopathy floor — the curator applied a gene-specific
    MAP2K2 rule below the VCEP-wide value; not captured generically.

    BS1 fires in the band (bs1_threshold, ba1_threshold]: above expected for
    the disorder but below the stand-alone benign cutoff (anything above
    ba1_threshold is covered by BA1, not double-counted as BS1). The upper
    bound is now the gene-specific BA1 (was a flat 0.05) so the window stays
    consistent with the tightened per-gene BA1 above."""
    if (
        gnomad_af is not None
        and gnomad_af >= bs1_threshold
        and gnomad_af < ba1_threshold
    ):
        return _hard_coded_entry(
            "BS1", "met", "BS1_Strong",
            f"gnomAD {af_basis}={gnomad_af:.5f} is at or above the expected threshold "
            f"for {inheritance or 'stated'} inheritance ({bs1_threshold}); "
            f"below BA1 stand-alone cutoff ({ba1_threshold})",
        )
    if gnomad_af is None:
        bs1_ev = "gnomAD data unavailable — BS1 cannot be evaluated"
    elif gnomad_af >= ba1_threshold and ba1_met:
        bs1_ev = (
            f"gnomAD {af_basis}={gnomad_af:.6f} exceeds BA1 stand-alone cutoff "
            f"({ba1_threshold}) — this frequency is already counted under BA1, "
            "which fired; counting it again under BS1 would double-count one "
            "observation"
        )
    elif gnomad_af >= ba1_threshold:
        bs1_ev = (
            f"gnomAD AF={gnomad_af:.6f} exceeds the BA1 cutoff "
            f"({ba1_threshold}) but BA1 was WITHHELD (see the BA1 entry for the "
            "reason). BS1 is the same frequency argument at lower strength, so "
            "it is withheld too rather than substituting for the criterion that "
            "was just declined — no frequency-based benign evidence is applied "
            "to this variant"
        )
    else:
        bs1_ev = (
            f"gnomAD AF={gnomad_af:.6f} below expected threshold "
            f"({bs1_threshold}; {freq_src})"
        )
    return _hard_coded_entry("BS1", "not_met", None, bs1_ev)


def _eval_pm2(gnomad_af: float | None, pm2_threshold: float, pm2_src: str,
              ba1_met: bool, bs1_met: bool, gnomad_confirmed_absent: bool,
              gn: dict) -> dict:
    """── PM2 — Absent in controls ─────────────────────────────────────
    Gene-specific VCEP ceiling when available; pm2_threshold == 0 encodes the
    RASopathy absence-only rule (PM2 only when the variant is absent from
    gnomAD). Generic fallback (AR 0.01 / AD 1e-4) for uncovered genes."""
    if ba1_met or bs1_met:
        return _hard_coded_entry(
            "PM2", "not_met", None, "Precluded by BA1/BS1",
        )
    if gnomad_af is None and not gnomad_confirmed_absent:
        _reason = (
            "gnomAD frequency unresolved for this indel (no matching rsID) — "
            "PM2 not applied (frequency unknown, not confirmed absent)"
            if gn.get("indel_unresolved")
            else "gnomAD lookup unavailable — PM2 cannot be evaluated"
        )
        return _hard_coded_entry("PM2", "not_met", None, _reason)
    if gnomad_af is None or gnomad_af == 0:
        return _hard_coded_entry(
            "PM2", "met", "PM2_Supporting",
            f"Absent from gnomAD — PM2_Supporting ({pm2_src}; "
            "ClinGen SVI 2020 default strength)",
        )
    if pm2_threshold > 0 and gnomad_af <= pm2_threshold:
        return _hard_coded_entry(
            "PM2", "met", "PM2_Supporting",
            f"gnomAD FAF={gnomad_af:.6f} at/below PM2 ceiling {pm2_threshold} "
            f"({pm2_src}) — PM2_Supporting per ClinGen SVI 2020 default strength",
        )
    if pm2_threshold == 0:
        return _hard_coded_entry(
            "PM2", "not_met", None,
            f"gnomAD FAF={gnomad_af:.6f} — present in gnomAD; this gene's VCEP "
            "applies PM2 only when absent (absence-only rule)",
        )
    return _hard_coded_entry(
        "PM2", "not_met", None,
        f"gnomAD FAF={gnomad_af:.5f} exceeds PM2 ceiling "
        f"({pm2_threshold}, {pm2_src})",
    )


_REVEL_PP3 = ((0.932, "PP3_Strong"), (0.773, "PP3_Moderate"), (0.644, "PP3_Supporting"))
_REVEL_BP4 = (
    (0.003, "BP4_VeryStrong"), (0.016, "BP4_Strong"),
    (0.183, "BP4_Moderate"), (0.290, "BP4_Supporting"),
)


def _calibrated_revel_call(revel) -> tuple[str, str] | None:
    """Map a REVEL score to ``("PP3"|"BP4", criteria_strength)`` per ClinGen SVI
    (Pejaver 2022). Returns ``None`` when REVEL is missing/unparseable or falls
    in the indeterminate band (0.290 < REVEL < 0.644)."""
    if revel is None:
        return None
    try:
        r = float(revel)
    except (TypeError, ValueError):
        return None
    for thr, strength in _REVEL_PP3:
        if r >= thr:
            return ("PP3", strength)
    for thr, strength in _REVEL_BP4:
        if r <= thr:
            return ("BP4", strength)
    return None


_STRENGTH_RANK = {"Supporting": 1, "Moderate": 2, "Strong": 3, "VeryStrong": 4}


def _insilico_strength_ceiling(gene: str | None, ev: dict) -> str | None:
    """Max strength tier the PP3/BP4 REVEL ladder may reach. ACGS-2024 engine-
    wide posture: ``"Supporting"`` for every gene (see module comment). The
    ``gene``/``ev`` params are retained for a future per-VCEP re-enablement."""
    return "Supporting"


def _clamp_insilico_strength(strength: str, ceiling: str | None) -> str:
    """Lower a ``"PP3_Strong"``-style strength to ``ceiling`` (never raise).
    ``ceiling=None`` is a no-op."""
    if not ceiling:
        return strength
    code, _, tier = strength.partition("_")
    if _STRENGTH_RANK.get(tier, 1) > _STRENGTH_RANK.get(ceiling, 4):
        return f"{code}_{ceiling}"
    return strength


_TAVTIGIAN_POINTS = {"Supporting": 1, "Moderate": 2, "Strong": 4, "VeryStrong": 8}


def _strength_points(strength: str | None) -> int:
    """Magnitude (unsigned) of a ``"PP3_Strong"``-style strength, 0 if absent."""
    if not strength:
        return 0
    return _TAVTIGIAN_POINTS.get(strength.partition("_")[2], 0)


_SPEC_STRENGTH_RANK = {
    "Supporting": 1, "Moderate": 2, "Strong": 4,
    "Very Strong": 8, "Stand Alone": 8,
    "VeryStrong": 8, "StandAlone": 8,
}
_RANK_TO_SUFFIX = {1: "Supporting", 2: "Moderate", 4: "Strong", 8: "VeryStrong"}


def spec_strength_ceiling(gene: str | None, code: str) -> str | None:
    """Highest strength this gene's VCEP publishes for ``code``, e.g. "Moderate".

    None when the gene has no spec, the criterion is not applicable, or no
    strength row is published — in which case there is no ceiling to impose and
    the caller must not invent one.
    """
    entry = _crit_spec(gene, code)
    if entry.get("applicability") != "applicable":
        return None
    ranks = [_SPEC_STRENGTH_RANK[s] for s in (entry.get("strengths") or {})
             if s in _SPEC_STRENGTH_RANK]
    return _RANK_TO_SUFFIX.get(max(ranks)) if ranks else None


def cap_all_to_spec_strength(criteria: list[dict], gene: str | None) -> list[dict]:
    """Apply cap_to_spec_strength across a whole criteria list.

    ⚠ MUST RUN POST-MERGE AND POST-GATE, not only inside
    compute_hard_coded_criteria. The two real over-calls this exists for are
    invisible there:
      * PM1 on TNNT2 (spec publishes Supporting; we emit Moderate) is set by
        `_gate_criteria_applicability`, which runs AFTER the deterministic set is
        built — the same ordering trap that hid the PM1/PM5 collision;
      * BS3 on RRAS2 (spec publishes Supporting; base is Strong) is
        AI-evaluated, so it does not exist until the AI set is merged.
    Calling it only at build time therefore fixed neither. Only ever lowers, so
    running it at both points is safe and idempotent.
    """
    for entry in criteria:
        cap_to_spec_strength(entry, gene)
    return criteria


def cap_to_spec_strength(entry: dict, gene: str | None) -> dict:
    """Lower ``entry``'s strength to what its gene's VCEP actually publishes.

    Only ever LOWERS, and only for a `met` criterion. Mutates and returns the
    entry so it can be used inline where criteria are built.
    """
    if not entry or entry.get("status") != "met":
        return entry
    code = entry.get("code") or ""
    ceiling = spec_strength_ceiling(gene, code)
    if not ceiling:
        return entry
    cur = entry.get("criteria_strength") or code
    cur_pts = _strength_points(cur) or abs(_BARE_CODE_POINTS.get(code, 0))
    cap_pts = _SPEC_STRENGTH_RANK[ceiling]
    if cur_pts <= cap_pts:
        return entry
    _gn = (_VCEP_CRIT.get((gene or "").upper()) or {}).get("gn", "CSpec")
    entry["criteria_strength"] = f"{code}_{ceiling}"
    entry["evidence"] = (
        f"{gene} VCEP ({_gn}) publishes {code} at {ceiling} only, so it is "
        f"applied at {ceiling} rather than the ACMG 2015 base strength. "
        f"(was: {str(entry.get('evidence') or '')[:160]})"
    )
    return entry


def _eval_pp3_bp4(
    ev: dict, vep: dict, consequence: str, spliceai_max: float,
    pp3_ceiling: str | None = None,
) -> tuple[list[dict], bool, bool]:
    """── PP3 / BP4 — In-silico evidence (mutually exclusive) ──────────
    Returns ``([PP3_entry, BP4_entry], noncoding_splice, splice_clear)`` —
    the two boolean splice flags are computed here (PP3/BP4 + BP7 share them)
    and handed back so BP7 stays on every path without recomputation.

    For missense variants this is now REVEL-PRIMARY: a REVEL score in a
    calibrated pathogenic/benign band (see ``_calibrated_revel_call``) drives
    PP3/BP4 at the matching ClinGen SVI strength (Supporting→Strong for PP3,
    Supporting→VeryStrong for BP4). When REVEL is absent (e.g. indels,
    non-coding) or indeterminate, the function falls back to the legacy
    ≥2-tool CADD/AlphaMissense/SpliceAI consensus (Supporting only) and the
    non-coding/splice arms below — preserving the prior behaviour exactly."""
    revel = vep.get("revel_score")
    cadd = vep.get("cadd_phred")
    am_block = ev.get("alphamissense") or {}
    am_score = am_block.get("am_pathogenicity")
    if am_score is None:
        am_score = am_block.get("score")
    scores = {
        "revel": revel,
        "cadd": cadd,
        "alphamissense": am_score,
        "spliceai_max": spliceai_max,
    }
    thresholds_path = {
        "revel": 0.70, "cadd": 25, "alphamissense": 0.564, "spliceai_max": 0.5,
    }
    thresholds_benign = {
        "revel": 0.30, "cadd": 15, "alphamissense": 0.340,
    }
    tools_path: list[str] = []
    tools_benign: list[str] = []
    for k, v in scores.items():
        if v is None:
            continue
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if f >= thresholds_path.get(k, float("inf")):
            tools_path.append(k)
        if k in thresholds_benign and f <= thresholds_benign[k]:
            tools_benign.append(k)
    sa_ok = bool((ev.get("spliceai") or {}).get("ok"))
    noncoding_splice = (
        any(t in consequence for t in (
            "intron_variant", "splice_region", "splice_polypyrimidine",
            "non_coding_transcript", "3_prime_utr_variant", "5_prime_utr_variant",
            "splice_donor_region", "splice_acceptor_region",
        ))
        and "synonymous" not in consequence
        and not any(t in consequence for t in (
            "missense", "stop_gain", "stop_lost", "frameshift",
            "start_lost", "inframe", "splice_donor_variant",
            "splice_acceptor_variant",
        ))
    )
    is_pure_utr = (
        any(t in consequence for t in ("3_prime_utr_variant", "5_prime_utr_variant"))
        and not any(t in consequence for t in ("splice", "intron"))
    )
    splice_clear = (sa_ok and spliceai_max < 0.1) or (is_pure_utr and not sa_ok)
    out: list[dict] = []
    revel_call = _calibrated_revel_call(revel) if "missense" in consequence else None
    if revel_call is not None and not (
        revel_call[0] == "BP4" and sa_ok and spliceai_max >= 0.2
    ):
        code, strength = revel_call
        rv = float(revel)
        cite = "ClinGen SVI calibrated REVEL threshold (Pejaver 2022)"
        capped = _clamp_insilico_strength(strength, pp3_ceiling)
        cap_note = (
            f"; capped to {capped.replace('_', ' ')} — in-silico evidence not "
            "escalated above Supporting (ClinGen cardiac VCEPs / ACGS 2024)"
            if capped != strength else ""
        )
        strength = capped
        if code == "PP3":
            out.append(_hard_coded_entry(
                "PP3", "met", strength,
                f"REVEL {rv:.3f} → {strength.replace('_', ' ')} — {cite}{cap_note}",
            ))
            out.append(_hard_coded_entry(
                "BP4", "not_met", None,
                f"REVEL {rv:.3f} supports pathogenicity (PP3), not BP4",
            ))
        else:
            out.append(_hard_coded_entry(
                "PP3", "not_met", None,
                f"REVEL {rv:.3f} supports a benign call (BP4), not PP3",
            ))
            out.append(_hard_coded_entry(
                "BP4", "met", strength,
                f"REVEL {rv:.3f} → {strength.replace('_', ' ')} — {cite}{cap_note}",
            ))
        return out, noncoding_splice, splice_clear
    if len(tools_path) >= 2 and len(tools_path) > len(tools_benign):
        out.append(_hard_coded_entry(
            "PP3", "met", "PP3_Supporting",
            f"{len(tools_path)} tools above pathogenic threshold: "
            f"{', '.join(tools_path)}",
        ))
        out.append(_hard_coded_entry(
            "BP4", "not_met", None,
            f"PP3 outweighs BP4 ({len(tools_path)} pathogenic vs "
            f"{len(tools_benign)} benign tools)",
        ))
    elif len(tools_benign) >= 2 and len(tools_benign) > len(tools_path):
        out.append(_hard_coded_entry(
            "PP3", "not_met", None,
            f"BP4 outweighs PP3 ({len(tools_benign)} benign vs "
            f"{len(tools_path)} pathogenic tools)",
        ))
        out.append(_hard_coded_entry(
            "BP4", "met", "BP4_Supporting",
            f"≥2 tools below benign threshold: {', '.join(tools_benign)}",
        ))
    else:
        if len(tools_path) == 0 and len(tools_benign) == 0:
            pp3_ev = "No in-silico tools above pathogenic threshold"
            bp4_ev = "No in-silico tools below benign threshold"
        elif len(tools_path) == len(tools_benign):
            pp3_ev = bp4_ev = (
                f"Conflicting computational evidence — {len(tools_path)} "
                f"pathogenic vs {len(tools_benign)} benign tools"
            )
        elif len(tools_path) == 1:
            pp3_ev = (
                f"Only 1 tool above pathogenic threshold (requires ≥2): "
                f"{tools_path}"
            )
            bp4_ev = "Insufficient benign tools (requires ≥2)"
        else:
            pp3_ev = "Insufficient pathogenic tools (requires ≥2)"
            bp4_ev = (
                f"Only 1 tool below benign threshold (requires ≥2): "
                f"{tools_benign}"
            )
        if (
            len(tools_path) == 0
            and len(tools_benign) < 2
            and noncoding_splice
            and splice_clear
        ):
            out.append(_hard_coded_entry("PP3", "not_met", None, pp3_ev))
            sa_note = (
                f"SpliceAI max delta {spliceai_max:.3f} (<0.1)" if sa_ok
                else "in UTR with no splice context (SpliceAI not applicable)"
            )
            out.append(_hard_coded_entry(
                "BP4", "met", "BP4_Supporting",
                f"Non-coding variant ({consequence}) — {sa_note}; no predicted "
                "splice impact; missense in-silico tools not applicable",
            ))
        elif (
            len(tools_benign) == 0
            and noncoding_splice
            and sa_ok
            and spliceai_max >= _SPLICEAI_PP3_DELTA
        ):
            out.append(_hard_coded_entry(
                "PP3", "met", "PP3_Supporting",
                f"Non-coding/splice-region variant ({consequence}) with SpliceAI "
                f"max delta {spliceai_max:.3f} (≥{_SPLICEAI_PP3_DELTA}) — predicted "
                "splice impact; "
                "missense in-silico tools not applicable",
            ))
            out.append(_hard_coded_entry("BP4", "not_met", None, bp4_ev))
        else:
            out.append(_hard_coded_entry("PP3", "not_met", None, pp3_ev))
            out.append(_hard_coded_entry("BP4", "not_met", None, bp4_ev))
    return out, noncoding_splice, splice_clear


_SPLICEAI_PP3_DELTA = 0.2


def _eval_bp7(consequence: str, spliceai_max: float, sa_ok: bool,
              noncoding_splice: bool, splice_clear: bool,
              phylop: float | None = None) -> dict:
    """── BP7 — Synonymous / non-coding, no predicted splice impact ─────
    ClinGen SVI extends BP7 beyond synonymous to intronic / UTR / splice-
    region variants with no predicted splice impact. Require SpliceAI to be
    PRESENT (sa_ok) and confidently low — never fire on missing SpliceAI
    (spliceai_max defaults to 0.0 when unavailable, which would otherwise
    mis-fire BP7 on absent data). `noncoding_splice` (computed above) already
    excludes coding/canonical-splice consequences and synonymous.

    Per ACMG/AMP, BP7 ALSO requires the nucleotide to be NOT highly conserved.
    ACMG defines no numeric cutoff; we use phyloP100way > 2.0 as "highly
    conserved", matching the ClinGen Myeloid Malignancy VCEP operationalization
    (RUNX1 rules, Wu et al. Blood Adv 2022 — phyloP100way ≤ 2.0 counts as not
    conserved, chosen for >85% concordance with known-benign variants). A
    positively-conserved variant (phyloP > 2.0) is BLOCKED from BP7.

    FAIL-CLOSED on missing conservation: if phyloP is UNAVAILABLE (fetch failed,
    pyBigWig/bigWig absent, or offline-strict without a local track), BP7 is
    WITHHELD (conservative not_met) rather than fired — we cannot confirm the
    nucleotide is not highly conserved, and firing BP7 ungated risks a
    false-benign call. (⚠ This changes BP7 firing on variants whose phyloP
    could not be read — re-validate before relying on it.)"""
    is_synonymous = "synonymous" in consequence
    bp7_eligible = is_synonymous or noncoding_splice
    bp7_clear = (sa_ok and spliceai_max < 0.1) if is_synonymous else splice_clear
    kind = "Synonymous" if is_synonymous else "Non-coding/intronic/UTR"
    conservation_known = phylop is not None
    is_conserved = conservation_known and float(phylop) > 2.0
    if bp7_eligible and bp7_clear and conservation_known and not is_conserved:
        sa_note = (
            f"SpliceAI max delta {spliceai_max:.3f} (<0.1)" if sa_ok
            else "in UTR with no splice context (SpliceAI not applicable)"
        )
        cons_note = f"phyloP100way {float(phylop):.2f} (≤2.0, not highly conserved)"
        return _hard_coded_entry(
            "BP7", "met", "BP7_Supporting",
            f"{kind} variant ({consequence}) — {sa_note}; {cons_note}; "
            "no predicted splice impact",
        )
    if bp7_eligible and bp7_clear and is_conserved:
        return _hard_coded_entry(
            "BP7", "not_met", None,
            f"{kind} variant ({consequence}) with no predicted splice impact, "
            f"but phyloP100way {float(phylop):.2f} (>2.0) — nucleotide is "
            "highly conserved, so BP7 (not highly conserved) does not apply",
        )
    if bp7_eligible and bp7_clear and not conservation_known:
        return _hard_coded_entry(
            "BP7", "not_met", None,
            f"{kind} variant ({consequence}) with no predicted splice impact, but "
            "phyloP100way conservation could not be determined — BP7 requires the "
            "nucleotide to be NOT highly conserved, so it is withheld (fail-closed)",
        )
    if not bp7_eligible:
        return _hard_coded_entry(
            "BP7", "not_met", None,
            f"Not a synonymous/non-coding variant (consequence: "
            f"{consequence or 'unknown'})",
        )
    if not sa_ok:
        return _hard_coded_entry(
            "BP7", "not_met", None,
            "SpliceAI unavailable for a splice-region/intronic variant — cannot "
            "confirm no splice impact for BP7 (conservative not_met)",
        )
    return _hard_coded_entry(
        "BP7", "not_met", None,
        f"SpliceAI max delta {spliceai_max:.3f} ≥0.1 — possible splice impact",
    )


def _pvs1_nmd_escape_strength(vep: dict, ev: dict) -> tuple[str, str]:
    """Abou-Tayoun / ClinGen SVI PVS1 strength for a nonsense/frameshift variant
    predicted to ESCAPE NMD (PTC in the last exon or last ~50 bp of the
    penultimate exon). Per the SVI decision tree such a truncation is
    PVS1_Strong when it removes a functionally critical region OR >10% of the
    protein, and PVS1_Moderate otherwise. "Critical region" is not
    machine-assessable here, so we apply the >10%-protein-removed backstop
    (identical rule to the splice arm's ``_pvs1_splice_strength``). Returns
    ``(tier, rationale)``; can only LOWER strength (Strong or Moderate, never
    the full +8). Falls back to PVS1_Strong when the removed fraction cannot be
    quantified (conservative — preserves the prior behaviour)."""
    pos = vep.get("protein_start") or _hgvsp_int_position(vep.get("hgvsp") or "")
    total = None
    up = ev.get("uniprot") or {}
    if _as_float(up.get("length")):
        total = float(up["length"])
    if not total:
        el = vep.get("exon_lengths")
        if el and sum(el) > 0:
            total = sum(el) / 3.0
    try:
        pos = int(pos)
    except (TypeError, ValueError):
        pos = None
    if not pos or not total or pos > total:
        return ("PVS1_Strong",
                "predicted to escape NMD; fraction of protein removed could not "
                "be quantified — PVS1_Strong (conservative)")
    frac_removed = (total - pos) / total
    if frac_removed > 0.10:
        return ("PVS1_Strong",
                f"predicted to escape NMD; truncation at residue {pos}/{int(total)} "
                f"removes {frac_removed:.1%} of the protein (>10%) — PVS1_Strong")
    return ("PVS1_Moderate",
            f"predicted to escape NMD; truncation at residue {pos}/{int(total)} "
            f"removes only {frac_removed:.1%} of the protein (≤10%, C-terminal) — "
            "PVS1_Moderate per the Abou-Tayoun/ClinGen SVI PVS1 decision tree")


def _next_inframe_met(sequence: str | None) -> int | None:
    """1-based residue of the next in-frame Methionine AFTER the initiator.

    Every downstream in-frame ATG translates to an M in the protein sequence,
    so "next putative in-frame start codon" is just the first M at residue > 1.
    Returns None when no sequence is available or there is no downstream M.
    """
    seq = (sequence or "").strip().upper()
    if len(seq) < 2:
        return None
    idx = seq.find("M", 1)
    return idx + 1 if idx >= 0 else None


def _pvs1_start_loss_strength(ev: dict) -> tuple[str | None, str]:
    """ClinGen SVI / Abou Tayoun 2018 initiation-codon arm.

    Verbatim from the recommendation: the SVI "generally does not recommend
    assigning PVS1 or PVS1_Strong for start loss variants". Then a three-way:

      1. an alternative functional transcript uses an alternative start codon
         -> do NOT apply PVS1 at any strength level;
      2. else, >=1 pathogenic variant reported 5' of the next downstream
         putative in-frame start codon (Methionine) -> PVS1_Moderate;
      3. else -> PVS1_Supporting.

    Returns ``(tier, rationale)``; tier None means not_met (branch 1).

    SUPPORTING IS THE DEFAULT AND MODERATE IS THE UPGRADE, deliberately.
    Moderate is the branch that requires positive evidence — a *reported*
    pathogenic variant in the window between the lost ATG and the next
    in-frame Met. So an absent protein sequence or a failed ClinVar lookup
    must land on Supporting (the paper's else-branch), never on Moderate.
    Getting this the other way round is how "we could not check" turns into
    an extra point of pathogenic evidence.

    Branch 1 is NOT implemented: deciding that an alternative *functional*
    transcript uses a different start codon needs per-transcript translation
    starts across the gene's transcript set. Its absence can only ever leave
    PVS1 too HIGH (Supporting/Moderate where the paper says none), so it is
    disclosed in the evidence string rather than silently assumed away.
    """
    payload = (
        ev.get("pvs1_start_loss")
        or ((ev.get("domain_plp") or {}).get("pvs1_start_loss"))
        or {}
    )
    disclose = (
        " Not assessed: whether an alternative functional transcript uses a "
        "different start codon, which would withdraw PVS1 entirely "
        "(Abou Tayoun 2018 initiation-codon branch 1)."
    )
    if payload.get("alt_start_transcript") is True:
        return (None, (
            "an alternative functional transcript uses a different start codon, "
            "so the ClinGen SVI / Abou Tayoun 2018 decision tree withholds PVS1 "
            "at every strength level for an initiation-codon variant"
        ))
    met = payload.get("next_inframe_met")
    upstream = payload.get("upstream_plp")
    if met and upstream:
        examples = ", ".join(
            str(v.get("name")) for v in (payload.get("top_variants") or [])[:2]
            if v.get("name")
        )
        return ("PVS1_Moderate", (
            f"{upstream} pathogenic/likely-pathogenic ClinVar record(s) lie 5' of "
            f"the next in-frame Met at residue {met}"
            + (f" (e.g. {examples})" if examples else "")
            + " — PVS1_Moderate per the ClinGen SVI / Abou Tayoun 2018 "
            "initiation-codon arm." + disclose
        ))
    if met:
        return ("PVS1_Supporting", (
            f"no pathogenic/likely-pathogenic ClinVar record lies 5' of the next "
            f"in-frame Met at residue {met}, so the Moderate branch is not met — "
            "PVS1_Supporting per the ClinGen SVI / Abou Tayoun 2018 "
            "initiation-codon arm." + disclose
        ))
    return ("PVS1_Supporting", (
        "PVS1_Supporting per the ClinGen SVI / Abou Tayoun 2018 initiation-codon "
        "arm. The Moderate branch requires a pathogenic variant reported 5' of "
        "the next in-frame Met, and that could not be established here ("
        + str(payload.get("reason") or "no protein sequence available")
        + "), so the paper's default applies rather than the upgrade." + disclose
    ))


def _eval_pvs1(ev: dict, vep: dict, consequence: str,
               spliceai_max: float) -> dict:
    """── PVS1 — Null variant in established LoF gene ──────────────────"""
    gene_lof = _gene_lof_mechanism(
        ev,
        gene_symbol=vep.get("gene_symbol") or "",
        consequence=consequence,
    )
    nmd_escape = bool(vep.get("nmd_escape"))
    is_lof = any(c in consequence for c in _LOF_CONSEQUENCE_TYPES)
    if not is_lof and gene_lof and spliceai_max >= 0.2:
        _entry = _hard_coded_entry(
            "PVS1", "not_met", None,
            f"Non-canonical splice-region/near-splice variant ({consequence}) "
            f"with SpliceAI max delta {spliceai_max:.3f} in an established LoF "
            "gene — a SpliceAI prediction alone does not meet the ClinGen PVS1 "
            "splice criteria (canonical ±1,2 or RNA-confirmed only); routed to "
            "PP3 as predictive evidence",
        )
        _entry["_pvs1_rna_routable"] = True
        return _entry
    if is_lof and gene_lof and (vep.get("gene_symbol") or "").upper() == "KCNQ1":
        _k = _pvs1_kcnq1_strength(vep)
        if _k:
            return _hard_coded_entry("PVS1", "met", _k[0], _k[1])

    _other_lof = ("stop_gained", "frameshift", "splice_acceptor",
                  "splice_donor", "transcript_ablation", "exon_loss")
    if (
        is_lof and gene_lof
        and "start_lost" in consequence
        and not any(c in consequence for c in _other_lof)
    ):
        tier, why = _pvs1_start_loss_strength(ev)
        if tier is None:
            return _hard_coded_entry("PVS1", "not_met", None, why)
        return _hard_coded_entry(
            "PVS1", "met", tier,
            f"Initiation-codon variant ({consequence}) in established LoF gene "
            f"— {why}",
        )
    if is_lof and gene_lof and not nmd_escape:
        is_canonical_splice = (
            "splice_acceptor" in consequence or "splice_donor" in consequence
        )
        splice = _pvs1_splice_strength(vep, ev) if is_canonical_splice else None
        if splice is None and is_canonical_splice:
            splice = _pvs1_splice_cryptic_strength(vep, ev)
        if splice and splice[0] != "PVS1":
            return _hard_coded_entry(
                "PVS1", "met", splice[0],
                f"Canonical splice variant ({consequence}) in established LoF "
                f"gene — {splice[1]}; {splice[0]} per ClinGen PVS1 splice arm",
            )
        return _hard_coded_entry(
            "PVS1", "met", "PVS1",
            f"Null variant ({consequence}) in established LoF gene",
        )
    if is_lof and gene_lof and nmd_escape:
        tier, why = _pvs1_nmd_escape_strength(vep, ev)
        return _hard_coded_entry(
            "PVS1", "met", tier,
            f"LoF variant ({consequence}) {why}",
        )
    if is_lof and not gene_lof:
        _gnomad_constraint = (
            ((ev.get("gnomad") or {}).get("gene") or {}).get("gnomad_constraint")
            or {}
        )
        return _hard_coded_entry(
            "PVS1", "not_met", None,
            f"LoF consequence but LoF mechanism not established: "
            f"CHDgene listed={(ev.get('chdgene') or {}).get('listed')}, "
            f"GenCC best={(ev.get('gencc') or {}).get('best_classification') or 'none'}, "
            f"gnomAD pLI={_gnomad_constraint.get('pLI')}, "
            f"LOEUF={_gnomad_constraint.get('oe_lof_upper')}",
        )
    return _hard_coded_entry(
        "PVS1", "not_met", None,
        f"Not a null variant (consequence: {consequence or 'unknown'})",
    )


_BP3_FUNCTIONAL_FEATURE_TYPES = frozenset({
    "Domain", "Motif", "Zinc finger", "DNA binding",
    "Active site", "Binding site", "Transmembrane", "Signal",
})

_BP3_FUNCTIONAL_REPEAT_RE = re.compile(
    r"\b("
    r"wd|ank|ankyrin|lrr|leucine-rich|tpr|tetratricopeptide|"
    r"arm|armadillo|kelch|ef-hand|heat|spectrin|fibronectin|"
    r"immunoglobulin|ig-like|cadherin|laminin|egf|sushi|kringle|"
    r"annexin|pentatricopeptide|nebulin|zinc.finger"
    r")\b",
    re.IGNORECASE,
)


def _bp3_functionless_repeat(ev: dict, vep: dict) -> tuple[bool, str]:
    """Is this variant inside an annotated repeat that has NO known function?

    Returns ``(qualifies, reason)`` — ``reason`` is the evidence string for
    either outcome. This replaces ``vep["in_repeat_region"]``, a key NOTHING in
    production ever populated (only tests set it), which made BP3 dead code and
    left PM4's repeat carve-out permanently disengaged.

    FAIL-CLOSED: BP3 is a benign criterion, so every unknown resolves to "does
    not qualify". Missing UniProt annotation or an unresolvable protein position
    means we cannot establish "without a known function", and asserting it
    anyway would be a false-benign. Mirrors BP7's fail-closed conservation gate.
    """
    pos = _hgvsp_int_position(_mane_hgvsp_from_vep(vep) or "") or (
        vep.get("protein_start") or _hgvsp_int_position(vep.get("hgvsp") or "")
    )
    try:
        pos = int(pos)
    except (TypeError, ValueError):
        pos = None
    if not pos:
        return False, (
            "BP3 requires the variant's protein position to locate an annotated "
            "repeat; it could not be resolved — BP3 withheld (fail-closed)"
        )
    up = ev.get("uniprot") or {}
    if not up.get("ok"):
        return False, (
            "UniProt feature annotation unavailable — BP3 requires a repeat "
            "region with NO known function, which cannot be established "
            "without it (fail-closed)"
        )
    feats = [f for f in (up.get("features") or []) if isinstance(f, dict)]

    def _span(f) -> tuple[int, int] | None:
        try:
            s, e = int(f.get("start")), int(f.get("end"))
        except (TypeError, ValueError):
            return None
        return (s, e) if s <= e else (e, s)

    hits = []
    for f in feats:
        if f.get("type") != "Repeat":
            continue
        sp = _span(f)
        if sp and sp[0] <= pos <= sp[1]:
            hits.append((f, sp))
    if not hits:
        return False, (
            f"Residue {pos} does not fall in any UniProt-annotated repeat "
            "region — BP3 not applicable"
        )

    functional: list[tuple[tuple[int, int], str]] = []
    for f in feats:
        if f.get("type") in _BP3_FUNCTIONAL_FEATURE_TYPES:
            sp = _span(f)
            if sp:
                functional.append((sp, f.get("description") or f.get("type") or ""))
    for d in (up.get("domains") or []):
        sp = _span(d)
        if sp:
            functional.append((sp, d.get("name") or "domain"))
    for b in (up.get("binding_sites") or []):
        sp = _span(b)
        if sp:
            functional.append((sp, b.get("description") or "binding site"))
    for a in (up.get("active_sites") or []):
        try:
            p = int(a.get("position"))
        except (TypeError, ValueError):
            continue
        functional.append(((p, p), a.get("description") or "active site"))

    for f, (rs, re_) in hits:
        desc = f.get("description") or "(unnamed)"
        named = _BP3_FUNCTIONAL_REPEAT_RE.search(desc)
        if named:
            return False, (
                f"Residue {pos} lies in UniProt repeat '{desc}' ({rs}-{re_}), but "
                f"this is a named functional repeat family ({named.group(1)}) — a "
                "structural fold, not a repetitive region without known function; "
                "BP3 does not apply"
            )
        overlap = [
            (sp, why) for sp, why in functional
            if sp[0] <= re_ and rs <= sp[1]
        ]
        if overlap:
            (os_, oe), why = overlap[0]
            return False, (
                f"Residue {pos} lies in UniProt repeat '{desc}' ({rs}-{re_}), but "
                f"that repeat overlaps an annotated functional feature "
                f"('{why}', {os_}-{oe}) — BP3 requires a repetitive region "
                "WITHOUT a known function"
            )
        return True, (
            f"In-frame change at residue {pos} in UniProt repeat '{desc}' "
            f"({rs}-{re_}) with no overlapping domain, motif, active/binding "
            "site, transmembrane or signal annotation — repetitive region "
            "without a known function"
        )
    return False, "BP3 not applicable"


_PM4_REDIRECT_RE = re.compile(
    r"PVS1 is not applicable.{0,200}?"
    r"truncating variants that do NOT undergo nonsense mediated decay"
    r"|see\s+PM4\s+for\s+truncating variants that do NOT undergo NMD",
    re.IGNORECASE | re.DOTALL)


_BP7_NEEDS_BP4_RE = re.compile(
    r"in conjunction with BP4|after assignment of BP4", re.IGNORECASE)


def _bp7_requires_bp4(gene: str | None) -> bool:
    """Whether this gene's VCEP makes BP4 a precondition for BP7.

    Fails open (False) on an unknown gene or a missing record, so BP7 keeps its
    current behaviour everywhere the mandate is not published. Never raises."""
    g = (gene or "").strip().upper()
    if not g:
        return False
    bp7 = _crit_spec(g, "BP7")
    if not bp7:
        return False
    blob = " ".join([
        str(bp7.get("comments") or ""),
        str(bp7.get("vcep_specifications") or ""),
        *[str(v) for v in (bp7.get("strengths") or {}).values()],
    ])
    return bool(_BP7_NEEDS_BP4_RE.search(blob))


def _pm4_nmd_redirect_gene(gene: str | None) -> bool:
    """Whether this gene's PVS1 record redirects NMD-escaping truncating
    variants to PM4, and PM4 is applicable to take them. Never raises."""
    g = (gene or "").strip().upper()
    if not g:
        return False
    if _criterion_applicable(g, "PVS1"):
        return False
    pm4 = _crit_spec(g, "PM4")
    if pm4.get("applicability") != "applicable":
        return False
    pvs1 = _crit_spec(g, "PVS1")
    blob = " ".join([
        str(pm4.get("comments") or ""),
        str(pm4.get("vcep_specifications") or ""),
        *[str(v) for v in (pm4.get("strengths") or {}).values()],
        str(pvs1.get("comments") or ""),
        str(pvs1.get("vcep_specifications") or ""),
    ])
    return bool(_PM4_REDIRECT_RE.search(blob))


def _eval_pm4_bp3(vep: dict, consequence: str, ev: dict | None = None,
                  gene: str | None = None) -> list[dict]:
    """── PM4 / BP3 — Length change in / out of repeat region ──────────
    ACMG/AMP 2015: PM4 = "Protein length changes as a result of in-frame
    deletions/insertions in a non-repeat region OR stop-loss variants." The
    non-repeat qualifier binds only to the in-frame-indel clause — a stop-loss
    changes protein length wherever it sits, and BP3 ("in-frame indels in a
    repetitive region") has no stop-loss arm to catch it. The gate used to AND
    `not repeat_region` across BOTH clauses, so a stop-loss inside a repeat
    region fell through to the terminal branch and was reported as "Not an
    in-frame indel or stop-loss (consequence: stop_lost)" — a message that
    contradicts its own input. Stop-loss is now evaluated independently of
    repeat status. See tests/test_criteria_gate_spec_audit.py."""
    in_frame = "inframe" in consequence
    stop_lost = "stop_lost" in consequence
    truncating = ("stop_gained" in consequence) or ("frameshift" in consequence)
    if (truncating and vep.get("nmd_escape") is True
            and _pm4_nmd_redirect_gene(gene)):
        return [
            _hard_coded_entry(
                "PM4", "met", "PM4_Moderate",
                f"Truncating variant that escapes NMD ({consequence}) — the "
                f"{(gene or '').upper()} VCEP withdraws PVS1 and directs these "
                f"to PM4 ('See PM4 for truncating variants that do NOT undergo "
                f"NMD'), where PM4 is applicable at Moderate",
            ),
            _hard_coded_entry(
                "BP3", "not_met", None,
                f"BP3 applies to in-frame indels only "
                f"(consequence: {consequence or 'unknown'})",
            ),
        ]
    if not (in_frame or stop_lost):
        return [
            _hard_coded_entry(
                "PM4", "not_met", None,
                f"Not an in-frame indel or stop-loss "
                f"(consequence: {consequence or 'unknown'})",
            ),
            _hard_coded_entry(
                "BP3", "not_met", None,
                f"BP3 applies to in-frame indels only "
                f"(consequence: {consequence or 'unknown'})",
            ),
        ]
    bp3_ok, bp3_why = (
        _bp3_functionless_repeat(ev or {}, vep) if in_frame else (False, "")
    )
    if bp3_ok:
        return [
            _hard_coded_entry(
                "PM4", "not_met", None,
                "In-frame indel in a repetitive region without known function "
                "(see BP3) — protein-length change carries no weight here",
            ),
            _hard_coded_entry("BP3", "met", "BP3_Supporting", bp3_why),
        ]
    return [
        _hard_coded_entry(
            "PM4", "met", "PM4_Moderate",
            "Stop-loss variant — protein length change"
            if stop_lost else
            "In-frame indel outside any functionless repetitive region — "
            "protein length change",
        ),
        _hard_coded_entry(
            "BP3", "not_met", None,
            bp3_why or "BP3 applies to in-frame indels only (stop-loss has no "
                       "BP3 arm)",
        ),
    ]


_DENOVO_TABLE1 = {
    "highly_specific":          (2.0, 1.0),
    "consistent_not_specific":  (1.0, 0.5),
    "consistent_heterogeneous": (0.5, 0.25),
    "not_consistent":           (0.0, 0.0),
}
_DENOVO_TABLE2 = ((4.0, "VeryStrong"), (2.0, "Strong"), (1.0, "Moderate"),
                  (0.5, "Supporting"))

_DENOVO_PHENOTYPE_ROW = {
    "RASopathy": "highly_specific",
    "PotassiumChannel": "consistent_not_specific",
}


def _denovo_points_tier(
    confirmed_n, unconfirmed_n, gene, code, pp4_met=False,
) -> str | None:
    """SVI Table 1 + Table 2 strength for `code` on `gene`, or None when this
    gene's VCEP does not publish a points ladder or does not state which
    Table 1 row applies.

    Returns a bare tier name clamped to the rungs the gene actually publishes:
    RASopathy PS2 has no Supporting rung and RASopathy PM6 has no VeryStrong
    rung, so a points total landing on an unpublished tier drops to the highest
    published rung at or below it. Never raises.
    """
    def _i(v):
        try:
            return int(v) if v is not None and not isinstance(v, bool) and int(v) > 0 else 0
        except (TypeError, ValueError):
            return 0

    strengths = _crit_spec(gene, code).get("strengths") or {}
    if len(strengths) < 2:
        return None
    row = _DENOVO_PHENOTYPE_ROW.get((_vcep_freq(gene) or {}).get("vcep"))
    if row is None:
        return None
    if pp4_met and row == "consistent_not_specific":
        row = "highly_specific"

    conf_pt, unconf_pt = _DENOVO_TABLE1[row]
    points = _i(confirmed_n) * conf_pt + _i(unconfirmed_n) * unconf_pt
    tier = next((t for threshold, t in _DENOVO_TABLE2 if points >= threshold), None)
    if tier is None:
        return None

    published = {k.replace(" ", ""): k for k in strengths}
    order = ["Supporting", "Moderate", "Strong", "VeryStrong"]
    for cand in reversed(order[: order.index(tier) + 1]):
        if cand in published:
            return cand
    return None


def _eval_ps2(denovo_confirmed: bool, gene: str | None = None,
              confirmed_n=0, unconfirmed_n=0, trio_stated: bool = True,
              denovo_status: str = "") -> dict:
    """── PS2 — De novo confirmed ──────────────────────────────────────
    Confirmed de novo (paternity + maternity) = PS2_Strong at base 2015
    strength, and it stays there. There is deliberately no escalation to
    PS2_VeryStrong from a count of independent occurrences: the ClinGen SVI
    de-novo SOP (Biesecker/Harrison 2018) awards Very Strong only via a POINTS
    total that weights phenotypic specificity (2 confirmed de novo reach Very
    Strong ONLY for a highly specific phenotype; a merely-consistent phenotype
    is Strong), and several VCEPs cap/condition it on absence of family history.
    That phenotype-points input is not collected, so escalating on count alone
    over-calls (a spurious +4 can flip LP→P) — see the cycle-3 PS2 ceiling note.

    A `denovo_count` parameter used to arrive here and did nothing but append a
    clause to the evidence string explaining why it was being ignored; the form
    field behind it has been removed. If phenotype-specificity scoring is ever
    collected, the escalation belongs here, driven by that input.

    🔴 RE-SCOPED IN FULL 2026-09-07 AND STILL NOT IMPLEMENTABLE. Do not reopen
    on the strength of "the occurrence counts are available" — they are, and the
    count was never the blocker. Measured across all 27 spec'd genes:

      * the 9 genes that STATE which phenotype-specificity row to use (the
        Cardiomyopathy 8 and KCNQ1) publish only ONE PS2 rung, Strong, so there
        is nothing to escalate to;
      * the 16 RASopathy genes that publish THREE rungs as point thresholds
        ("1 Point." / "2 Points." / "4 Points.") state no phenotype row and
        carry no points-per-proband table;
      * every PS2 and PM6 strength shortfall observed is on those 16 genes.

    So no gene in the harvest holds both the ladder and the input to climb it.
    Worse, this function's PS2_Strong for ONE confirmed de novo is only
    consistent with the table's "highly specific" row (2 points), which is the
    row the 9 stating VCEPs explicitly tell you not to default to. Adopting
    their "consistent but not highly specific" row instead makes one confirmed
    de novo 1 point = MODERATE, dropping every current PS2 by 2 points on 15
    rows. The two internally-consistent readings move the same rows ~30 points
    in opposite directions.

    Full scope, with the numbers and what would actually unblock it:
    Recorded in the project's own scoping notes, outside this repository.

    `denovo_confirmed` is trio-only (see app.py): a duo tests one parent and so
    cannot establish both maternity and paternity, which ACMG PS2 requires. A
    duo asserting de novo routes to PM6 instead of satisfying PS2."""
    if not denovo_confirmed:
        _why2 = {
            "inherited_affected":
                "Not de novo — reported as inherited from an AFFECTED parent. "
                "That is co-segregation evidence, so it is counted toward PP1 "
                "as one affected carrier relative rather than as a de novo "
                "criterion.",
            "inherited_unaffected":
                "Not de novo — reported as inherited from an UNAFFECTED "
                "parent. That is a non-penetrance observation, and it is NOT "
                "BS4 (which requires an AFFECTED relative who does NOT carry "
                "the variant). Its only home is BS2, which this gene's VCEP "
                "either marks Not Applicable for reduced penetrance and adult "
                "onset, or scores only from several phenotyped unaffected "
                "carriers — so a single unaffected parent earns nothing.",
            "unconfirmed":
                "Not PS2 — the variant is absent from both parents but "
                "parental relationships are not genetically confirmed; see "
                "PM6.",
            "": "De novo status was not assessed.",
        }.get((denovo_status or "").strip().lower(),
              "De novo not confirmed by parental testing")
        return _hard_coded_entry("PS2", "not_met", None, _why2)
    _n_conf = max(int(confirmed_n or 0), 1) if denovo_confirmed else 0
    _tier = _denovo_points_tier(_n_conf, unconfirmed_n, gene, "PS2")
    if _tier and _tier != "Strong":
        _pts = _DENOVO_TABLE1[
            _DENOVO_PHENOTYPE_ROW[(_vcep_freq(gene) or {}).get("vcep")]]
        return _hard_coded_entry(
            "PS2", "met", f"PS2_{_tier}",
            f"De novo confirmed by parental testing (full trio), "
            f"{_n_conf} occurrence(s) with confirmed parental relationships"
            + (f" and {int(unconfirmed_n)} with assumed relationships"
               if int(unconfirmed_n or 0) else "")
            + f". {gene} VCEP scores this at {_tier} under the ClinGen SVI de "
              f"novo points system (v1.1 Table 1: {_pts[0]} point(s) per "
              f"confirmed occurrence for this panel's phenotype category, "
              f"{_pts[1]} per assumed; Table 2 thresholds 1/2/4 = "
              f"Moderate/Strong/VeryStrong)",
        )
    return _hard_coded_entry(
        "PS2", "met", "PS2_Strong",
        "De novo status confirmed by parental testing (full trio)"
        if trio_stated else
        "De novo reported as CONFIRMED by the curator; parental testing detail "
        "was not specified. ACMG PS2 requires both maternity and paternity to "
        "be confirmed — recorded on the curator's assertion",
    )


def _eval_pm6(trio: str, denovo_status: str, inheritance: str,
              denovo_confirmed: bool, ps2_met: bool,
              gene: str | None = None, unconfirmed_n=0) -> dict:
    """── PM6 — De novo assumed (mutually exclusive with PS2) ──────────
    Assumed de novo (paternity/maternity not confirmed) = PM6_Moderate at base
    2015 strength. As with PS2 above, there is no escalation to PM6_Strong from
    a count of occurrences — the SVI points escalation needs a
    phenotype-specificity input HeartVar does not collect, and a count-only
    upgrade over-calls on the VCEP-covered genes.

    A DUO that asserts "confirmed de novo" lands here, not on PS2: one parent
    tested cannot establish both maternity and paternity, which is exactly the
    "without confirmation of paternity and maternity" case PM6 exists for. That
    state used to satisfy PS2's gate (+4 instead of +2, with a false "full trio"
    evidence string) — see the denovo_confirmed comment in app.py and
    tests/test_criteria_gate_spec_audit.py. Only trio + confirmed reaches PS2,
    and PS2 still precludes PM6 below."""
    pm6_eligible = (
        denovo_status in ("unconfirmed", "confirmed")
        or inheritance == "DN"
    )
    if trio == "duo" and denovo_status in ("unconfirmed", "confirmed"):
        return _hard_coded_entry(
            "PM6", "not_met", None,
            "Only ONE parent was tested (duo), so no de novo criterion is "
            "applied. Both PS2 and PM6 require the variant to be absent from "
            "BOTH parents — with one parent untested it may simply have been "
            "inherited from the untested parent. ACMG PS2: \"Confirmation of "
            "paternity only is insufficient.\" Testing the second parent "
            "would make this PM6 (parentage assumed) or PS2 (parentage "
            "confirmed).",
        )

    if pm6_eligible and not denovo_confirmed:
        _why6 = (
            "De novo observed — the variant was absent from both parents — but "
            "parental relationships are not genetically confirmed, so PM6 "
            "applies rather than PS2 (ClinGen SVI de novo v1.1 prices an "
            "unconfirmed-parentage occurrence at half a confirmed one)"
        )
        _n6 = max(int(unconfirmed_n or 0), 1)
        _tier6 = _denovo_points_tier(0, _n6, gene, "PM6")
        if _tier6 and _tier6 != "Moderate":
            return _hard_coded_entry(
                "PM6", "met", f"PM6_{_tier6}",
                f"{_why6}. {_n6} occurrence(s) with assumed parental "
                f"relationships; {gene} VCEP scores this at {_tier6} under the "
                f"ClinGen SVI de novo points system (v1.1 Tables 1 and 2)",
            )
        return _hard_coded_entry("PM6", "met", "PM6_Moderate", _why6)
    if ps2_met:
        return _hard_coded_entry(
            "PM6", "not_met", None,
            "Precluded by PS2 (confirmed de novo)",
        )
    return _hard_coded_entry(
        "PM6", "not_met", None,
        "De novo not asserted by inheritance/trio inputs",
    )


def _eval_pm3(gene: str | None, inheritance: str, zygosity: str,
              in_trans_pathogenic: str, denovo_status: str) -> dict:
    """── PM3 — In trans with pathogenic (recessive) ───────────────────
    PM3 is a recessive-disorder criterion: the variant is detected in trans
    with a pathogenic/likely-pathogenic allele in an affected proband. Two
    deterministic paths:
      (1) confirmed phase — curator stated in_trans_pathogenic == "yes"
          (a second P/LP allele confirmed IN TRANS). This is the stronger,
          phase-confirmed observation.
      (2) the existing inherited-from-affected-parent path, where phase is
          only inferred (kept for backward-compat; not phase-confirmed).
    VCEP applicability gate: PM3 is only marked Applicable for the recessive-
    mechanism genes in our CSpec table (LZTR1 / KCNQ1); the cardiomyopathy /
    RASopathy / FBN1 VCEPs mark PM3 Not Applicable. _criterion_applicable is
    default-open, so uncovered genes still evaluate PM3 (no silent loss)."""
    pm3_applicable = _criterion_applicable(gene, "PM3")
    pm3_in_trans = (
        inheritance == "AR"
        and zygosity == "het"
        and in_trans_pathogenic == "yes"
    )
    pm3_inherited = (
        inheritance == "AR"
        and zygosity == "het"
        and denovo_status == "inherited_affected"
    )
    if not pm3_applicable and (pm3_in_trans or pm3_inherited):
        return _hard_coded_entry(
            "PM3", "not_met", None,
            f"{gene} VCEP marks PM3 Not Applicable — in-trans / recessive "
            "evidence is not scored under PM3 for this gene (CSpec)",
        )
    if pm3_in_trans:
        return _hard_coded_entry(
            "PM3", "met", "PM3_Moderate",
            "Heterozygous variant in AR gene confirmed IN TRANS with a "
            "pathogenic/likely-pathogenic allele in an affected proband",
        )
    if pm3_inherited:
        return _hard_coded_entry(
            "PM3", "met", "PM3_Supporting",
            "Heterozygous variant in AR gene, inherited from affected "
            "parent — possible trans configuration, phase NOT confirmed. "
            "ClinGen SVI PM3 v1.0 Table 1 awards 0.5 points for a "
            "phase-unknown observation with a P/LP allele (vs 1.0 confirmed), "
            "and Table 2 puts 0.5 points at PM3_Supporting.",
        )
    return _hard_coded_entry(
        "PM3", "not_met", None,
        "PM3 conditions not satisfied (requires AR + heterozygous + "
        "in trans with a pathogenic allele, or inherited from an "
        "affected parent)",
    )


def _eval_bp2(inheritance: str, zygosity: str, in_trans_pathogenic: str) -> dict:
    """── BP2 — In trans / cis with pathogenic (dominant) ──────────────
    ACMG/AMP 2015 (Richards, Table 4): "Observed in trans with a pathogenic
    variant for a fully penetrant dominant gene/disorder; or observed in cis
    with a pathogenic variant in any inheritance pattern." BOTH arms require a
    SECOND, pathogenic allele in the same gene. The criterion turns on the
    PHASE of that second allele — not on which parent transmitted this one.

    This previously fired on (AD + het + denovo_status == "inherited_unaffected")
    and never consulted in_trans_pathogenic at all. That conflated BP2 with a
    penetrance/segregation argument ("an unaffected parent carries it too"),
    which under ACMG is BS2/BS4 territory and is NOT what BP2 asserts. The
    practical harm was a false-benign: BP2 is often the only benign criterion
    to fire, and a lone benign-Supporting (-1) IS Likely benign under the
    Tavtigian point system, so any AD het variant transmitted by an unaffected
    parent was auto-demoted to LB on evidence the curator never entered. See
    backend/tests/test_bp2_in_trans_spec.py.

    Only the in-trans arm is evaluable here: the app collects phase-confirmed
    in-trans status (in_trans_pathogenic, shared with PM3) but has no in-cis
    input, so the cis arm stays not_met rather than being guessed at.

    Inheritance splits the two in-trans criteria cleanly and mutually
    exclusively: AR + het + in trans → PM3 (pathogenic direction, the second
    hit completes a biallelic genotype); AD + het + in trans → BP2 (benign
    direction, the pathogenic allele alone already explains the phenotype).

    The "fully penetrant" qualifier is left to the curator: penetrance is not
    derivable per-gene from the data held at runtime, and reduced penetrance
    weakens the inference, so the caveat is surfaced in the evidence string
    instead of silently assumed."""
    if (
        inheritance == "AD"
        and zygosity == "het"
        and in_trans_pathogenic == "yes"
    ):
        return _hard_coded_entry(
            "BP2", "met", "BP2_Supporting",
            "Heterozygous variant in a dominant gene confirmed IN TRANS with "
            "a pathogenic/likely-pathogenic allele, which alone accounts for "
            "the phenotype (assumes full penetrance — reassess if the "
            "gene/disorder shows reduced penetrance)",
        )
    return _hard_coded_entry(
        "BP2", "not_met", None,
        "BP2 conditions not satisfied (requires a second pathogenic allele: "
        "AD + heterozygous + confirmed in trans with a pathogenic variant, or "
        "in cis with a pathogenic variant — no in-cis input is collected). "
        "Inheritance from an unaffected parent is not BP2 evidence",
    )


def _eval_bs2(gene: str | None, inheritance: str,
              gnomad_hom: int, gnomad_hemi: int,
              ba1_exception: str | None = None) -> dict:
    """── BS2 — Observed in healthy adults ─────────────────────────────
    ACMG/AMP 2015: "Observed in a healthy adult individual for a recessive
    (HOMOZYGOUS), dominant (HETEROZYGOUS), or X-linked (HEMIZYGOUS) disorder,
    with full penetrance expected at an early age." Base strength Strong (-4).

    VCEP-governed for covered genes (verified per CSpec 2026-06-05): the
    Cardiomyopathy (GN002 etc.), Long-QT/KCNQ1 (GN112) and FBN1 (GN022) VCEPs
    mark BS2 **Not Applicable** (reduced/incomplete penetrance, variable
    expressivity, adult onset). The RASopathy VCEP keeps BS2 but restricts it to
    phenotyped unaffected FAMILY members and EXPLICITLY forbids general
    population/gnomAD data. So a gnomAD-based BS2 must never fire for a covered
    gene — suppressed below. Out-of-table genes use the calibrated ladder.

    COUNT LADDER (was: full Strong at `gnomad_hom > 5`, nothing below). Modern
    VCEP practice tiers BS2 by observation count instead of firing Strong off a
    single threshold. We adopt the ClinGen Pulmonary Arterial Hypertension /
    BMPR2 specification — an autosomal-dominant, incomplete-penetrance disease,
    the closest published analogue to the cardiac genes here:
        >= 3 healthy homozygotes/hemizygotes → BS2_Strong (-4)
        <= 2                                 → not_met
    (Eichstaedt et al., Hum Mutat 2025, Table 1: "≥ 3 counts: BS2_strong;
    ≥ 2 counts: BS2_supp"; CSpec BMPR2.) The old >5 rule missed the 3-5 band
    entirely and jumped straight to Strong above it.

    DELIBERATE DEVIATION: BMPR2 publishes BS2_Supporting at 2 counts; we withhold
    below 3. BS2 is a BENIGN criterion and a false-benign closes the diagnostic
    question, so the project errs toward VUS at the weakest count. This deviates
    from the published ladder in the CONSERVATIVE direction only — the >=3 Strong
    bar and the proband-zygosity fix both stand.

    REDUCED-PENETRANCE EXCLUSION: ACMG BS2 requires "full penetrance expected at
    an early age". The nine variants on ClinGen's BA1 exception list (Ghosh 2018,
    see _BA1_EXCEPTION_VARIANTS) are all common REDUCED-PENETRANCE / hypomorphic
    alleles — that is precisely why they are common yet not benign (HFE C282Y
    haemochromatosis, GJB2 V37I mild hearing loss, BTD D444H partial biotinidase
    deficiency, ACADS R171W largely biochemical SCAD). BS2's own precondition
    therefore fails for them no matter how many healthy homozygotes gnomAD holds,
    so no population-frequency benign evidence (BA1, BS1, BS2) is applied. Found
    in practice: HFE c.845G>A has 3,248 gnomAD homozygotes and was reaching
    Likely benign against a ClinGen PATHOGENIC classification.

    HETEROZYGOUS COUNTS ARE DELIBERATELY UNUSED. BMPR2 is explicit that BS2
    "cannot be used for heterozygotes due to incomplete penetrance" — the same
    reasoning that makes the cardiac VCEPs disable BS2 altogether. Only
    homozygous (autosomal) and hemizygous (X-linked) healthy observations count.

    NOT GATED ON THE PROBAND'S ZYGOSITY. The previous gate required
    `zygosity == "hom"` before consulting the gnomAD homozygote count, but BS2
    is a statement about the healthy CONTROL's genotype, not the proband's: a
    heterozygous proband whose variant is homozygous in several healthy adults
    is precisely the BS2 observation, and it was being dropped.

    CAVEAT (not machine-checkable): gnomAD is not phenotyped and is not screened
    for adult-onset cardiac disease, so "healthy adult ... full penetrance
    expected at an early age" is an assumption for out-of-table genes. This is
    exactly why the covered VCEPs forbid population-based BS2; the suppression
    below is the guard. See tests/test_bs2_calibration.py."""
    if ba1_exception is not None:
        return _hard_coded_entry(
            "BS2", "not_met", None,
            f"On the ClinGen SVI BA1 exception list ({ba1_exception}) — a common "
            "but REDUCED-PENETRANCE allele that ClinGen holds Pathogenic/VUS "
            "despite a high population frequency. ACMG BS2 requires the healthy "
            "adult observation to come from a disorder with FULL PENETRANCE "
            "EXPECTED AT AN EARLY AGE, which is not the case here, so healthy "
            "homozygotes are expected and are not benign evidence — BS2 withheld "
            "(no population-frequency benign evidence is applied to this variant)",
        )
    _bs2_spec = _vcep_freq(gene)
    _bs2_vcep = (_bs2_spec or {}).get("vcep")
    if _bs2_vcep == "RASopathy":
        return _hard_coded_entry(
            "BS2", "not_met", None,
            "RASopathy VCEP: BS2 requires phenotyped unaffected FAMILY "
            "members and forbids general-population/gnomAD data — not "
            "assessable from population frequency data",
        )
    # for correctness of attribution, not for the answer: the frequency table
    if not _criterion_applicable(gene, "BS2"):
        return _hard_coded_entry(
            "BS2", "not_met", None,
            f"{gene} VCEP marks BS2 Not Applicable (reduced/incomplete "
            "penetrance, variable expressivity, adult onset)",
        )
    if _bs2_spec:
        return _hard_coded_entry(
            "BS2", "not_met", None,
            f"{_bs2_vcep} VCEP marks BS2 Not Applicable "
            "(reduced/incomplete penetrance, variable expressivity, adult onset)",
        )

    if inheritance == "MT":
        return _hard_coded_entry(
            "BS2", "not_met", None,
            "Mitochondrial inheritance: expression depends on HETEROPLASMY "
            "LEVEL, not zygosity, and gnomAD's mtDNA homoplasmy/heteroplasmy "
            "counts cannot establish the healthy-adult homozygous observation "
            "BS2 requires — withheld (fail-closed)",
        )

    n_hom = max(0, int(gnomad_hom or 0))
    n_hemi = max(0, int(gnomad_hemi or 0))
    n = max(n_hom, n_hemi)
    kind = "homozygous" if n_hom >= n_hemi else "hemizygous"
    if n >= 3:
        return _hard_coded_entry(
            "BS2", "met", "BS2_Strong",
            f"Observed {kind} in {n} gnomAD control(s) — >=3 healthy "
            "homozygotes/hemizygotes is BS2 at Strong (ClinGen BMPR2 VCEP "
            "count ladder; base ACMG BS2 strength)",
        )
    return _hard_coded_entry(
        "BS2", "not_met", None,
        f"Only {n} healthy homozygous/hemizygous gnomAD observation(s) — BS2 "
        "requires at least 3. Heterozygous controls are not counted (incomplete "
        "penetrance), and the 2-count Supporting tier the BMPR2 VCEP publishes "
        "is deliberately not applied here (benign evidence is held to the "
        "stronger bar)",
    )


def compute_hard_coded_criteria(
    ev: dict, clinical_context: dict, gene: str | None = None,
) -> list[dict]:
    """Deterministically evaluate the 18 ACMG/AMP criteria this function owns.

    NOT the whole ``hard_coded_criteria_codes`` set, which is 21. PM1, PP1 and
    BS4 joined that set on 2026-09-08 but are derived elsewhere, by
    ``no_ai.infer_supplementary_criteria``, because they read the curator's
    structured segregation inputs rather than the evidence blob. Callers that
    need all 21 must add them; ``app._python_authoritative_supplementary`` is
    the one place that does, and it fills any that did not fire with a
    ``not_assessed`` placeholder. Reading this docstring as "one entry per code
    in the JSON set" is what shipped an 18-entry PRECOMPUTED block under prompt
    text promising 21.

    Returns normalised criterion dicts, each carrying ``source="hard_coded"``. Missing or ambiguous inputs yield
    ``status="not_met"`` with an evidence string explaining why, never raise.

    The caller combines these results with the AI's returned criteria for
    the remaining interpretive codes, applies cross-criterion mutual
    exclusion, and computes the final points + classification."""
    ev = ev or {}
    clinical_context = clinical_context or {}
    inheritance = (clinical_context.get("inheritance_input") or "").upper()
    zygosity = (clinical_context.get("zygosity") or "").lower()
    trio = (clinical_context.get("trio_status") or "").lower()
    denovo_status = (clinical_context.get("denovo_status") or "").lower()
    denovo_confirmed = bool(clinical_context.get("denovo_confirmed"))
    denovo_conf_n = clinical_context.get("denovo_confirmed_count") or 0
    denovo_unconf_n = clinical_context.get("denovo_unconfirmed_count") or 0
    in_trans_pathogenic = (clinical_context.get("in_trans_pathogenic") or "").lower()

    vep = ev.get("vep") or {}
    consequence = (vep.get("most_severe_consequence") or "").lower()

    gnomad_af, gnomad_af_basis, gnomad_an = _vcep_freq_af(ev, gene)
    _gn = ev.get("gnomad") or {}
    gnomad_confirmed_absent = (
        bool(_gn.get("ok"))
        and _gn.get("variant_found") is False
        and not _gn.get("indel_unresolved")
    )
    gnomad_hom, gnomad_hemi = _gnomad_hom_hemi(ev)
    spliceai_max = _spliceai_per_score_max(ev)

    out: list[dict] = []


    ba1_threshold = _ba1_threshold(gene)
    _ba1_basis = "VCEP-specific" if _vcep_freq(gene) else "ACMG/SVI 5% default"
    ba1 = _eval_ba1(
        gnomad_af, ba1_threshold, _ba1_basis,
        gnomad_an=gnomad_an,
        ba1_exception=_ba1_exception_hit(gene, vep),
        af_basis=gnomad_af_basis,
    )
    out.append(ba1)
    ba1_met = ba1["status"] == "met"

    bs1_threshold = _bs1_threshold(ev, gene, inheritance)
    freq_src = _freq_source_label(ev, gene, inheritance)
    if not _vcep_freq(gene):
        log.info(
            "freq fallback: gene=%s inh=%s basis=%s (BS1 floor=%s, PM2 ceiling=%s)",
            gene, inheritance or "AD/unknown", freq_src,
            bs1_threshold, _pm2_threshold(gene, inheritance, ev),
        )
    bs1 = _eval_bs1(
        gnomad_af, bs1_threshold, ba1_threshold, inheritance, freq_src,
        ba1_met=ba1_met, af_basis=gnomad_af_basis,
    )
    out.append(bs1)
    bs1_met = bs1["status"] == "met"

    pm2_threshold = _pm2_threshold(gene, inheritance, ev)
    _pm2_src = freq_src
    out.append(_eval_pm2(
        gnomad_af, pm2_threshold, _pm2_src, ba1_met, bs1_met,
        gnomad_confirmed_absent, _gn,
    ))

    pp3_ceiling = _insilico_strength_ceiling(gene, ev)
    pp3_bp4, noncoding_splice, splice_clear = _eval_pp3_bp4(
        ev, vep, consequence, spliceai_max, pp3_ceiling,
    )
    out.extend(pp3_bp4)
    sa_ok = bool((ev.get("spliceai") or {}).get("ok"))

    _bp7 = _eval_bp7(
        consequence, spliceai_max, sa_ok, noncoding_splice, splice_clear,
        vep.get("phylop100way"),
    )
    _bp7_cons_ok = (
        vep.get("phylop100way") is not None
        and _as_float(vep.get("phylop100way")) is not None
        and float(_as_float(vep.get("phylop100way"))) <= 2.0
    )
    _bp7["_bp7_rna_routable"] = (
        any(t in consequence for t in ("synonymous", "intron"))
        and "splice_region" not in consequence
        and _bp7_cons_ok
    )
    out.append(_bp7)

    out.append(_eval_pvs1(ev, vep, consequence, spliceai_max))

    out.extend(_eval_pm4_bp3(vep, consequence, ev, gene))

    ps2 = _eval_ps2(denovo_confirmed, gene, denovo_conf_n, denovo_unconf_n,
                    trio_stated=bool(clinical_context.get("denovo_trio_stated", True)),
                    denovo_status=denovo_status)
    out.append(ps2)
    ps2_met = ps2["status"] == "met"

    out.append(_eval_pm6(
        trio, denovo_status, inheritance, denovo_confirmed, ps2_met,
        gene, denovo_unconf_n,
    ))

    out.append(_eval_pm3(
        gene, inheritance, zygosity, in_trans_pathogenic, denovo_status,
    ))

    out.append(_eval_bp2(inheritance, zygosity, in_trans_pathogenic))

    out.append(_eval_bs2(
        gene, inheritance, gnomad_hom, gnomad_hemi,
        ba1_exception=_ba1_exception_hit(gene, vep),
    ))

    out.extend(
        _clinvar_pp5_bp6_criteria(
            ev.get("clinvar"), ba1_met=ba1_met, bs1_met=bs1_met,
        )
    )

    out.append(_clinvar_ps1_criterion(ev.get("clinvar_pm5_candidates"), gene))

    out.append(_clinvar_pm5_criterion(ev.get("clinvar_pm5_candidates"), gene))

    for _entry in out:
        cap_to_spec_strength(_entry, gene)

    return out


def apply_pm1_collision_rules(
    criteria: list[dict], gene: str | None = None,
) -> list[dict]:
    """The PS1/PM5 and PM1/PM5 mutual exclusions, as their own callable pass.

    WHY THIS IS A SEPARATE FUNCTION. `_gate_criteria_applicability` does not
    only suppress criteria — it ASSERTS PM1 to `met` for a variant inside a
    gene's published hotspot range or exon (see the "scopes PM1 to enumerated
    hotspots" branch). That gate runs AFTER `apply_cross_criterion_exclusions`
    at every call site, so a PM1 the gate asserted never had these collisions
    resolved: the deterministic / no-AI path could score PM1 AND PM5 together
    (+4) where the AI path, whose PM1 arrives before the exclusions, scores PM5
    alone (+2).

    Found 2026-09-01 on TNNI3 NM_000363.5:c.575G>A. ai_mode=none returned
    PM1_Moderate + PM5_Moderate = 6 points = Likely pathogenic; ai_mode=server
    returned PM5 alone = 5 points = VUS. Same variant, same phenotype, and the
    no-AI answer violates the Cardiomyopathy CSpec's explicit "PM5 should not be
    combined with PM1". Direction of the defect is an OVER-call, and it was one
    of the cases where adding AI made the tier worse.

    Defined once and called twice — inside apply_cross_criterion_exclusions at
    its original position, and again after the applicability gate — because a
    second copy of these rules would be free to drift from the first.

    Only ever demotes, so calling it again on an already-resolved list is a
    no-op.
    """
    by_code = {c["code"]: c for c in criteria}

    def _demote(code: str, reason: str) -> None:
        c = by_code.get(code)
        if not c or c.get("status") != "met":
            return
        c["status"] = "not_met"
        c["criteria_strength"] = None
        c["evidence"] = reason

    _bp7 = by_code.get("BP7")
    _bp4 = by_code.get("BP4")
    if (
        _bp7 and _bp7.get("status") == "met"
        and _bp7_requires_bp4(gene)
        and not (_bp4 and _bp4.get("status") == "met")
    ):
        _demote("BP7", (
            f"BP7 not applied: the {(gene or '').upper()} VCEP requires it to be "
            f"used in conjunction with BP4, and BP4 is not met. A curator who "
            f"judges the computational evidence benign can apply both by hand."
        ))

    ps1 = by_code.get("PS1")
    pm5 = by_code.get("PM5")
    if (
        ps1 and pm5
        and ps1.get("status") == "met"
        and pm5.get("status") == "met"
    ):
        _demote("PM5", "Precluded by PS1 (same-AA pathogenic at residue)")

    pm1 = by_code.get("PM1")
    if (
        pm1 and pm5
        and pm1.get("status") == "met"
        and pm5.get("status") == "met"
    ):
        pm5_text = _crit_text(gene, "PM5")
        _pm5_spec = _vcep_freq(gene)
        _pm5_is_rasopathy = bool(_pm5_spec) and _pm5_spec.get("vcep") == "RASopathy"
        _cm_clause = "should not be combined with pm1" in pm5_text.lower()
        if _cm_clause:
            _pm1_pts = _strength_points(pm1.get("criteria_strength")) or 2
            _pm5_pts = _strength_points(pm5.get("criteria_strength")) or 2
            if _pm5_pts >= _pm1_pts:
                _demote(
                    "PM1",
                    "Precluded by PM5 (Cardiomyopathy CSpec: \"PM5 should not "
                    "be combined with PM1. If both are applicable at MODERATE "
                    "weight, use of PM5 is most appropriate since it is variant "
                    "specific\") — PM5 retained as the variant-specific "
                    "evidence. NOTE: this removes PM1's points, so the tier CAN "
                    "fall; the CSpec forbids counting both.",
                )
            else:
                _demote(
                    "PM5",
                    "Precluded by PM1 (Cardiomyopathy CSpec does not combine "
                    "PM5 with PM1). The CSpec prefers PM5 when both are "
                    f"MODERATE, but PM1 here is stronger (+{_pm1_pts} vs "
                    f"+{_pm5_pts}), so the stronger criterion is kept — "
                    "HeartVar's tie-break, not a CSpec instruction",
                )
        elif _pm5_is_rasopathy:
            _demote(
                "PM5",
                "Precluded by PM1 (RASopathy VCEP does not co-apply them — its "
                "PM5 needs ≥2 pathogenic residue changes at the codon in ≥5 "
                "probands, and where PM1 fires the residue region is already "
                "covered); PM1 retained. NOTE: this removes PM5's points, so the "
                "tier CAN fall; the VCEP does not co-apply them.",
            )
    return criteria


def apply_cross_criterion_exclusions(
    criteria: list[dict],
    gene: str | None = None,
    inheritance: str | None = None,
    has_curator_segregation: bool = False,
    verify_source_evidence: bool = True,
    zygosity: str | None = None,
    in_trans_pathogenic: str | None = None,
) -> tuple[list[dict], str | None]:
    """Apply the cross-criterion mutual-exclusion rules after the
    hard-coded set and AI set have been combined. Returns
    ``(cleaned_criteria, forced_classification)`` — if BA1 fires the
    forced classification is "Benign" and every pathogenic met criterion
    is demoted to not_met regardless of source.

    ``gene``/``inheritance`` (CSpec applicability + the request's short
    inheritance code) drive the recessive PS4 guard below; both default to
    None so the rule is a no-op when context is unavailable.

    ``has_curator_segregation`` — True when the request carried structured
    curator segregation fields (``clinical_context.seg_affected_carriers``).
    It distinguishes a curator-sourced PP1 (left entirely untouched — the
    existing informed/no-key path is authoritative) from a literature-sourced
    PP1 (uninformed mode), which is hardened by the facts channel: PMID +
    genotyped-carrier gate, Kelly ladder, unaffected temper, cap-at-Moderate,
    and cap-only vs the LLM's emitted strength.

    ``verify_source_evidence`` — False when the caller's criteria carry no
    ``evidence`` prose and no ``facts``, which is the case for POST
    /api/acmg/rescore (the browser posts code/status/strength/direction/
    curator_override and nothing else). The five gates that decide from those
    two fields — the recessive-PS4 case-control escape hatch, the PS4
    proband-count cap, the PP1 literature hardening, and the PS3/BS3 assay+PMID
    gates — are then SKIPPED rather than run blind, because "no evidence in the
    payload" is not the same fact as "no evidence exists" and reading it that
    way demotes criteria that were already policed on the way out (a curator
    editing PP1 was silently dropping an engine PS4 from Moderate to Supporting
    and switching an engine PS3 off entirely). Every ACMG LOGIC rule still runs:
    those read status/strength, which the client does send. Default True keeps
    the curate path — the only path where the LLM's evidence is actually present
    and worth policing — byte-identical.

    Rules applied:
      - BA1 met → Benign + zero all pathogenic met criteria
      - PS2 + PM6 both met → keep PS2 only
      - PP3 + BP4 both met → keep neither (conflicting evidence)
      - PVS1 met → demote PP3 (predictive evidence already in PVS1)
      - PS1 + PM5 both met → keep PS1 only (Strong wins)
      - PM1 + PM5 both met where the VCEP does not combine them. The
        Cardiomyopathy CSpec names which one survives — "If both are
        applicable at MODERATE weight, use of PM5 is most appropriate since it
        is variant specific" → keep PM5, demote PM1. The RASopathy VCEP gives
        no such instruction, so PM1 is kept there. Tier-neutral either way
        (both +2)
      - PM4 + BP3 both met → keep PM4 only (non-repeat wins)
      - PS4 (proband-count, no case-control OR) on a recessive gene whose
        VCEP places in-trans evidence under PM3 → demote (those occurrences
        are PM3, not PS4 case-enrichment; ACMG forbids double-counting)
      - BS3 met without a cited functional assay + PMID → demote (benign
        ClinVar / frequency is not BS3 evidence; safe-direction)
      - PP1 + BS4 both met → keep NEITHER (contradictory segregation readings of
        one pedigree; conflict surfaced for the curator rather than arbitrated)
      - PP1 met from published literature (uninformed mode, no curator
        segregation) → verify facts (PMID + genotyped affected carriers),
        re-derive strength on the Kelly ladder, temper for unaffected carriers,
        cap at Moderate, and CAP-ONLY vs the LLM's strength (never raise)
    """
    by_code = {c["code"]: c for c in criteria}

    def _demote(code: str, reason: str) -> None:
        c = by_code.get(code)
        if not c or c.get("status") != "met":
            return
        c["status"] = "not_met"
        c["criteria_strength"] = None
        c["evidence"] = reason

    def _curator_set(code: str) -> bool:
        """True when the CURATOR switched this criterion on, not the engine.

        Used only by the three VERIFICATION gates below — the PS3 and BS3
        functional-assay prose/PMID gates, and the PP1 literature-facts
        hardening. Those exist to police LLM OUTPUT: the model fired PS3 off an
        in-silico score, or PP1 off a family history with no genotyped carriers.
        They ask "did whatever emitted this actually cite evidence for it?".

        A human ticking the box in the Criteria tab IS the citation — they read
        the paper, and the override is recorded in the report — so running an
        anti-fabrication gate over their choice would silently discard the
        judgement they just entered and hand back a tier that ignores their edit.

        The ACMG LOGIC rules are deliberately NOT skipped: PP1-vs-BS4,
        PP3-vs-BP4, BA1 forcing, PS1-over-PM5 and the rest are properties of the
        framework rather than guesses about a source's honesty, so an override
        obeys them exactly like an engine-applied criterion.

        This exempts the EDITED criterion only. The criteria around it need
        ``verify_source_evidence=False`` — see the parameter docs above.
        See tests/test_acmg_rescore_endpoint.py.
        """
        c = by_code.get(code)
        return bool(c and c.get("curator_override"))

    forced_classification: str | None = None
    ba1 = by_code.get("BA1")
    if ba1 and ba1.get("status") == "met":
        forced_classification = "Benign"
        for c in criteria:
            if c.get("direction") == "pathogenic" and c.get("status") == "met":
                _demote(
                    c["code"],
                    "Precluded by BA1 (gnomAD popmax above stand-alone benign "
                    "cutoff — stand-alone benign)",
                )

    ps2 = by_code.get("PS2")
    pm6 = by_code.get("PM6")
    if (
        ps2 and pm6
        and ps2.get("status") == "met"
        and pm6.get("status") == "met"
    ):
        _demote("PM6", "Precluded by PS2 (confirmed de novo)")

    pp3 = by_code.get("PP3")
    bp4 = by_code.get("BP4")
    if (
        pp3 and bp4
        and pp3.get("status") == "met"
        and bp4.get("status") == "met"
    ):
        _demote("PP3", "Conflicting computational evidence with BP4")
        _demote("BP4", "Conflicting computational evidence with PP3")

    pvs1 = by_code.get("PVS1")
    if (
        pvs1 and pp3
        and pvs1.get("status") == "met"
        and pp3.get("status") == "met"
    ):
        _demote(
            "PP3",
            "Precluded by PVS1 (splice/predictive evidence already counted in "
            "the PVS1 null-variant call — SVI PVS1 tree, no double-counting)",
        )

    pm1 = by_code.get("PM1")
    if (
        pp3 and pm1
        and pp3.get("status") == "met"
        and pm1.get("status") == "met"
    ):
        pp3_pts = _strength_points(pp3.get("criteria_strength"))
        pm1_pts = _strength_points(pm1.get("criteria_strength"))
        allowed = 4 - pp3_pts
        if pp3_pts + pm1_pts > 4:
            reason = (
                "PP3+PM1 correlated-evidence cap (combined ≤ Strong, Pejaver "
                "2022) — PP3 already counts the in-silico signal"
            )
            if allowed >= 2 and pm1_pts > 2:
                pm1["criteria_strength"] = "PM1_Moderate"
                pm1["evidence"] = f"{pm1.get('evidence', '')} [reduced: {reason}]"
            elif allowed >= 1 and pm1_pts > 1:
                pm1["criteria_strength"] = "PM1_Supporting"
                pm1["evidence"] = f"{pm1.get('evidence', '')} [reduced: {reason}]"
            else:
                _demote("PM1", f"Precluded — {reason}")

    apply_pm1_collision_rules(criteria, gene)

    pm4 = by_code.get("PM4")
    bp3 = by_code.get("BP3")
    if (
        pm4 and bp3
        and pm4.get("status") == "met"
        and bp3.get("status") == "met"
    ):
        _demote("BP3", "Precluded by PM4 (variant not in repeat region)")

    bs4 = by_code.get("BS4")
    _zyg = (zygosity or "").strip().lower()
    _trans = (in_trans_pathogenic or "").strip().lower()
    if (
        bs4
        and bs4.get("status") == "met"
        and (inheritance or "").upper() == "AR"
        and _zyg == "het"
        and _trans == "yes"
    ):
        _demote("BS4", (
            "Precluded: the proband is a compound heterozygote in a recessive "
            "disorder, where non-segregation in relatives cannot distinguish "
            "WHICH of the two alleles is benign — only one of them need be. Per "
            "Biesecker 2024 (ClinGen PP1/BS4 guidance) such an observation "
            "provides little to no evidence against pathogenicity, so the -4.0 "
            "BS4 argument is withdrawn. (BS4 is retained for AR with "
            "homozygosity, where it remains informative.)"
        ))

    pvs1 = by_code.get("PVS1")
    if (
        pm4 and pvs1
        and pm4.get("status") == "met"
        and pvs1.get("status") == "met"
    ):
        _demote("PM4", (
            "Precluded by PVS1 ("
            + str(pvs1.get("criteria_strength") or "PVS1")
            + "). Abou Tayoun 2018 forbids PM4 alongside PVS1 at any strength — "
            "both score the same protein-length observation, so applying them "
            "together double-counts one piece of evidence"
        ))

    ps4 = by_code.get("PS4")
    if (
        verify_source_evidence
        and ps4
        and ps4.get("status") == "met"
        and (inheritance or "").upper() == "AR"
        and _criterion_applicable(gene, "PM3")
    ):
        _rec_facts = ps4.get("facts") if isinstance(ps4.get("facts"), dict) else {}
        _rec_or = _rec_facts.get("or_lower_ci")
        has_case_control = (
            isinstance(_rec_or, (int, float)) and not isinstance(_rec_or, bool)
            and _rec_or > 0
        )
        if not has_case_control:
            _demote(
                "PS4",
                f"{gene} is recessive and its VCEP places in-trans evidence "
                "under PM3 — biallelic / in-trans proband occurrences are PM3, "
                "not PS4 case-enrichment (ACMG: no double-counting). PS4 for a "
                "recessive variant requires a case-control odds ratio, which "
                "was not quoted.",
            )

    # CITATION CORRECTED 2026-08-25. This ladder was attributed in-comment (and
    # A clinical tool that mis-attributes its own scoring rule to the guidance a
    ps4 = by_code.get("PS4")
    if verify_source_evidence and ps4 and ps4.get("status") == "met":
        _ps4_spec = _vcep_freq(gene)
        _ps4_is_raso = bool(_ps4_spec) and _ps4_spec.get("vcep") == "RASopathy"
        ev4 = ps4.get("evidence") or ""
        _facts = ps4.get("facts") if isinstance(ps4.get("facts"), dict) else {}
        _fact_or = _facts.get("or_lower_ci")
        _fact_or_present = (
            isinstance(_fact_or, (int, float)) and not isinstance(_fact_or, bool)
            and _fact_or > 0
        )
        _has_or = _fact_or_present
        _ps4_pub = _ps4_spec_occurrence_ladder(gene)
        if _ps4_pub and not _has_or:
            _fact_n = _facts.get("proband_count")
            _from_facts = (
                isinstance(_fact_n, (int, float)) and not isinstance(_fact_n, bool)
                and _fact_n >= 0
            )
            n_prob = int(_fact_n) if _from_facts else 0
            if not _from_facts:
                for _m in re.finditer(
                    r"(\d+)\s+(?:unrelated\s+)?(?:probands?|patients?|individuals?"
                    r"|index\s+cases?|cases?|families|kindreds?|occurrences?)",
                    ev4, re.IGNORECASE,
                ):
                    n_prob = max(n_prob, int(_m.group(1)))
            cap = None
            for _lab in ("Strong", "Moderate", "Supporting"):
                if _lab in _ps4_pub and n_prob >= _ps4_pub[_lab]:
                    cap = _lab
                    break
            if cap is None:
                _demote(
                    "PS4",
                    f"PS4 not applied: {gene} VCEP publishes a minimum of "
                    f">={min(_ps4_pub.values())} independent occurrence(s) and "
                    f"{'the evidence quantifies ' + str(n_prob) if n_prob else 'no count was quantified'}"
                    f". A curator who has counted the occurrences can apply it "
                    f"by hand.",
                )
            else:
                _rank = {"Supporting": 1, "Moderate": 2, "Strong": 4,
                         "Very Strong": 8, "VeryStrong": 8}
                _cur = (ps4.get("criteria_strength") or "PS4").split("_", 1)
                _cur_tier = _cur[1] if len(_cur) > 1 else "Strong"
                if _rank.get(_cur_tier, 4) > _rank[cap]:
                    ps4["criteria_strength"] = f"PS4_{cap}"
                    ps4["evidence"] = (
                        ev4.rstrip()
                        + f" [HeartVar: PS4 held at {cap} — {gene} VCEP "
                          f"publishes PS4_{cap} at >={_ps4_pub[cap]} independent "
                          f"occurrence(s) and the evidence quantifies "
                          f"{n_prob}.]"
                    )
        elif not _ps4_is_raso and not _has_or:
            _fact_n = _facts.get("proband_count")
            _from_facts = (
                isinstance(_fact_n, (int, float)) and not isinstance(_fact_n, bool)
                and _fact_n >= 0
            )
            if _from_facts:
                n_prob = int(_fact_n)
                if "pmid" in ev4.lower():
                    cap = ("Strong" if n_prob >= 15
                           else "Moderate" if n_prob >= 6 else "Supporting")
                else:
                    cap = "Supporting"
            else:
                n_prob = 0
                for _m in re.finditer(
                    r"(\d+)\s+(?:unrelated\s+)?(?:probands?|patients?|individuals?"
                    r"|index\s+cases?|cases?|families|kindreds?)",
                    ev4, re.IGNORECASE,
                ):
                    n_prob = max(n_prob, int(_m.group(1)))
                cap = ("Strong" if n_prob >= 15
                       else "Moderate" if n_prob >= 6 else "Supporting")
            _rank = {"Supporting": 1, "Moderate": 2, "Strong": 4,
                     "Very Strong": 8, "VeryStrong": 8}
            _cur = (ps4.get("criteria_strength") or "PS4").split("_", 1)
            _cur_tier = _cur[1] if len(_cur) > 1 else "Strong"
            if _rank.get(_cur_tier, 4) > _rank[cap]:
                ps4["criteria_strength"] = f"PS4_{cap}"
                _src = "facts.proband_count" if _from_facts else "evidence text"
                ps4["evidence"] = (
                    ev4.rstrip()
                    + f" [HeartVar: PS4 capped to {cap} — {_src} quantifies "
                    + (f"{n_prob} proband(s)" if n_prob else "no proband count")
                    + " and no case-control odds ratio; graded on the Kelly et al."
                    " 2018 proband ladder ≥2/≥6/≥15 = Supporting/Moderate/Strong."
                    + (f" Note: the {gene} VCEP's published rungs are"
                       " odds-ratio bounds, not proband counts, so this"
                       " strength follows panel practice rather than the"
                       " registry text."
                       if _ps4_case_control_only(gene) else "")
                    + "]"
                )

    pp1 = by_code.get("PP1")
    if (
        pp1 and pp1.get("status") == "met"
        and verify_source_evidence
        and not has_curator_segregation
        and pp1.get("source") != "inferred"
        and not _curator_set("PP1")
    ):
        _pf = pp1.get("facts") if isinstance(pp1.get("facts"), dict) else {}
        _pp1_ev = pp1.get("evidence") or ""
        _pf_pmid = _pf.get("pmid")
        _has_pmid = bool(
            (_pf_pmid is not None
             and re.fullmatch(r"\d{4,9}", str(_pf_pmid).strip()))
            or re.search(r"pmid\s*[:#]?\s*\d+|\bpmid\b", _pp1_ev, re.IGNORECASE)
        )
        _carriers = _pf.get("seg_affected_carriers")
        _carriers_n = (
            int(_carriers) if isinstance(_carriers, (int, float))
            and not isinstance(_carriers, bool) else 0
        )
        if not _has_pmid:
            _demote(
                "PP1",
                "PP1 from published literature requires a cited PMID for the "
                "co-segregation report (facts.pmid or a PMID in the evidence); "
                "none was provided — suppressed to prevent inference from "
                "pedigree / AD-label / family-history phrasing.",
            )
        elif _carriers_n < 1:
            _demote(
                "PP1",
                "PP1 from published literature requires ≥1 GENOTYPED affected "
                "carrier (facts.seg_affected_carriers); none was quoted — "
                "affected-but-untested relatives are PP4/phenotype context, not "
                "PP1. Suppressed.",
            )
        else:
            _base = _pp1_strength_from_lod_meioses(
                _pf.get("lod"), _pf.get("seg_meioses"),
                _carriers_n, _pf.get("seg_unaffected_carriers"),
                gene=gene,
            )
            if not _base:
                _demote(
                    "PP1",
                    f"PP1 co-segregation facts do not meet {gene or 'this gene'}'s "
                    f"published minimum (≥{_pp1_spec_ladder(gene)[0]} affected "
                    f"relatives carrying the variant or informative meioses, an "
                    f"LOD ≥0.9, or ≥2 genotyped affected carriers, after any "
                    f"unaffected-carrier temper) — one genotyped carrier alone is "
                    f"not co-segregation. Suppressed.",
                )
            else:
                _rank = {"Supporting": 1, "Moderate": 2, "Strong": 4}
                _inv = {1: "Supporting", 2: "Moderate", 4: "Strong"}
                if _rank[_base] > _rank["Moderate"]:
                    _base = "Moderate"
                _cur = (pp1.get("criteria_strength") or "PP1").split("_", 1)
                _cur_tier = _cur[1] if len(_cur) > 1 else "Supporting"
                _final_rank = min(_rank.get(_cur_tier, 1), _rank[_base])
                _final = _inv[_final_rank]
                _meioses = _pf.get("seg_meioses")
                _lod = _pf.get("lod")
                pp1["criteria_strength"] = f"PP1_{_final}"
                pp1["evidence"] = (
                    _pp1_ev.rstrip()
                    + f" [HeartVar: PP1 {_final} — literature co-segregation "
                    + f"({_carriers_n} genotyped affected carrier(s)"
                    + (f", {int(_meioses)} informative meioses"
                       if isinstance(_meioses, (int, float))
                       and not isinstance(_meioses, bool) and _meioses else "")
                    + (f", LOD {_lod}"
                       if isinstance(_lod, (int, float))
                       and not isinstance(_lod, bool) else "")
                    + "); ClinGen Cardiomyopathy VCEP (Kelly) ladder, "
                    "capped at Moderate for literature-sourced PP1, and never "
                    "raised above the emitted strength.]"
                )

    _ps3_splice_only = False
    _bs3_splice_only = False

    bs3 = by_code.get("BS3")
    if (
        verify_source_evidence
        and bs3 and bs3.get("status") == "met"
        and not _curator_set("BS3")
    ):
        ev = (bs3.get("evidence") or "").lower()
        has_pmid = re.search(r"pmid\s*[:#]?\s*\d+|\bpmid\b", ev)
        has_assay = _FUNCTIONAL_ASSAY_RE.search(ev)
        _bs3_splice_only = _is_splicing_only_assay(ev)
        if _bs3_splice_only:
            _demote(
                "BS3",
                "The only assay cited measures RNA/splicing. ClinGen SVI "
                "(Walker 2023) applies BS3 only to well-established assays of "
                "functional impact NOT directly captured by RNA-splicing "
                "assays — an RNA result showing no splicing impact belongs on "
                "BP7. Suppressed here and considered for BP7.",
            )
        elif not (has_pmid and has_assay):
            _demote(
                "BS3",
                "BS3 requires a variant-specific functional assay with a PMID; "
                "a benign/likely-benign ClinVar record or a frequency signal is "
                "not BS3 evidence (those are BA1/BS1/BP4). No functional assay + "
                "PMID was cited — suppressed.",
            )

    ps3 = by_code.get("PS3")
    if (
        verify_source_evidence
        and ps3 and ps3.get("status") == "met"
        and not _curator_set("PS3")
    ):
        ev = (ps3.get("evidence") or "").lower()
        _psf = ps3.get("facts") if isinstance(ps3.get("facts"), dict) else {}
        _pmids = [
            str(p).strip() for p in (_psf.get("assay_pmids") or [])
            if re.fullmatch(r"\d{4,9}", str(p).strip())
        ]
        _atype = _psf.get("assay_type")
        _type_ok = (
            isinstance(_atype, str)
            and bool(_FUNCTIONAL_ASSAY_RE.search(_atype.lower()))
        )
        _facts_ok = bool(_pmids) and _type_ok and any(p in ev for p in _pmids)
        _ps3_splice_only = _is_splicing_only_assay(
            (_atype or "").lower() if _facts_ok else ev
        )
        if _ps3_splice_only:
            _demote(
                "PS3",
                "The only assay cited measures RNA/splicing. ClinGen SVI "
                "(Walker 2023) applies PS3 only to well-established assays of "
                "functional impact NOT directly captured by RNA-splicing "
                "assays — splicing results belong on PVS1 (loss of function) or "
                "BP7, not PS3. Suppressed here and, where the variant's "
                "consequence allows it, applied to PVS1 instead.",
            )
        elif not _facts_ok:
            has_pmid = re.search(r"pmid\s*[:#]?\s*\d+|\bpmid\b", ev)
            has_assay = _FUNCTIONAL_ASSAY_RE.search(ev)
            if not (has_pmid and has_assay):
                _demote(
                    "PS3",
                    "PS3 requires a variant-specific, well-established functional "
                    "assay with a PMID; an in-silico/ProtVar prediction, a "
                    "same-residue pathogenic neighbour (PM5/PM1), a splice "
                    "prediction (PP3), or a pathogenic ClinVar record is not PS3 "
                    "evidence. No functional assay + PMID was cited — suppressed.",
                )

    _pvs1_route = by_code.get("PVS1")
    if (
        _ps3_splice_only
        and _pvs1_route
        and _pvs1_route.get("_pvs1_rna_routable")
        and _pvs1_route.get("status") != "met"
        and not _curator_set("PVS1")
    ):
        _pvs1_route["status"] = "met"
        _pvs1_route["criteria_strength"] = "PVS1_Strong"
        _pvs1_route["evidence"] = (
            "An RNA/splicing assay experimentally confirms an abnormal "
            "transcript for this variant, which is the 'RNA-confirmed' arm of "
            "the ClinGen PVS1 splice criteria. Per ClinGen SVI (Walker 2023) "
            "splicing assay data belongs on PVS1 rather than PS3, so the "
            "evidence is applied here at PVS1_Strong rather than full PVS1 — "
            "the variant is not at a canonical +/-1,2 site and the transcript's "
            "frame and NMD competence are not established from the assay. "
            f"(was: {str(_pvs1_route.get('evidence') or '')[:150]})"
        )

    _bp7_route = by_code.get("BP7")
    if (
        _bs3_splice_only
        and _bp7_route
        and _bp7_route.get("_bp7_rna_routable")
        and _bp7_route.get("status") != "met"
        and not _curator_set("BP7")
    ):
        _bp7_route["status"] = "met"
        _bp7_route["criteria_strength"] = "BP7_Supporting"
        _bp7_route["evidence"] = (
            "An RNA/splicing assay experimentally shows NO splicing impact for "
            "this synonymous/non-coding variant. Per ClinGen SVI (Walker 2023) "
            "that result belongs on BP7 rather than BS3, so it is applied here. "
            "The nucleotide was also confirmed NOT highly conserved "
            "(phyloP100way <=2.0) — the RNA result substitutes for the "
            "splice-impact test only, not for BP7's conservation requirement. "
            f"(was: {str(_bp7_route.get('evidence') or '')[:150]})"
        )

    ps3_cap = by_code.get("PS3")
    if ps3_cap and ps3_cap.get("status") == "met" and not _curator_set("PS3"):
        _ceiling = _ps3_spec_ceiling(gene)
        _earned, _earned_why = _ps3_validated_strength(ps3_cap.get("facts"))
        _gn3 = (_VCEP_CRIT.get((gene or "").upper()) or {}).get("gn", "CSpec")

        _floor = _ps3_spec_floor(gene) if _earned is None else None
        if _floor is not None:
            _earned = _floor
            _earned_why = (
                f"{gene} VCEP ({_gn3}) publishes the strength to use when an "
                f"assay reports no known variant validation controls: \"if no "
                f"known variant validation controls (i.e. established "
                f"pathogenic and benign variants) were used, then score at the "
                f"{_floor.lower()} strength\". The assay's own validation is "
                f"not described, so that published floor applies rather than "
                f"the ClinGen SVI (Brnich 2019) no-evidence starting point"
            )

        if _earned is None:
            _demote(
                "PS3",
                f"PS3 not applied: {_earned_why}. A curator who has assessed the "
                f"assay can apply it by hand.",
            )
        else:
            _limits = [_earned] + ([_ceiling] if _ceiling else [])
            _final = min(_limits, key=lambda s: _PS3_RANK[f"PS3_{s}"])
            _cur = ps3_cap.get("criteria_strength") or "PS3"
            if _PS3_RANK.get(_cur, 4) > _PS3_RANK[f"PS3_{_final}"]:
                if _ceiling and _final == _ceiling and _ceiling != _earned:
                    _why3 = (
                        f"{gene} VCEP ({_gn3}) publishes no PS3 row above "
                        f"{_ceiling}, so PS3 is held there even though "
                        f"{_earned_why}."
                    )
                else:
                    _why3 = (
                        f"PS3 applied at {_final} because {_earned_why}."
                    )
                ps3_cap["criteria_strength"] = f"PS3_{_final}"
                ps3_cap["evidence"] = (
                    f"{_why3} (was: {str(ps3_cap.get('evidence') or '')[:140]})"
                )


    pp1_seg = by_code.get("PP1")
    bs4_seg = by_code.get("BS4")
    if (
        pp1_seg and bs4_seg
        and pp1_seg.get("status") == "met"
        and bs4_seg.get("status") == "met"
    ):
        _demote(
            "PP1",
            "Contradictory segregation evidence — BS4 (lack of segregation in "
            "affected family members) was also applied. PP1 and BS4 are "
            "opposite readings of one pedigree and cannot both hold; both are "
            "withheld and the conflict is surfaced for curator review.",
        )
        _demote(
            "BS4",
            "Contradictory segregation evidence — PP1 (cosegregation with "
            "disease in affected family members) was also applied. BS4 and PP1 "
            "are opposite readings of one pedigree and cannot both hold; both "
            "are withheld and the conflict is surfaced for curator review.",
        )

    def _is_strong(c: dict, codes: tuple[str, ...]) -> bool:
        if c.get("status") != "met" or c["code"] not in codes:
            return False
        return not str(c.get("criteria_strength") or "").endswith(
            ("Moderate", "Supporting")
        )

    bp6 = by_code.get("BP6")
    if bp6 and bp6.get("status") == "met":
        strong_path = next(
            (c for c in criteria if _is_strong(c, ("PVS1", "PS1", "PS3", "PS4"))),
            None,
        )
        if strong_path is not None:
            _demote(
                "BP6",
                "BP6 (benign ClinVar assertion) withheld: HeartVar fired "
                f"{strong_path['code']} (strong pathogenic evidence) "
                "independently — the engine's pathogenic evidence overrides "
                "the benign ClinVar assertion. Conflict surfaced for curator "
                "review.",
            )

    pp5 = by_code.get("PP5")
    if pp5 and pp5.get("status") == "met":
        strong_benign = next(
            (c for c in criteria if _is_strong(c, ("BS2", "BS3", "BS4"))),
            None,
        )
        if strong_benign is not None:
            _demote(
                "PP5",
                "PP5 (pathogenic ClinVar assertion) withheld: HeartVar fired "
                f"{strong_benign['code']} (strong benign evidence) "
                "independently — the engine's benign evidence overrides the "
                "pathogenic ClinVar assertion. Conflict surfaced for curator "
                "review.",
            )

    return criteria, forced_classification


def merge_hard_coded_and_ai(
    hard_coded: list[dict], ai_criteria: list[dict],
) -> list[dict]:
    """Combine the hard-coded set with the AI-evaluated set, preserving
    the canonical 28-code order so the frontend renderer keeps stable
    ordering. Hard-coded entries always win — if the AI accidentally
    returns a code we own, we drop the AI version. Missing entries are
    filled with a placeholder not_met record so downstream consumers
    always see the full 28; the placeholder's wording depends on WHO was
    silent (see below).

    ``hard_coded`` is the deterministic set plus, since 2026-09-08, the
    PM1/PP1/BS4 entries from no_ai.infer_supplementary_criteria, which the
    caller concatenates on (see app._python_authoritative_supplementary)."""
    by_code: dict[str, dict] = {c["code"]: c for c in hard_coded}
    for c in ai_criteria:
        code = c.get("code") or ""
        if not code or code in by_code:
            continue
        by_code[code] = c
    ordered: list[dict] = []
    for code in CANONICAL_CRITERIA_ORDER:
        entry = by_code.get(code)
        if entry is None:
            _server_owned = code in HARD_CODED_CRITERIA_CODES
            entry = {
                "code": code,
                "name": _CRITERION_NAMES.get(code, code),
                "status": "not_met",
                "direction": "benign" if code.startswith("B") else "pathogenic",
                "criteria_strength": None,
                "evidence": (
                    "Evaluated deterministically by HeartVar. The rule did "
                    "not fire on the inputs supplied."
                    if _server_owned else
                    "AI did not return this criterion — defaulted to not_met"
                ),
                "source": "hard_coded" if _server_owned else "ai",
            }
        ordered.append(entry)
    return ordered


_RELATIVE_WORDS = (
    r"relative|sibling|brother|sister|cousin|aunt|uncle|nephew|niece|"
    r"grandparent|grandfather|grandmother|grandson|granddaughter|"
    r"parent|father|mother|child|son|daughter|twin|proband's\s+\w+"
)
_AFFECTION_WORDS = (
    r"affected|diagnos\w+|symptomatic|died|deceased|passed\s+away|"
    r"sudden\s+death|transplant\w*"
)

_FAMILY_HISTORY_NEGATIVE_RE = re.compile(
    r"\b("
    r"no family history|no\s+(?:known\s+)?family history|"
    r"unaffected parents?|parents?\s+(?:are\s+)?unaffected|"
    r"no\s+relatives?\s+affected|no\s+affected\s+relatives?|"
    r"no\s+history\s+of|nil|n/?a|"
    r"none\s+reported|negative\s+family\s+history|"
    rf"no\s+(?:\w+\s+){{0,2}}(?:{_RELATIVE_WORDS})s?\s+(?:are\s+|is\s+|were\s+|was\s+)?affected|"
    rf"(?:{_RELATIVE_WORDS})s?\s+(?:are\s+|is\s+|were\s+|was\s+|both\s+)?unaffected|"
    r"denies\s+(?:any\s+)?family\s+history|"
    r"no\s+other\s+affected|"
    r"(?:apparently\s+)?sporadic|"
    r"family\s+history\s+(?:is\s+)?(?:negative|unremarkable|non[-\s]?contributory)"
    r")\b",
    re.IGNORECASE,
)
_FAMILY_HISTORY_SEGREGATION_RE = re.compile(
    r"\b("
    r"co[-\s]?segregat\w*|segregat\w+\s+(?:with|in)|"
    r"multiple\s+generations?\s+tested|"
    r"variant\s+confirm\w+\s+in\s+(?:affected\s+)?relatives?|"
    r"tested\s+positive\s+in|carry\s+the\s+same\s+variant|"
    r"carries\s+the\s+same\s+variant|"
    r"both\s+(?:parents|relatives)\s+(?:tested|carriers?)|"
    r"test(?:ed|s|ing)?\s+positive|"
    r"informative\s+meios\w*|lod\s*(?:score|=|:|of|\u2265|>)|"
    r"also\s+(?:carries|carry|carrie[sd]|a\s+carrier)|"
    r"carr(?:y|ies|ied)\s+(?:the\s+)?(?:same\s+|this\s+)?variant|"
    rf"variant\s+(?:was\s+)?(?:also\s+)?(?:found|detected|identified|present)\s+in\s+(?:the\s+|an?\s+)?(?:affected\s+)?(?:{_RELATIVE_WORDS})s?|"
    r"genotyp\w+|"
    r"segregat\w+"
    r")\b",
    re.IGNORECASE,
)
_FAMILY_HISTORY_POSITIVE_RE = re.compile(
    r"\b("
    r"affected\s+(?:relative|sibling|brother|sister|cousin|aunt|uncle|"
    r"grandparent|grandfather|grandmother|parent|father|mother|child|son|daughter)s?|"
    r"(?:father|mother|brother|sister|cousin|aunt|uncle|sibling|son|daughter|"
    r"grandparent|grandfather|grandmother|child|parent)s?\s+(?:has|had|with|"
    r"diagnosed|carries|carrier)|"
    r"family\s+history\s+of|family\s+history\s+notable|"
    r"(?:positive|known)\s+family\s+history|"
    r"obligate\s+carriers?|known\s+carriers?|carrier\s+(?:father|mother|parent)"
    r")\b",
    re.IGNORECASE,
)


_FAMILY_HISTORY_RELATIVE_RE = re.compile(rf"\b(?:{_RELATIVE_WORDS})s?\b", re.IGNORECASE)
_FAMILY_HISTORY_AFFECTION_RE = re.compile(rf"\b(?:{_AFFECTION_WORDS})\b", re.IGNORECASE)


def classify_family_history(family: str | None) -> str:
    """Categorise the free-text family-history field into one of four
    discrete buckets so the Summary tab's Query Details panel can show
    a short label instead of the raw paragraph.

    Returns one of: "no_family_history", "positive_family_history",
    "segregation_data", "unknown". The full free-text is still passed
    through to the AI prompt unchanged — this classifier exists purely
    for the display row.

    Detection order is intentional: segregation evidence is the most
    specific signal, so it wins over a generic "positive family history"
    when both are present (e.g. "Father carries the same variant"
    matches segregation, not just positive). Negation phrases are
    checked first so "no family history" can't be miscategorised as
    "positive" by a stray "history of" later in the sentence.
    """
    text = (family or "").strip()
    if not text:
        return "unknown"
    if _FAMILY_HISTORY_NEGATIVE_RE.search(text):
        return "no_family_history"
    if _FAMILY_HISTORY_SEGREGATION_RE.search(text):
        return "segregation_data"
    if _FAMILY_HISTORY_POSITIVE_RE.search(text):
        return "positive_family_history"
    if _FAMILY_HISTORY_RELATIVE_RE.search(text) and \
            _FAMILY_HISTORY_AFFECTION_RE.search(text):
        return "positive_family_history"
    return "unknown"


def _hgvsp_int_position(hgvsp: str) -> int | None:
    """Parse the integer residue position out of a VEP ``hgvsp`` suffix
    (e.g. ``NP_004324.2:p.Phe595Leu`` → 595). Returns None when no digits
    are present. The two position extractors below apply their own
    consequence gates first and share only this final parse step."""
    m = re.search(r"\d+", (hgvsp or "").split(":")[-1])
    return int(m.group(0)) if m else None


_AA3_TO_AA1: dict[str, str] = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Sec": "U", "Pyl": "O",
}


def _hgvsp_ref_alt_aa(hgvsp: str | None) -> tuple[str | None, str | None]:
    """Reference and alternate amino acids from an hgvsp, as one-letter codes.

    ``NP_000129.3:p.Cys570Arg`` -> ("C", "R"). Returns (None, None) for anything
    that is not a simple single-residue substitution — frameshifts, in-frame
    indels, ``p.Cys570=`` — because the identity rules that call this are only
    defined for substitutions, and guessing on a frameshift would invent one.
    """
    tail = (hgvsp or "").split(":")[-1]
    m = re.fullmatch(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", tail.strip())
    if not m:
        return None, None
    return _AA3_TO_AA1.get(m.group(1)), _AA3_TO_AA1.get(m.group(3))


def _any_protein_position_from_vep(vep: dict) -> int | None:
    """Position extractor for PM1's domain query — accepts any
    protein-altering variant, not just missense.

    Differs from ``_protein_position_from_vep`` (which is PM5-scoped and
    missense-only) because PM1 applies to any variant inside a critical
    functional domain: stop-gain, in-frame indel, missense, etc. We
    skip non-protein consequences (splice, intronic, synonymous, 5'/3'
    UTR) where there's no defined residue.
    """
    if not isinstance(vep, dict) or not vep.get("ok"):
        return None
    consequence = (vep.get("most_severe_consequence") or "").lower()
    if any(skip in consequence for skip in (
        "synonymous", "splice_region", "intron", "intergenic",
        "upstream", "downstream", "5_prime_utr", "3_prime_utr",
        "non_coding", "regulatory",
    )):
        return None
    hgvsp = vep.get("hgvsp") or ""
    if not hgvsp:
        return None
    return _hgvsp_int_position(hgvsp)


def _mane_hgvsp_from_vep(vep: dict) -> str | None:
    """Return the ``hgvsp`` of the MANE Select transcript in a VEP result.

    WHY THIS EXISTS. Every residue number HeartVar compares against ClinVar
    has to be in ClinVar's own frame of reference. ClinVar writes its ``Name``
    field on its preferred transcript, which is the MANE Select one — so
    ``NM_001308093.3(GATA4):c.910G>A (p.Gly304Arg)``. But
    ``_pick_transcript_consequence`` deliberately prefers the transcript the
    CURATOR supplied, so ``vep["hgvsp"]`` for the same variant entered as
    ``NM_002052.5:c.907G>A`` is ``NP_002043.2:p.Gly303Arg``. NP_002043.2 is
    442 aa and NP_001295022.1 is 443 aa (one extra Val at position 206), so
    every residue past 205 is numbered +1 in the MANE isoform: residue 303 of
    the supplied transcript and residue 304 of ClinVar's are the SAME residue.

    Comparing the supplied-transcript integer against ClinVar's integer
    therefore both MISSED genuine same-residue evidence (PM5 false negative)
    and MATCHED a different residue that happened to carry the same integer
    (PS1/PM5 false positive). The MANE number is the one to match on; the
    curator's number is the one to display.

    ``transcript_consequences_all`` is built by ``_build_transcript_table``
    and is populated by both the live and the offline VEP clients, so the
    MANE row is already in memory wherever a VEP result is. Falls back to
    MANE Plus Clinical, then to the picked transcript when it is itself MANE.
    Returns None when no MANE annotation is available — callers must then
    keep using the supplied-transcript number rather than guess.
    """
    if not isinstance(vep, dict):
        return None
    rows = vep.get("transcript_consequences_all") or []
    for flag in ("is_mane_select", "is_mane_plus_clinical"):
        for row in rows:
            if isinstance(row, dict) and row.get(flag) and row.get("hgvsp"):
                return row["hgvsp"]
    if vep.get("is_mane_select") or vep.get("is_mane_clinical"):
        return vep.get("hgvsp") or None
    return None


def _mane_hgvsc_from_vep(vep: dict) -> str | None:
    """Bare ``c.…`` token of the MANE Select transcript in a VEP result.

    Companion to ``_mane_hgvsp_from_vep``, for the ClinVar variant-level
    lookup rather than the residue-level one. ClinVar's ``Name`` carries the
    coding change of its preferred (MANE Select) transcript, so this is the
    token that actually matches: GATA4 entered as ``NM_002052.5:c.886G>A``
    has to be looked up as ``c.889G>A``. Returns the bare token (no
    transcript prefix) because ``_hgvs_match`` matches inside the Name.
    """
    if not isinstance(vep, dict):
        return None
    hgvsc = None
    rows = vep.get("transcript_consequences_all") or []
    for flag in ("is_mane_select", "is_mane_plus_clinical"):
        for row in rows:
            if isinstance(row, dict) and row.get(flag) and row.get("hgvsc"):
                hgvsc = row["hgvsc"]
                break
        if hgvsc:
            break
    if not hgvsc and (vep.get("is_mane_select") or vep.get("is_mane_clinical")):
        hgvsc = vep.get("hgvsc")
    if not hgvsc:
        return None
    token = str(hgvsc).split(":")[-1].strip()
    return token or None


def _any_mane_protein_position_from_vep(vep: dict) -> int | None:
    """MANE-transcript counterpart of ``_any_protein_position_from_vep`` —
    the residue number to MATCH ClinVar Names against for PM1's domain query.
    Applies the same non-protein-consequence gate; returns None when the
    variant has no protein consequence or when no MANE annotation exists."""
    if _any_protein_position_from_vep(vep) is None:
        return None
    return _hgvsp_int_position(_mane_hgvsp_from_vep(vep) or "")


def _pm1_assessment_from_domain_plp(
    payload: dict, most_severe: str | None = None,
) -> tuple[str, str]:
    """Derive a PM1 verdict + one-sentence summary from the domain-P/LP
    payload returned by ``get_domain_plp_evidence``.

    Returns ``(status, sentence)`` where status is one of
    ``met / supporting / not_met / not_applicable``. Thresholds come
    straight from the user spec:
      - ≥3 P/LP in domain OR any ≥2-star P/LP → met (Moderate, +2)
      - 1-2 P/LP                              → supporting (Supporting, +1)
      - 0 P/LP                                → not_met
      - variant outside any annotated domain  → not_applicable

    ``most_severe`` is the variant's VEP most-severe consequence. When it is a
    loss-of-function consequence this verdict is demoted to ``not_met`` — a
    hard mirror of the ``_PM1_LOF_CONSEQUENCES`` guard in
    ``_gate_criteria_applicability`` (see lines above). Without this, the
    informational "Domain pathogenicity context" panel would claim "PM1 met"
    on a LoF variant the final criteria list (and the chat) correctly
    suppress — the exact JAG1 p.Arg235* contradiction this reconciles.
    """
    if not isinstance(payload, dict) or not payload.get("ok"):
        return "not_applicable", "PM1 not assessed (domain lookup unavailable)."
    if payload.get("not_applicable"):
        pre = payload.get("pm1_assessment")
        if isinstance(pre, str) and pre.strip():
            return "not_applicable", pre
        reason = (payload.get("reason") or "").lower()
        if "lookup unavailable" in reason:
            sentence = "PM1 not assessed — UniProt lookup unavailable."
        elif "no protein position" in reason or "non-coding" in reason:
            sentence = (
                "Variant has no protein consequence (splice / non-coding) — "
                "PM1 not assessed from domain evidence."
            )
        else:
            sentence = (
                "Variant not within an annotated UniProt domain — "
                "PM1 evidence not applicable."
            )
        return "not_applicable", sentence
    cons = (most_severe or "").lower()
    if cons in _PM1_LOF_CONSEQUENCES:
        dom_name = payload.get("domain_name") or "the containing domain"
        total_lof = int(payload.get("total_plp") or 0)
        ctx = (
            f" ({dom_name} does contain {total_lof} P/LP variant"
            f"{'s' if total_lof != 1 else ''}, but that reflects other "
            f"loss-of-function variants scored by PVS1, not a missense hotspot)"
            if total_lof else ""
        )
        return (
            "not_met",
            f"PM1 not applied — this variant is {cons} (loss-of-function, "
            f"covered by PVS1); PM1 is a missense / in-frame criterion{ctx}.",
        )
    dom = payload.get("domain_name") or "the containing domain"
    span = (
        f"aa {payload.get('domain_start')}-{payload.get('domain_end')}"
        if payload.get("domain_start") is not None else ""
    )
    total = int(payload.get("total_plp") or 0)
    has_2star = bool(payload.get("has_two_star_plp"))
    pieces_loc = f"{dom}" + (f" ({span})" if span else "")
    if total == 0:
        return (
            "not_met",
            f"0 P/LP variants in {pieces_loc} — domain is not an established "
            "mutational hotspot.",
        )
    star_clause = (
        f"; {1 if has_2star else 0}+ with ≥2-star review"
        if has_2star else " (no ≥2-star records)"
    )
    if total >= 3 or has_2star:
        return (
            "met",
            f"PM1 met: {total} P/LP variants in {pieces_loc}{star_clause}.",
        )
    return (
        "supporting",
        f"PM1 supporting (downgraded to Supporting): {total} P/LP variant"
        f"{'s' if total != 1 else ''} in {pieces_loc}.",
    )
