"""ACMG point summation and tier classification.

Pure functions that turn a criteria list into a signed points_total and the
Tavtigian-2020 tier label. Drives off the constant tables in
``backend.acmg.constants`` so the cutoffs stay single-sourced with the browser.
"""
from __future__ import annotations

from .constants import (
    _TIER_POINTS,
    _BARE_CODE_POINTS,
    _CLASSIFICATION_THRESHOLDS,
)


_TIER_LOOKUP = {k.replace(" ", "").lower(): v for k, v in _TIER_POINTS.items()}
_BARE_LOOKUP = {k.replace(" ", "").upper(): v for k, v in _BARE_CODE_POINTS.items()}


def _norm_tier(tier: str) -> str:
    """Fold a strength tier token to its _TIER_LOOKUP key: "Very Strong",
    "Very_Strong" and "very strong" all become "verystrong"."""
    return tier.replace(" ", "").replace("_", "").lower()


def _points_for(
    criteria_strength: str | None, code: str, curator_override: bool = False,
) -> int:
    """Map a criteria_strength string ("PVS1_Strong", "PM2_Supporting",
    "BA1", "PP3", ...) to its signed point value per Tavtigian 2020.

    PP5/BP6 ("reputable source") contribute ZERO points *when the engine applies
    them*. The ClinGen SVI Working Group recommended these criteria be retired
    (Biesecker & Harrison, Genet Med 2018) — they cannot be applied reproducibly
    and import an unverifiable third-party classification (a circularity when
    the reputable source is a ClinVar expert-panel record that is itself the
    reference call). No ClinGen VCEP specification uses them. The PP5/BP6
    entries are still surfaced to the curator (status "met", with the ClinVar
    assertion string) as *context*, but they do not move the point total or the
    classification.

    ``curator_override`` lifts that zeroing, and ONLY that zeroing. It is set
    exclusively by the Criteria-tab override (POST /api/acmg/rescore) when a
    human has deliberately switched a criterion on and chosen its strength. The
    distinction is the point: the engine silently reusing ClinVar's verdict is
    circular, whereas a curator electing to count it is a documented judgement
    they own and the report records. Nothing the engine produces carries the
    flag, so the automatic path is unchanged — see
    tests/test_acmg_rescore_endpoint.py.
    """
    if (code or "").strip().upper() in ("PP5", "BP6") and not curator_override:
        return 0
    if not criteria_strength:
        return 0
    s = criteria_strength.strip()
    base, sep, tier = s.partition("_")
    direction_sign = -1 if base[:1].upper() == "B" else 1
    if sep == "_":
        pts = _TIER_LOOKUP.get(_norm_tier(tier))
        if pts is not None:
            return direction_sign * pts
    return _BARE_LOOKUP.get(s.replace(" ", "").upper(), 0)


def compute_points_total(criteria: list[dict]) -> int:
    """Sum signed points across met criteria using their criteria_strength."""
    total = 0
    for c in criteria or []:
        if c.get("status") != "met":
            continue
        total += _points_for(
            c.get("criteria_strength"), c.get("code", ""),
            bool(c.get("curator_override")),
        )
    return total


def classification_for(points: int) -> str:
    for band in _CLASSIFICATION_THRESHOLDS:
        if points >= band["min"]:
            return band["label"]
    return "Benign"


_BENIGN_STRONG_POINTS = -4


def apply_benign_combining_floor(tier: str, criteria: list[dict] | None) -> str:
    """ACMG/AMP-2015 combining-rule floor on the benign side: demote a
    "Likely benign" that rests on a SINGLE sub-Strong benign criterion back to
    "VUS". Every other tier is returned unchanged.

    The Tavtigian-2020 bands this engine scores against put Likely benign at
    -1 to -6, so one Supporting-strength criterion (-1, in practice almost
    always BP4 off an in-silico predictor) was enough to leave VUS. ACMG/AMP
    2015 (Richards, Genet Med 2015, Table 5) requires 1 Strong + 1 Supporting
    or >= 2 Supporting for Likely benign, so a lone BP4 is a VUS there. This is
    the best-known divergence between the verbal rules and the point system and
    it runs in the false-benign direction, so the stricter rule wins.

    Criteria worth zero points (PP5/BP6, retired per SVI 2018) are not counted
    — BP4+BP6 is still a single line of scoring benign evidence. A BP6 a curator
    has deliberately switched on DOES score, so it counts here too: once a human
    owns it, it is a second line of benign evidence in exactly the sense the
    2015 combining rules mean.

    DELIBERATE CARVE-OUT: a lone Strong benign criterion (BS1 alone, -4) keeps
    its point-system tier of Likely benign. Strict 2015 wants a Supporting
    alongside it, but a strong frequency argument is real benign evidence and
    demoting it is a much larger calibration change than this guard is for."""
    if tier != "Likely benign":
        return tier
    benign_points = []
    for c in criteria or []:
        if c.get("status") != "met":
            continue
        pts = _points_for(
            c.get("criteria_strength"), c.get("code", ""),
            bool(c.get("curator_override")),
        )
        if pts < 0:
            benign_points.append(pts)
    if len(benign_points) == 1 and benign_points[0] > _BENIGN_STRONG_POINTS:
        return "VUS"
    return tier


_GV_UNCAPPED = frozenset({"Definitive", "Strong"})


def apply_gene_validity_ceiling(tier: str, validity: str | None) -> str:
    """Cap a Pathogenic call at Likely pathogenic on a weakly-validated gene.

    ClinGen SVI does not recommend classifying variants beyond Likely pathogenic
    where the gene-disease relationship is Moderate or weaker. The RASopathy
    VCEP states it on MRAS c.67G>C:

        "Given the Moderate strength of gene-disease relationship between MRAS
         and autosomal dominant RASopathy, ClinGen's sequence variant
         interpretation working group does not recommend the classification of
         variants in this gene beyond likely pathogenic."

    ``validity`` is the governing tier from
    ``hard_coded.gene_validity_ceiling(gene)``; pass None when it is unknown and
    nothing happens. It is a CEILING, so only "Pathogenic" can move — LP stays
    LP and the benign side is never touched. Fails open on every unrecognised
    input, because a spurious demotion is the dangerous direction.
    """
    if tier != "Pathogenic":
        return tier
    v = (validity or "").strip()
    if not v or v in _GV_UNCAPPED:
        return tier
    return "Likely pathogenic"


def classification_for_criteria(
    points: int,
    criteria: list[dict] | None,
    forced_classification: str | None = None,
    validity: str | None = None,
) -> str:
    """The tier for a scored criteria set: a stand-alone forced call (BA1 →
    Benign) if there is one, else the Tavtigian band with the benign combining
    floor and the gene-validity ceiling applied. This — not bare
    ``classification_for`` — is what the request paths should call, so neither
    adjustment can be forgotten at a call site.

    ``validity`` is the gene's ClinGen Gene-Disease Validity tier. Omit it and
    the ceiling is inert, so every existing caller keeps its behaviour."""
    if forced_classification is not None:
        return forced_classification
    tier = apply_benign_combining_floor(classification_for(points), criteria)
    return apply_gene_validity_ceiling(tier, validity)


def tier_without_clinvar_assertion(
    criteria: list[dict], forced_classification: str | None = None,
    validity: str | None = None,
) -> tuple[int, str]:
    """The (points, tier) the engine reaches WITHOUT the PP5/BP6 ClinVar
    reputational assertion. This is the honest proxy for novel-variant
    performance (where no ClinVar record exists), reported alongside the
    with-PP5 tier so the benchmark can show what tiered PP5/BP6 actually
    buys. A forced_classification (e.g. BA1→Benign) still applies."""
    filtered = [c for c in criteria if c.get("code") not in ("PP5", "BP6")]
    pts = compute_points_total(filtered)
    tier = classification_for_criteria(pts, filtered, forced_classification, validity)
    return pts, tier
