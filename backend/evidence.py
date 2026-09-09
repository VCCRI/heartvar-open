"""Evidence-gathering layer for HeartVar.

Extracted verbatim from backend/app.py (the SSE evidence orchestrator
``_gather_evidence_sse`` and its supporting helpers). This module owns the
parallel DB-gather phase + the variant-input → VEP resolution + the ClinVar
landscape keyword builder + the PM1 domain-PLP composer. It is imported BACK
into backend.app (app.py does ``from .evidence import ...``) — the dependency
direction is app → evidence, NEVER evidence → app, so there is no cycle.

Pure code-move: no logic was changed relative to the pre-split app.py.
"""
from __future__ import annotations

import asyncio
import json
import logging
import os
import re
from time import perf_counter

from .models import CurationRequest
from .acmg.hard_coded import (
    _pm1_assessment_from_domain_plp,
    _any_mane_protein_position_from_vep,
    _any_protein_position_from_vep,
    _hgvsp_int_position,
    _mane_hgvsc_from_vep,
    _mane_hgvsp_from_vep,
    _next_inframe_met,
)
from .clients.alphafold import fetch_alphafold
from .clients.alphamissense import fetch_alphamissense
from .clients.biogrid import fetch_biogrid
from .clients.chdgene import lookup as chdgene_lookup
from .clients.clinvar import (
    fetch_clinvar,
    get_domain_plp_evidence,
    get_gene_phenotype_strings,
    get_gene_variant_landscape,
    get_pm5_evidence,
    get_same_residue_records,
)
from .clients.ensembl_vep import (
    codon_genomic_positions,
    fetch_transcript_exons,
    fetch_variant_recoder_rsid,
    fetch_vep,
    fetch_vep_by_coordinates,
    gnomad_variant_id_with_provenance,
    liftover_grch37_to_grch38,
    parse_variant_input,
    recode_rsid_to_coords,
)
from .clients.gencc import fetch_gencc
from .clients.erepo_client import fetch_erepo
from .clients.hgnc_alias import canonicalise_gene_symbol
from .clients.gnomad import fetch_gnomad, fetch_same_site_frequencies
from .clients.spliceai import fetch_spliceai
from .clients.fetal_heart import fetch_fetal_heart
from .clients.gtex import fetch_gtex
from .clients.hpo_labels import resolved_token_map
from .clients.medgen import fetch_medgen
from .clients.mgi import fetch_mgi
from .clients.opentargets_client import fetch_opentargets
from .clients.panelapp import audit_phenotype_terms, fetch_panelapp
from .clients.pmcoa import fetch_pmc_fulltext
from .clients.protvar import fetch_protvar
from .clients.pubmed import (
    fetch_pubmed,
    get_gene_disease_literature,
)
from .clients.pubtator3 import fetch_pubtator3
from .clients.uniprot import fetch_uniprot, process_natural_variants
from .task_timing import track_task_timing

log = logging.getLogger("heartvar.evidence")

from .logredact import add_values as _log_add_values  # noqa: E402


def _sse(event: str, data: dict) -> str:
    """Encode one SSE message: event name + JSON payload."""
    return f"event: {event}\ndata: {json.dumps(data)}\n\n"


async def _fetch_clinvar_with_vep_coords(
    vep_task: asyncio.Task, gene: str, hgvs_c: str,
) -> dict:
    """Variant-level ClinVar lookup, chained on the VEP task so it can be
    matched in ClinVar's own frame of reference.

    ClinVar writes its ``Name`` on its preferred (MANE Select) transcript.
    Matching only the curator's bare ``c.`` token therefore reported a
    variant as absent from ClinVar whenever the curator used a non-MANE
    transcript with different numbering — GATA4 ``NM_002052.5:c.886G>A`` is
    ClinVar's ``NM_001308093.3:c.889G>A``, VCV000009030, Pathogenic/Likely
    pathogenic at 2★, and HeartVar said "not found", silently removing PP5.

    Awaiting VEP costs nothing in wall-clock (the ClinVar query is a local
    ~14 ms SQLite read and VEP dominates), but it does mean the ClinVar
    ``db_done`` event now arrives after VEP's rather than in parallel. VEP is
    already the lynchpin — a VEP failure aborts the curation — so nothing is
    lost when it fails: we fall back to the curator's token alone, i.e. the
    previous behaviour exactly.
    """
    # with vep at 1.3-3.0s that is a material share of a 16s figure attributed
    _t_wait = perf_counter()
    try:
        vep_data = await vep_task
    except Exception:  # noqa: BLE001 - VEP failure is handled by the caller
        vep_data = None
    _d_wait = perf_counter() - _t_wait
    mane_c = None
    chrom = pos = end = None
    if isinstance(vep_data, dict) and vep_data.get("ok"):
        mane_c = _mane_hgvsc_from_vep(vep_data)
        chrom = vep_data.get("seq_region_name")
        pos = vep_data.get("start")
        _log_add_values(
            pos, mane_c, vep_data.get("id"),
            f"{chrom}-{pos}" if chrom and pos else None,
            f"{chrom}:{pos}" if chrom and pos else None,
        )
        end = vep_data.get("end")
    _t_call = perf_counter()
    result = await fetch_clinvar(
        gene, hgvs_c, hgvs_c_mane=mane_c, chrom=chrom, pos=pos, end=end,
    )
    log.info(
        "[timing] clinvar-task %s: waited_for_vep=%.3fs then lookup=%.3fs",
        gene, _d_wait, perf_counter() - _t_call,
    )
    return result


def _recoder_can_resolve(hgvs_c: str) -> bool:
    """Whether Ensembl's variant_recoder could possibly resolve this string.

    ⚠ IT COULD NOT, AND IT WAS ASKED ON EVERY CURATION. The rsID booster below
    called ``fetch_variant_recoder_rsid(hgvs_c)``, and by that point
    ``hgvs_c`` has had its transcript stripped off — evidence.py logs
    "Stripped transcript prefix NM_005159.5 from HGVS input (raw=… → c.=…)"
    and keeps the accession separately. So the request that actually went out
    was ``/variant_recoder/human/c.301G%3EA``: a bare coding change with no
    transcript and no gene, which names no locus in any genome and cannot be
    recoded by anything.

    Ensembl answered 500/503, the call burned its whole 8 s retry deadline
    across three attempts, and it did so INSIDE the pubtator3 task — the single
    largest contributor to gather time (measured 2026-08-31: pubtator3 true
    13.11 s on ACTC1, of which 6.5 s was this; 38.02 s on MYH7). A guaranteed
    failure on the critical path of every curation.

    An rsID is explicitly a best-effort retrieval booster and never
    load-bearing: PubTator3 still queries on HGVS and AA terms. So the fix is
    to stop asking, not to repair the argument — repairing it would restore a
    live Ensembl dependency that HEARTVAR_VEP_OFFLINE exists to remove.

    NOTE ON RECALL. Dropping the rsID arm narrows PubTator3's variant matching
    a little (the "PS4 FIND gap" the task comment describes). The proper fix is
    to serve the rsID locally: offline VEP could emit it via ``--check_existing``
    if the data build kept the variation cache (``VEP_CACHE_KEEP_VARIATION=1``,
    see vep_offline.extra_args_needing_variation_cache), or the ClinVar build
    could keep variant_summary's ``RS# (dbSNP)`` column, which it currently
    drops. Either is a data-build change; this is not.
    """
    s = (hgvs_c or "").strip()
    if not s:
        return False
    if s.lower().startswith("rs"):
        return True
    return ":" in s


async def _fetch_literature_with_vep_aa(
    fetch_fn, vep_task: asyncio.Task, gene: str, hgvs_c: str,
    supplied_aa: str | None, pass_rsid: bool = False,
) -> dict:
    """Run a variant-literature client (``fetch_pubmed`` / ``fetch_pubtator3``),
    deriving the amino-acid token from VEP when the curator didn't supply one.

    The web form doesn't expose an amino-acid input (only the benchmark harness
    does), so the variant-specific pre-fetch used to only fire the strictest
    query (HGVS-c + pathogenic gate) and miss papers that cite the protein
    change (e.g. "Arg235Ter"). Chaining on the VEP task lets us pull the AA
    token from ``vep.hgvsp`` and run the broader queries — for PubMed this
    aligns the displayed list with the curator's "Search PubMed for this
    variant" link; for PubTator3 the one-letter form derived inside the client
    is the highest-value term the strict [tiab] query structurally can't match.

    When ``pass_rsid`` is set (PubTator3 ONLY — PubMed's exact-match [tiab]
    display contract must stay rsID-free), the dbSNP rsID is threaded into the
    client as its #1 tmVar3 entity term. Retrieval quick-win Q1: this arg used
    to be silently dropped — the call passed only (gene, hgvs_c, aa), so
    ``fetch_pubtator3``'s ``rsid`` param was always None and the single most
    reliable term for recurrent founder variants never fired. The rsID is taken
    from VEP's colocated-variant record when present (free — already fetched)
    and falls back to one Ensembl variant_recoder lookup otherwise.
    """
    aa = (supplied_aa or "").strip()
    rsid: str | None = None
    if not aa or pass_rsid:
        try:
            vep_data = await vep_task
        except Exception:
            vep_data = None
        if isinstance(vep_data, dict):
            if not aa:
                hgvsp = (vep_data.get("hgvsp") or "").strip()
                if ":p." in hgvsp:
                    aa = hgvsp.split(":p.", 1)[1]
                elif hgvsp.startswith("p."):
                    aa = hgvsp[2:]
            if pass_rsid:
                rsid = (vep_data.get("dbsnp_rsid") or "").strip() or None
    if pass_rsid:
        if not rsid and _recoder_can_resolve(hgvs_c):
            try:
                recoded = await fetch_variant_recoder_rsid(hgvs_c)
                rsid = (recoded or {}).get("rsid") or None
            except Exception:
                rsid = None
        return await fetch_fn(gene, hgvs_c, aa or None, rsid=rsid)
    return await fetch_fn(gene, hgvs_c, aa or None)


_PAREN_RE = re.compile(r"\s*\([^)]*\)")
_WITH_OR_WITHOUT_RE = re.compile(r"\s*,?\s*with or without\b.*")
_TYPE_SUFFIX_RE = re.compile(r"\s*,?\s*type[s]?\s+\w+\s*$")
_GENE_SUFFIX_RE = re.compile(r"\s*[-–]\s*\w+\d+\s*$")
_COMPOUND_SPLIT_RE = re.compile(r"\s+and\s+|\s*,\s*")
_MIN_KEYWORD_LEN = 4

CARDIAC_ROOT_TERMS: tuple[str, ...] = (
    "cardiac", "cardio", "heart", "aort", "valv", "coronary",
    "arrhythmi", "atrial", "ventricular", "atrioventricular",
    "septal", "conotruncal", "congenital", "vascular", "myopath",
    "cardiomyopath", "channelopathy", "conduction", "pericardi",
    "endocardi", "myocardi", "pulmonary", "tricuspid", "mitral",
    "ebstein", "tetralogy", "truncus", "heterotaxy", "situs",
    "dextrocardia", "hypoplastic", "transposition", "stenosis",
    "coarctation", "patent ductus", "fontan", "norwood",
    "charge", "noonan", "marfan", "loeys-dietz",
    "alagille", "kabuki", "holt-oram", "rasopathy",
)


def _contains_cardiac_root(text: str | None) -> bool:
    """Case-insensitive substring check: does ``text`` contain any
    cardiac root term?"""
    if not text:
        return False
    lower = text.lower()
    return any(root in lower for root in CARDIAC_ROOT_TERMS)


def _normalise_gencc_disease(title: str) -> list[str]:
    """Normalise a single GenCC disease title to ≥0 keywords suitable
    for ClinVar PhenotypeList substring matching.

    Steps (in order):
      a. lowercase
      b. strip parenthetical suffixes — "cardiomyopathy (hypertrophic)"
         → "cardiomyopathy"
      c. strip ", with or without X" tails (and everything after) —
         "testicular anomalies with or without CHD" → "testicular anomalies"
      d. strip end-of-string "type N" / "types N" — "Long QT syndrome
         type 1" → "Long QT syndrome"
      e. strip end-of-string gene-name suffix appended by GenCC —
         "structural heart disease - gata4" → "structural heart disease"
      f. split on " and " and "," into separate keywords
      g. drop keywords shorter than 4 chars after stripping

    Returns deduplicated list of normalised keyword parts (lowercased,
    case-insensitive dedup). Empty list when the input is blank or
    every part is too short.
    """
    if not title:
        return []
    t = title.lower().strip()
    t = _PAREN_RE.sub("", t)
    t = _WITH_OR_WITHOUT_RE.sub("", t)
    t = _TYPE_SUFFIX_RE.sub("", t)
    t = _GENE_SUFFIX_RE.sub("", t)
    parts = [p.strip() for p in _COMPOUND_SPLIT_RE.split(t) if p.strip()]
    seen: set[str] = set()
    out: list[str] = []
    for p in parts:
        if len(p) < _MIN_KEYWORD_LEN:
            continue
        if p in seen:
            continue
        seen.add(p)
        out.append(p)
    return out


def build_cardiac_keywords(
    phenotype_strings: list[str] | None,
    gencc_disease_names: list[str] | None,
    gene_symbol: str,
) -> list[str]:
    """Build the ClinVar landscape's ``condition_keywords`` filter from
    the gene's actual ClinVar phenotype vocabulary + its GenCC disease
    names.

    Data-driven replacement for the prior anticipate-every-string
    pattern-expansion approach. The previous logic curated patterns
    like "atrial septal" → ["atrial septal defect", "ASD"] and missed
    GATA4's dominant ClinVar phenotype "atrioventricular septal defect
    4" because "atrial septal" is not a substring of "atrioventricular
    septal". The new approach queries ClinVar directly for the gene's
    phenotype strings, so whatever vocabulary ClinVar uses is in the
    keyword list by construction.

    Inputs:
      - ``phenotype_strings`` — lowercased, individual disease titles
        already split out of ClinVar PhenotypeList for this gene by
        ``get_gene_phenotype_strings``.
      - ``gencc_disease_names`` — raw disease titles from the gene's
        GenCC submissions; passed through ``_normalise_gencc_disease``
        before filtering.
      - ``gene_symbol`` — used to construct the gene-eponymous
        fallbacks ("GENE-related disorder", "GENE syndrome") that
        ClinVar uses for vague entries.

    Logic:
      1. Keep ClinVar phenotype strings containing a CARDIAC_ROOT_TERM.
      2. Keep normalised GenCC names containing a CARDIAC_ROOT_TERM.
      3. Always add a small gene-eponymous + "cardiovascular phenotype"
         fallback set so vague ClinVar entries match even when nothing
         else does.

    Returns a deduplicated lowercased list. Empty input yields just
    the fallback set — never an empty list.

    Spec note: the user's signature was
    ``build_cardiac_keywords(phenotype_strings, gencc_disease_names)``
    but the gene-eponymous fallback needs the gene symbol, so a 3rd
    positional argument was added.
    """
    out: list[str] = []
    seen: set[str] = set()

    def _add(keyword: str | None) -> None:
        if not keyword:
            return
        kl = keyword.strip().lower()
        if not kl or kl in seen:
            return
        seen.add(kl)
        out.append(kl)

    for p in phenotype_strings or []:
        if _contains_cardiac_root(p):
            _add(p)
    for raw in gencc_disease_names or []:
        for part in _normalise_gencc_disease(raw):
            if _contains_cardiac_root(part):
                _add(part)
    gene_lower = (gene_symbol or "").strip().lower()
    if gene_lower:
        _add(f"{gene_lower}-related disorder")
        _add(f"{gene_lower} syndrome")
    _add("cardiovascular phenotype")
    return out


def _gene_disease_categories(evidence: dict) -> set[str]:
    """The cardiovascular disease categories this gene is actually curated for.

    Union of the category slugs already attached to the gene's GenCC
    submissions, its ClinGen curations and its PanelApp panels, plus the
    congenital-heart bucket when CHDgene lists the gene.

    CHDgene publishes no per-disease categorisation, so a listing can only
    say "this gene causes congenital heart disease" — enough to let a CHD
    phenotype count as matched, and deliberately nothing wider.

    Used by the validity card's unmatched-phenotype line. An empty set means
    no source has a curated disease for this gene, in which case the caller
    stays silent rather than reporting every term as unmatched.
    """
    cats: set[str] = set()
    gencc = evidence.get("gencc") or {}
    for key in ("submissions", "clingen_submissions"):
        for sub in (gencc.get(key) or []):
            if isinstance(sub, dict) and sub.get("category"):
                cats.add(sub["category"])
    panelapp = evidence.get("panelapp") or {}
    for panel in (panelapp.get("panels_found") or []):
        if isinstance(panel, dict) and panel.get("category"):
            cats.add(panel["category"])
    chdgene = evidence.get("chdgene") or {}
    if chdgene.get("ok") and chdgene.get("listed"):
        cats.add("Congenital heart disease")
    return cats


def _protein_position_from_vep(vep: dict) -> int | None:
    """Extract the integer residue position from VEP's ``hgvsp`` for PM5
    lookup. Returns None when VEP failed, the variant isn't missense, or
    the protein consequence carries no parseable position (synonymous,
    splice, in-frame indel, etc.). Missense_variant is the only VEP
    consequence PM5 applies to under the 2015 Richards framework."""
    if not isinstance(vep, dict) or not vep.get("ok"):
        return None
    consequence = (vep.get("most_severe_consequence") or "").lower()
    if "missense_variant" not in consequence:
        return None
    return _hgvsp_int_position(vep.get("hgvsp") or "")


def _mane_protein_position_from_vep(vep: dict) -> int | None:
    """MANE-transcript counterpart of ``_protein_position_from_vep`` — the
    residue number the ClinVar PM5/PS1 retrieval must MATCH on, because
    ClinVar's ``Name`` field is written on its preferred (MANE Select)
    transcript while ``vep["hgvsp"]`` is on the transcript the curator
    supplied. See ``_mane_hgvsp_from_vep`` for the full rationale and the
    GATA4 worked example.

    Same missense gate as ``_protein_position_from_vep``. Returns None when
    no MANE annotation is available, in which case the caller keeps matching
    on the supplied-transcript number (unchanged, pre-fix behaviour) rather
    than guessing an offset."""
    if _protein_position_from_vep(vep) is None:
        return None
    return _hgvsp_int_position(_mane_hgvsp_from_vep(vep) or "")


_PROBAND_AA_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")


def _proband_aa_from_vep(vep: dict) -> tuple[str | None, str | None]:
    """Extract the proband's ``(original_aa, alt_aa)`` 3-letter tokens from
    VEP's ``hgvsp`` for the same-amino-acid PS1/PM5 self-exclusion guard.

    Returns ``(None, None)`` when VEP failed, the variant is not a simple
    missense (no ``p.XxxNNNYyy`` token — synonymous, splice, in-frame indel,
    stop-gain etc.), or the token can't be parsed. Both criteria need the
    alt AA: PM5 drops any candidate encoding the SAME substitution as the
    proband (its own ClinVar record); PS1 keeps ONLY same-substitution
    comparison variants (ACMG-2015 "same amino acid change").
    """
    if not isinstance(vep, dict) or not vep.get("ok"):
        return None, None
    consequence = (vep.get("most_severe_consequence") or "").lower()
    if "missense_variant" not in consequence:
        return None, None
    hgvsp = vep.get("hgvsp") or ""
    m = _PROBAND_AA_RE.search(hgvsp.split(":")[-1])
    if not m:
        return None, None
    return m.group(1), m.group(3)


def _smallest_containing_domain(
    domains: list[dict] | None, aa_pos: int | None,
) -> dict | None:
    """Pick the most-specific UniProt domain containing ``aa_pos``.

    When multiple domains overlap the residue we want the tightest
    annotation (e.g. the SH3 subdomain inside a larger kinase domain).
    Strict containment only — the user spec is explicit that near-boundary
    residues (within ±5 aa) should NOT be treated as in-domain because
    PM1 calls for the variant to sit inside the established hotspot, not
    next to it.
    """
    if not isinstance(aa_pos, int) or aa_pos <= 0 or not domains:
        return None
    hits = [
        d for d in domains
        if isinstance(d.get("start"), int)
        and isinstance(d.get("end"), int)
        and d["start"] <= aa_pos <= d["end"]
    ]
    if not hits:
        return None
    hits.sort(key=lambda d: (d["end"] - d["start"], d.get("name") or ""))
    return hits[0]


_UNSTRUCTURED_REGION_HINTS = (
    "disordered", "low complexity", "low-complexity",
    "compositionally biased", "polar residues", "basic residues",
    "acidic residues",
)


def _is_unstructured_region(feat: dict) -> bool:
    if (feat.get("type") or "").lower() != "region":
        return False
    desc = (feat.get("description") or "").lower()
    return any(h in desc for h in _UNSTRUCTURED_REGION_HINTS)


def _describe_outside_domain_context(
    features: list[dict] | None,
    domains: list[dict] | None,
    aa_pos: int,
) -> tuple[str, str]:
    """Build a richer (reason_short, sentence) when the variant is at a
    known residue but not inside any UniProt Domain feature.

    The goal is for the PM1 box to read coherently against the Protein
    tab's "Domain context" row (which lists every overlapping feature,
    not just Domain-type). The PM1 verdict box would otherwise just say
    "not in any annotated UniProt domain", leaving the curator to
    reconcile it against a visible "Disordered" label on the same tab.

    Returns ("reason_short", "user-facing sentence") so the caller can
    log the short form and surface the sentence.
    """
    feats = features or []
    doms = domains or []
    overlap_unstructured = next(
        (f for f in feats if _is_unstructured_region(f)
         and f.get("start") is not None and f.get("end") is not None
         and f["start"] <= aa_pos <= f["end"]),
        None,
    )
    before = sorted(
        (d for d in doms if isinstance(d.get("end"), int) and d["end"] < aa_pos),
        key=lambda d: aa_pos - d["end"],
    )
    after = sorted(
        (d for d in doms if isinstance(d.get("start"), int) and d["start"] > aa_pos),
        key=lambda d: d["start"] - aa_pos,
    )
    nearest_before = before[0] if before else None
    nearest_after = after[0] if after else None

    def _dom_label(d: dict) -> str:
        return f"{d.get('name') or '(unnamed)'} (aa {d['start']}-{d['end']})"

    if overlap_unstructured:
        span = f"aa {overlap_unstructured['start']}-{overlap_unstructured['end']}"
        if nearest_before and nearest_after:
            sentence = (
                f"Variant lies in a disordered region ({span}) between "
                f"{_dom_label(nearest_before)} and {_dom_label(nearest_after)} — "
                "PM1 does not apply to unstructured regions."
            )
        elif nearest_after:
            sentence = (
                f"Variant lies in the disordered N-terminal region ({span}), "
                f"before {_dom_label(nearest_after)} — PM1 does not apply to "
                "unstructured regions."
            )
        elif nearest_before:
            sentence = (
                f"Variant lies in the disordered C-terminal region ({span}), "
                f"after {_dom_label(nearest_before)} — PM1 does not apply to "
                "unstructured regions."
            )
        else:
            sentence = (
                f"Variant lies in a disordered region ({span}) — PM1 does "
                "not apply to unstructured regions."
            )
        return "in_disordered_region", sentence

    if nearest_before and nearest_after:
        sentence = (
            f"Variant lies in the linker between {_dom_label(nearest_before)} "
            f"and {_dom_label(nearest_after)} — not within any annotated "
            "UniProt domain, PM1 evidence not applicable."
        )
        return "in_inter_domain_linker", sentence
    if nearest_after:
        sentence = (
            f"Variant lies N-terminal of the first annotated domain "
            f"({_dom_label(nearest_after)}) — PM1 evidence not applicable."
        )
        return "n_terminal_of_first_domain", sentence
    if nearest_before:
        sentence = (
            f"Variant lies C-terminal of the last annotated domain "
            f"({_dom_label(nearest_before)}) — PM1 evidence not applicable."
        )
        return "c_terminal_of_last_domain", sentence
    if not doms:
        return (
            "no_domains_annotated",
            "UniProt has no Domain-type features annotated for this gene — "
            "PM1 cannot be assessed from domain evidence.",
        )
    return (
        "outside_all_domains",
        "Variant not within an annotated UniProt domain — "
        "PM1 evidence not applicable.",
    )


async def _build_domain_plp_evidence(
    gene: str, vep_data: dict, uniprot_data: dict,
) -> dict:
    """``domain_plp`` payload, plus the PVS1 initiation-codon sub-payload.

    Thin wrapper over ``_build_domain_plp_core``. The start-loss lookup rides
    on this task rather than becoming its own SSE source ON PURPOSE: the
    frontend's DB_SOURCES list drives both the loading rows and the
    `dbCompletion` gate, and `setDbRowState` on a source with no row is the
    same class of failure as the D0 blank page. Nesting it here changes no
    SSE contract at all.

    It also cannot live inside the core: that function early-returns
    `not_applicable` when the residue is outside any annotated domain, which
    is exactly what happens for a start-loss variant at residue 1.
    """
    payload = await _build_domain_plp_core(gene, vep_data, uniprot_data)
    most_severe = (
        (vep_data.get("most_severe_consequence") or "")
        if isinstance(vep_data, dict) else ""
    )
    if "start_lost" in most_severe.lower():
        payload["pvs1_start_loss"] = await _build_pvs1_start_loss_evidence(
            gene, vep_data, uniprot_data,
        )
    return payload


async def _build_pvs1_start_loss_evidence(
    gene: str, vep_data: dict, uniprot_data: dict,
) -> dict:
    """Inputs for the ClinGen SVI / Abou Tayoun 2018 initiation-codon arm.

    The Moderate branch needs one fact: is there a reported pathogenic variant
    5' of the next downstream in-frame Methionine? So this resolves that Met
    from the UniProt residue string and counts P/LP ClinVar records in
    ``[1, met - 1]``, reusing ``get_domain_plp_evidence`` — which already does
    residue-window P/LP counting in ClinVar's MANE numbering frame.

    Every failure path returns ``next_inframe_met: None`` with a reason, which
    ``_pvs1_start_loss_strength`` reads as "Moderate not demonstrated" and
    resolves to PVS1_Supporting — the paper's default. Nothing here can push
    the strength UP by failing.
    """
    seq = (uniprot_data or {}).get("sequence") if isinstance(uniprot_data, dict) else None
    met = _next_inframe_met(seq)
    if met is None:
        return {
            "ok": True,
            "next_inframe_met": None,
            "upstream_plp": 0,
            "reason": (
                "no downstream in-frame Met in the UniProt sequence"
                if seq else
                "UniProt sequence unavailable (cache predates the `sequence` field)"
            ),
        }
    if met < 2:
        return {"ok": True, "next_inframe_met": None, "upstream_plp": 0,
                "reason": "next in-frame Met is the initiator itself"}
    res = await get_domain_plp_evidence(
        gene, aa_start=1, aa_end=met - 1, numbering_offset=0,
    )
    if not res.get("ok") or res.get("not_applicable"):
        return {"ok": True, "next_inframe_met": met, "upstream_plp": 0,
                "reason": res.get("reason") or res.get("error")
                or "ClinVar window query unavailable"}
    return {
        "ok": True,
        "next_inframe_met": met,
        "upstream_plp": int(res.get("total_plp") or 0),
        "count_P": res.get("count_P"),
        "count_LP": res.get("count_LP"),
        "has_two_star_plp": res.get("has_two_star_plp"),
        "top_variants": res.get("top_variants") or [],
        "numbering_offset": 0,
    }


async def _build_domain_plp_core(
    gene: str, vep_data: dict, uniprot_data: dict,
) -> dict:
    """Compose the ``domain_plp`` evidence payload: pick the smallest
    UniProt domain containing the variant residue and query ClinVar for
    P/LP variants in that range. Bundles the PM1 verdict so the
    downstream prompt + frontend don't have to repeat the threshold
    logic.
    """
    most_severe = (vep_data.get("most_severe_consequence") or "") if isinstance(vep_data, dict) else ""
    if not isinstance(uniprot_data, dict) or not uniprot_data.get("ok"):
        payload = {
            "ok": True,
            "not_applicable": True,
            "reason": "UniProt lookup unavailable",
        }
        status, sentence = _pm1_assessment_from_domain_plp(payload, most_severe)
        payload["pm1_status"] = status
        payload["pm1_assessment"] = sentence
        return payload
    aa_pos = _any_protein_position_from_vep(vep_data)
    if aa_pos is None:
        payload = {
            "ok": True,
            "not_applicable": True,
            "reason": (
                "no protein position (variant is non-coding or VEP "
                "returned no hgvsp)"
            ),
            "variant_position": None,
            "in_domain": False,
        }
        status, sentence = _pm1_assessment_from_domain_plp(payload, most_severe)
        payload["pm1_status"] = status
        payload["pm1_assessment"] = sentence
        return payload
    mane_aa_pos = _any_mane_protein_position_from_vep(vep_data)
    _up_len = uniprot_data.get("length")
    if (
        mane_aa_pos is not None
        and isinstance(_up_len, int) and _up_len > 0
        and mane_aa_pos > _up_len >= aa_pos
    ):
        mane_aa_pos = None
    match_pos = mane_aa_pos if mane_aa_pos is not None else aa_pos
    domain = _smallest_containing_domain(uniprot_data.get("domains"), match_pos)
    if domain is None:
        short, sentence = _describe_outside_domain_context(
            uniprot_data.get("features"),
            uniprot_data.get("domains"),
            match_pos,
        )
        payload = {
            "ok": True,
            "not_applicable": True,
            "reason": short,
            "pm1_assessment": sentence,
            "variant_position": aa_pos,
            "matched_protein_position": match_pos,
            "numbering_differs": match_pos != aa_pos,
            "in_domain": False,
        }
        status, sentence = _pm1_assessment_from_domain_plp(payload, most_severe)
        payload["pm1_status"] = status
        payload["pm1_assessment"] = sentence
        return payload
    base = await get_domain_plp_evidence(
        gene, domain["start"], domain["end"], exclude_position=match_pos,
        numbering_offset=0,
    )
    if not base.get("ok"):
        return base
    base["domain_name"] = domain.get("name") or "(unnamed)"
    base["variant_position"] = aa_pos
    base["matched_protein_position"] = match_pos
    base["numbering_differs"] = match_pos != aa_pos
    base["in_domain"] = True
    status, sentence = _pm1_assessment_from_domain_plp(base, most_severe)
    base["pm1_status"] = status
    base["pm1_assessment"] = sentence
    return base


async def build_same_site_evidence(
    evidence: dict, gene: str, exon_data: dict | None = None,
) -> dict:
    """Other variants reported at the SAME SITE as the proband's.

    Answers the curator question this was built for — "in the gnomAD section,
    is it possible to include if different variants at the same site were
    reported?" — with the three senses of "same site" kept apart, because
    they carry different evidentiary weight:

      * ``same_nucleotide`` — other gnomAD alleles at the variant's own base.
        Population-frequency context; also the honest place to note that the
        alleles at one base can have different rsIDs (chr8:11750234 G>T has
        none while G>A there is rs1205549216), which is why nothing here is
        keyed on rsID.
      * ``same_codon``      — the other two bases of the codon. Still
        frequency context, but a change here can alter the same residue.
      * ``clinvar_same_residue`` — ClinVar records at the same RESIDUE. This
        is the ACMG-relevant one (PM5 is defined at the residue), and every
        record carries an explicit ``pm5_eligible`` flag plus a reason,
        so a synonymous or VUS record can never read as PM5 support.

    Everything is served from the local DBs; nothing new is fetched. Returns
    ``available: False`` with a reason rather than raising, so a missing codon
    span or an off-panel gene degrades to a quiet absence in the UI.
    """
    vep = evidence.get("vep") or {}
    out: dict = {"ok": True, "available": False}
    if not vep.get("ok"):
        out["reason"] = "VEP unavailable"
        return out
    exons = None
    codon_source = None
    codon_note = None
    if isinstance(exon_data, dict) and exon_data.get("ok"):
        exons = exon_data.get("exons")
        codon_source = "picked transcript"
        codon_note = None
    else:
        mane_enst = (vep.get("mane_select_accession") or "")
        if not mane_enst:
            for row in (vep.get("transcript_consequences_all") or []):
                if isinstance(row, dict) and row.get("is_mane_select"):
                    mane_enst = row.get("mane_select_accession") or ""
                    break
        mane_enst = mane_enst.split(".")[0]
        if not mane_enst.upper().startswith("ENST"):
            mane_enst = ""
        if mane_enst:
            ladder = await fetch_transcript_exons(mane_enst)
            if isinstance(ladder, dict) and ladder.get("ok"):
                exons = ladder.get("exons")
                codon_source = f"MANE {mane_enst}"
                codon_note = (
                    f"Ensembl serves an exon model only for an unversioned ENST "
                    f"accession, so the codon span was proven contiguous against "
                    f"the MANE transcript {mane_enst} rather than against "
                    f"{vep.get('transcript_id') or 'the picked transcript'}."
                )
    codon = codon_genomic_positions(vep, exons)
    chrom = vep.get("seq_region_name")
    pos = vep.get("start")
    positions = codon or ([pos] if isinstance(pos, int) else [])
    if not positions:
        out["reason"] = "no genomic position"
        return out
    freq = await fetch_same_site_frequencies(chrom, positions)
    disp_pos = _protein_position_from_vep(vep)
    mane_pos = _mane_protein_position_from_vep(vep)
    _, proband_alt_aa = _proband_aa_from_vep(vep)
    proband_c = (vep.get("hgvsc") or "").split(":")[-1] or None
    same_residue = await get_same_residue_records(
        gene, disp_pos, mane_pos,
        proband_alt_aa=proband_alt_aa, proband_hgvs_c=proband_c,
        proband_hgvs_c_mane=_mane_hgvsc_from_vep(vep),
    )
    own_id = gnomad_variant_id_with_provenance(vep)[0]
    alleles = freq.get("alleles") or []
    same_nt = [
        a for a in alleles
        if a["position"] == pos and a["variant_id"] != own_id
    ]
    same_codon = [a for a in alleles if a["position"] != pos]
    return {
        "ok": True,
        "available": bool(freq.get("available")) or bool(same_residue.get("records")),
        "codon_positions": codon,
        "codon_span_available": codon is not None,
        "codon_span_source": codon_source if codon else None,
        "codon_span_note": codon_note if codon else None,
        "codon_span_confirmed_on_picked_transcript": codon_source == "picked transcript",
        "variant_position": pos,
        "chrom": chrom,
        "gnomad_available": bool(freq.get("available")),
        "gnomad_reason": freq.get("reason"),
        "gnomad_note": freq.get("note"),
        "same_nucleotide": same_nt,
        "same_codon": same_codon,
        "clinvar_same_residue": same_residue,
        "display_residue": disp_pos,
        "matched_residue": mane_pos,
        "numbering_differs": bool(
            mane_pos is not None and disp_pos is not None and mane_pos != disp_pos
        ),
    }


async def _vep_for_input(
    parsed: dict,
    gene: str,
    user_transcript: str | None,
    build: str,
) -> dict:
    """Run the VEP call appropriate for the parsed input.

    HGVS input goes through ``fetch_vep`` exactly as before. Coordinate
    input always lands on the GRCh38 VEP endpoint, lifting GRCh37 input
    through Ensembl's ``/map`` service first:

      - The downstream pipeline (gnomAD, SpliceAI, ClinVar local DB,
        AlphaMissense tabix) is GRCh38-only, so unifying annotation on
        GRCh38 keeps every coordinate comparison consistent.
      - The GRCh38 VEP endpoint returns MANE Select / Plus Clinical
        flags on each transcript_consequence (with ``mane=1`` in
        VEP_PARAMS); the GRCh37 endpoint does not, which would leave
        the Summary tab's "⚠ not MANE Select" warning firing on
        every GRCh37 submission.

    The original GRCh37 coordinates are preserved on the result as
    ``input_coords`` so the Summary header can render
    "Input: chr14:23898487:C:T (GRCh37 → lifted to GRCh38)".
    """
    if parsed["format"] == "hgvs":
        return await fetch_vep(gene, parsed["hgvs"], user_transcript)

    chrom = parsed["chrom"]
    pos = parsed["pos"]
    ref = parsed["ref"]
    alt = parsed["alt"]

    if build == "GRCh37":
        lifted = await liftover_grch37_to_grch38(chrom, pos)
        if lifted is None:
            return {
                "ok": False,
                "error": (
                    f"Liftover from GRCh37 to GRCh38 failed for "
                    f"chr{chrom}:{pos}. Downstream services require "
                    f"GRCh38 coordinates — please verify the input or "
                    f"resubmit in GRCh38."
                ),
                "failure_kind": "liftover",
                "input_format": "coordinates",
                "input_build": "GRCh37",
                "input_coords": f"{chrom}-{pos}-{ref}-{alt}",
            }
        lifted_chrom, lifted_pos = lifted
        lifted_coords = f"{lifted_chrom}-{lifted_pos}-{ref}-{alt}"
        vep_data = await fetch_vep_by_coordinates(
            lifted_chrom, lifted_pos, ref, alt, build="GRCh38"
        )
        if not vep_data.get("ok"):
            return vep_data
        vep_data["input_build"] = "GRCh37"
        vep_data["input_coords"] = f"{chrom}-{pos}-{ref}-{alt}"
        vep_data["lifted_coords_grch38"] = lifted_coords
        return vep_data

    return await fetch_vep_by_coordinates(chrom, pos, ref, alt, build="GRCh38")


def _result(task: asyncio.Task, name: str = "?") -> dict:
    """Resolve a finished task to a JSON-safe dict, including failures."""
    exc = task.exception()
    if exc is not None:
        log.error("DB task %s raised: %r", name, exc, exc_info=exc)
        return {"ok": False, "error": repr(exc)}
    value = task.result()
    if isinstance(value, dict):
        return value
    return {"ok": False, "error": f"unexpected return type: {type(value).__name__}"}


def _vep_failure_message(vep: dict) -> tuple[str, str]:
    """Map a failed VEP result to (user-message, kind) for a fatal SSE error,
    distinguishing a transient outage (retry) from an unresolvable variant.
    VEP is the lynchpin for every variant-level source + PVS1, so a failure can
    never yield a meaningful classification — we surface it instead of silently
    degrading to a misleading 'VUS'.

    The raw VEP error (long, JSON-y developer noise) is LOGGED server-side but
    kept OUT of the curator-facing message — they get one short, actionable
    sentence with the most-likely cause first instead."""
    kind = (vep.get("failure_kind") or "not_found") if isinstance(vep, dict) else "not_found"
    err = (vep.get("error") if isinstance(vep, dict) else None) or "unknown error"
    if kind == "transient":
        return (
            "The variant-annotation service (Ensembl VEP) is temporarily "
            "unavailable, so this variant could not be resolved. Please retry "
            "in a moment — no classification was produced.",
            "vep_transient",
        )
    if kind == "liftover":
        return (f"{err} No classification was produced.", "vep_unresolved")
    log.info("VEP unresolved (%s): %s", kind, err)
    if "reference allele" in err.lower():
        return (
            "The gene and variant don’t match: at this position the reference "
            "base isn’t the one in your variant. The most likely cause is the "
            "wrong gene symbol; also check the transcript (an unusual isoform "
            "can shift the c. position) and the genome build, then retry. "
            "No classification was produced.",
            "vep_gene_mismatch",
        )
    return (
        "This variant couldn’t be resolved by Ensembl VEP. Check the gene "
        "symbol, transcript, HGVS punctuation, and genome build, then retry. "
        "No classification was produced.",
        "vep_unresolved",
    )


_CURATION_DEADLINE_S = float(os.environ.get("HEARTVAR_CURATION_DEADLINE_S", "90") or "90")


async def _gather_evidence_sse(req: CurationRequest, state: dict):
    """Run the parallel DB-gather phase, yielding SSE strings as each source
    resolves. Populates state['evidence'], state['variant_id'],
    state['task_elapsed'], state['task_true_elapsed'], state['gather_elapsed']
    in-place.

    Two per-source timings are recorded, and they measure different things —
    see backend/task_timing.py. 'task_true_elapsed' is the source's own duration;
    'task_elapsed' is when the orchestrator got round to collecting it.
    """
    timing: dict = state.setdefault("timing", {})
    t_input_parse_start = perf_counter()

    gene = (req.gene or "").strip()
    if gene:
        _canon = canonicalise_gene_symbol(gene)
        if _canon.get("is_alias") and _canon.get("approved"):
            state["gene_canonicalised"] = {"from": gene, "to": _canon["approved"]}
            log.info("Canonicalised gene alias %s → %s", gene, _canon["approved"])
            gene = _canon["approved"]
        elif not _canon.get("recognized"):
            state["gene_unrecognized"] = gene
            log.info("Gene symbol %s not recognised by HGNC alias map", gene)
    raw_hgvs = (req.hgvs_c or "").strip()
    parsed = parse_variant_input(raw_hgvs)
    timing["input_parse"] = perf_counter() - t_input_parse_start
    if parsed["format"] == "protein":
        yield _sse("error", {
            "message": (
                "Protein-level (p.) HGVS can't be mapped to a unique coding "
                "change (codon degeneracy), so it can't be annotated. Please "
                "supply the coding (c.) HGVS or genomic coordinates."
            ),
            "kind": "input_protein",
        })
        return
    if parsed["format"] == "rsid":
        rsid = parsed["rsid"]
        recoded = await recode_rsid_to_coords(rsid)
        if not recoded:
            yield _sse("error", {
                "message": (
                    f"Could not resolve {rsid} to GRCh38 coordinates via Ensembl. "
                    "Check the rsID, or enter the variant as HGVS or coordinates."
                ),
                "kind": "rsid_unresolved",
            })
            return
        state["input_rsid"] = rsid
        parsed = {
            "format": "coordinates", "build": None,
            "chrom": recoded["chrom"], "pos": recoded["pos"],
            "ref": recoded["ref"], "alt": recoded["alt"],
        }
        rsid_build_override = "GRCh38"
    else:
        rsid_build_override = None
    if parsed["format"] == "unknown":
        yield _sse("error", {
            "message": (
                "Could not parse variant input. Expected HGVS notation "
                "(e.g. 'c.1208G>A', 'NM_000257.4:c.1208G>A'), genomic "
                "coordinates (e.g. 'chr7:117548628:C:T'), or a dbSNP rsID "
                "(e.g. 'rs727503113')."
            ),
        })
        return

    state["input_format"] = parsed["format"]

    pre_resolved_vep: dict | None = None
    if parsed["format"] == "coordinates":
        build = rsid_build_override or (req.genome_build or "GRCh38").strip() or "GRCh38"
        if build not in ("GRCh38", "GRCh37"):
            yield _sse("error", {
                "message": f"Unsupported genome_build {build!r} — expected GRCh38 or GRCh37.",
            })
            return
        t_coord_pre_vep = perf_counter()
        pre_resolved_vep = await _vep_for_input(parsed, gene, None, build)
        timing["coord_pre_vep"] = perf_counter() - t_coord_pre_vep
        if not pre_resolved_vep.get("ok"):
            msg, kind = _vep_failure_message(pre_resolved_vep)
            yield _sse("error", {"message": msg, "kind": kind})
            return
        if not gene:
            derived = (pre_resolved_vep.get("derived_gene_symbol") or "").strip()
            if not derived:
                yield _sse("error", {
                    "message": (
                        "VEP returned no gene symbol for these "
                        "coordinates — variant may be intergenic. "
                        "Please supply a gene symbol explicitly."
                    ),
                })
                return
            gene = derived
            state["gene_source"] = "derived"
        else:
            state["gene_source"] = "user"
        hgvs_c = pre_resolved_vep.get("derived_hgvs") or ""
        user_transcript = None
    else:
        if not gene:
            yield _sse("error", {
                "message": "Gene symbol is required when the variant is supplied as HGVS.",
            })
            return
        state["gene_source"] = "user"
        hgvs_c = parsed["hgvs"]
        user_transcript = parsed["transcript"]
        if user_transcript:
            log.info(
                "Stripped transcript prefix %s from HGVS input (raw=%r → c.=%r)",
                user_transcript, raw_hgvs, hgvs_c,
            )

    state["user_transcript"] = user_transcript
    state["hgvs_c"] = hgvs_c
    state["gene"] = gene
    evidence: dict[str, dict] = {}
    state["evidence"] = evidence
    evidence["hpo_resolved"] = resolved_token_map(req.hpo or "")
    task_started: dict[str, float] = {}
    task_elapsed: dict[str, float] = {}
    true_elapsed: dict[str, float] = {}

    def _track(name: str, task: asyncio.Task) -> asyncio.Task:
        return track_task_timing(name, task, task_started, true_elapsed)

    t_chdgene = perf_counter()
    evidence["chdgene"] = chdgene_lookup(gene)
    task_elapsed["chdgene"] = perf_counter() - t_chdgene
    true_elapsed["chdgene"] = task_elapsed["chdgene"]
    timing["chdgene_sync"] = task_elapsed["chdgene"]
    yield _sse("db_pending", {"source": "chdgene"})
    yield _sse("db_done", {"source": "chdgene", "data": evidence["chdgene"]})

    t_gather_start = perf_counter()
    pheno_keywords = list(resolved_token_map(req.hpo or "").values())
    if pre_resolved_vep is not None:
        async def _return_pre_vep():
            return pre_resolved_vep
        vep_task = asyncio.create_task(_return_pre_vep())
    else:
        vep_task = asyncio.create_task(fetch_vep(gene, hgvs_c, user_transcript))
    tasks: dict[str, asyncio.Task] = {
        "vep": vep_task,
        "clinvar": asyncio.create_task(
            _fetch_clinvar_with_vep_coords(vep_task, gene, hgvs_c)
        ),
        "uniprot": asyncio.create_task(fetch_uniprot(gene)),
        "gtex": asyncio.create_task(fetch_gtex(gene)),
        "alphafold": asyncio.create_task(fetch_alphafold(gene)),
        "fetal_heart": asyncio.create_task(fetch_fetal_heart(gene)),
        "panelapp": asyncio.create_task(fetch_panelapp(gene, req.hpo)),
        "pubmed": asyncio.create_task(
            _fetch_literature_with_vep_aa(
                fetch_pubmed, vep_task, gene, hgvs_c, req.amino_acid or None)
        ),
        "pubtator3": asyncio.create_task(
            _fetch_literature_with_vep_aa(
                fetch_pubtator3, vep_task, gene, hgvs_c, req.amino_acid or None,
                pass_rsid=True)
        ),
        "gencc": asyncio.create_task(fetch_gencc(gene, req.hpo)),
        "medgen": asyncio.create_task(fetch_medgen(gene, req.hpo)),
        "mgi": asyncio.create_task(fetch_mgi(gene)),
        "biogrid": asyncio.create_task(fetch_biogrid(gene)),
        "gene_literature": asyncio.create_task(
            get_gene_disease_literature(gene, pheno_keywords)
        ),
    }
    if req.enable_erepo:
        tasks["erepo"] = asyncio.create_task(fetch_erepo(gene, hgvs_c))
    for _name, _task in tasks.items():
        _track(_name, _task)
    sources = [
        "vep", "gnomad", "spliceai", "alphamissense", "clinvar",
        "uniprot", "gtex", "alphafold", "fetal_heart", "panelapp", "pubmed",
        "pubtator3", "pmcoa",
        "gencc", "medgen", "protvar", "mgi", "biogrid",
        "clinvar_gene_landscape", "clinvar_pm5_candidates",
        "domain_plp",
        "opentargets_evidence", "gene_literature",
    ]
    if req.enable_erepo:
        sources.append("erepo")
    for src in sources:
        yield _sse("db_pending", {"source": src})

    variant_id: str | None = None
    assembly: str | None = None
    post_vep_started = False
    _LITERATURE_SOURCES = frozenset({
        "pubmed", "pubtator3", "pmcoa", "gene_literature",
    })
    scoring_ready_sent = False

    pmcoa_started = False
    pubmed_lit_done = False
    pubtator3_lit_done = False
    pubmed_pmids: list[str] = []
    pubtator3_pmids: list[str] = []

    def _maybe_start_pmcoa():
        nonlocal pmcoa_started
        if pmcoa_started or not (pubmed_lit_done and pubtator3_lit_done):
            return
        pmcoa_started = True
        seen: set[str] = set()
        union: list[str] = []
        for pid in pubmed_pmids + pubtator3_pmids:
            if pid and pid not in seen:
                seen.add(pid)
                union.append(pid)
        tasks["pmcoa"] = _track("pmcoa", asyncio.create_task(fetch_pmc_fulltext(union)))

    landscape_keyword_source: str | None = None
    vep_done = False
    uniprot_done = False
    domain_plp_started = False

    def _maybe_start_domain_plp():
        nonlocal domain_plp_started
        if domain_plp_started or not (vep_done and uniprot_done):
            return
        domain_plp_started = True
        tasks["domain_plp"] = _track(
            "domain_plp",
            asyncio.create_task(
                _build_domain_plp_evidence(
                    gene,
                    evidence.get("vep") or {},
                    evidence.get("uniprot") or {},
                )
            ),
        )

    _deadline = perf_counter() + _CURATION_DEADLINE_S
    while tasks:
        remaining = _deadline - perf_counter()
        if remaining <= 0:
            cut = list(tasks.items())
            for src, t in cut:
                if not t.done():
                    t.cancel()
                yield _sse("db_done", {"source": src, "data": {
                    "ok": False, "failure_kind": "deadline",
                    "error": "High demand, try again soon",
                }})
            await asyncio.gather(*[t for _, t in cut], return_exceptions=True)
            tasks.clear()
            log.warning(
                "Curation deadline %.0fs exceeded (%s %s) — %d source(s) cut",
                _CURATION_DEADLINE_S, gene, hgvs_c, len(cut),
            )
            cut_names = {src for src, _ in cut}
            if not (evidence.get("vep") or {}).get("ok"):
                state.pop("evidence", None)
                yield _sse("error", {
                    "message": (
                        "This curation timed out before the variant could be "
                        "annotated (Ensembl VEP). Please retry — no classification "
                        "was produced."
                    ),
                    "kind": "vep_transient",
                })
                return
            critical_cut = cut_names & {"gnomad", "clinvar"}
            if critical_cut:
                state.pop("evidence", None)
                yield _sse("error", {
                    "message": (
                        "This curation timed out before key frequency / clinical "
                        f"evidence resolved ({', '.join(sorted(critical_cut))}). "
                        "Please retry — no classification was produced, to avoid a "
                        "result computed on incomplete evidence."
                    ),
                    "kind": "vep_transient",
                })
                return
            break
        done, _pending = await asyncio.wait(
            tasks.values(), return_when=asyncio.FIRST_COMPLETED, timeout=remaining,
        )
        if not done:
            continue
        for finished in done:
            src = next(k for k, t in tasks.items() if t is finished)
            task_elapsed[src] = perf_counter() - task_started.get(src, t_gather_start)
            data = _result(finished, src)
            if src == "erepo":
                state["erepo"] = data
                yield _sse("db_done", {"source": "erepo", "data": data})
                del tasks["erepo"]
                continue
            if src == "gencc" and landscape_keyword_source is None:
                async def _landscape_with_keywords(gencc_data=data):
                    phenotype_strings = await get_gene_phenotype_strings(gene)
                    gencc_disease_names = [
                        s.get("disease") or ""
                        for s in (gencc_data.get("submissions") or [])
                    ] if (isinstance(gencc_data, dict) and gencc_data.get("ok")
                          and gencc_data.get("found")) else []
                    keywords = build_cardiac_keywords(
                        phenotype_strings, gencc_disease_names, gene,
                    )
                    state["landscape_keyword_source"] = (
                        "gencc" if (phenotype_strings or gencc_disease_names)
                        else "generic")
                    state["landscape_condition_keywords"] = keywords
                    return await get_gene_variant_landscape(
                        gene, condition_keywords=keywords or None,
                    )

                landscape_keyword_source = "pending"
                tasks["clinvar_gene_landscape"] = _track(
                    "clinvar_gene_landscape",
                    asyncio.create_task(_landscape_with_keywords()),
                )
            if src == "pubmed":
                pubmed_lit_done = True
                v_papers = (
                    (data.get("variant_papers") or [])
                    if isinstance(data, dict) else []
                )
                pubmed_pmids = [
                    str(p.get("pmid")) for p in v_papers if p.get("pmid")
                ]
                _maybe_start_pmcoa()
            if src == "pubtator3":
                pubtator3_lit_done = True
                pt_papers = (
                    (data.get("papers") or [])
                    if isinstance(data, dict) else []
                )
                pubtator3_pmids = (
                    [str(p.get("pmid")) for p in pt_papers if p.get("pmid") and p.get("pmcid")]
                    + [str(p.get("pmid")) for p in pt_papers if p.get("pmid") and not p.get("pmcid")]
                )
                _maybe_start_pmcoa()
            if src == "vep" and isinstance(data, dict) and user_transcript:
                data["user_transcript"] = user_transcript
            if src == "vep" and isinstance(data, dict):
                data["gene_source"] = state.get("gene_source", "user")
                if state.get("gene_canonicalised"):
                    data["gene_canonicalised"] = state["gene_canonicalised"]
                if state.get("gene_unrecognized"):
                    data["gene_unrecognized"] = state["gene_unrecognized"]
            if src == "clinvar_gene_landscape" and isinstance(data, dict):
                data["keyword_source"] = (
                    state.get("landscape_keyword_source") or "generic"
                )
            evidence[src] = data
            yield _sse("db_done", {"source": src, "data": data})
            del tasks[src]
            if not scoring_ready_sent:
                outstanding = set(tasks) - _LITERATURE_SOURCES
                if not outstanding:
                    scoring_ready_sent = True
                    yield _sse("scoring_ready", {
                        "waiting_on": sorted(set(tasks) & _LITERATURE_SOURCES),
                    })
            if src == "vep" and not data.get("ok"):
                msg, kind = _vep_failure_message(data)
                for pending in tasks.values():
                    if not pending.done():
                        pending.cancel()
                await asyncio.gather(*tasks.values(), return_exceptions=True)
                state.pop("evidence", None)
                yield _sse("error", {"message": msg, "kind": kind})
                return
            if src == "vep" and not post_vep_started:
                variant_id, variant_id_canonical = gnomad_variant_id_with_provenance(data)
                assembly = data.get("assembly_name") if isinstance(data, dict) else None
                _allele = (data.get("allele_string") or "") if isinstance(data, dict) else ""
                _indel_hgvs = (
                    (data.get("hgvsc") or "") if isinstance(data, dict) else ""
                ) if "-" in _allele else None
                tasks["gnomad"] = _track(
                    "gnomad",
                    asyncio.create_task(
                        fetch_gnomad(
                            variant_id, gene, indel_hgvs=_indel_hgvs or None,
                            variant_id_canonical=variant_id_canonical,
                        )
                    ),
                )
                tasks["spliceai"] = _track(
                    "spliceai", asyncio.create_task(fetch_spliceai(variant_id, assembly)))
                tasks["protvar"] = _track(
                    "protvar", asyncio.create_task(fetch_protvar(data)))
                tasks["alphamissense"] = _track(
                    "alphamissense", asyncio.create_task(fetch_alphamissense(data)))
                protein_position = _protein_position_from_vep(data)
                mane_protein_position = _mane_protein_position_from_vep(data)
                mane_hgvs_c = _mane_hgvsc_from_vep(data)
                _proband_orig_aa, _proband_alt_aa = _proband_aa_from_vep(data)
                tasks["clinvar_pm5_candidates"] = _track(
                    "clinvar_pm5_candidates",
                    asyncio.create_task(
                        get_pm5_evidence(
                            gene,
                            protein_position,
                            proband_alt_aa=_proband_alt_aa,
                            proband_hgvs_c=hgvs_c,
                            mane_protein_position=mane_protein_position,
                            proband_hgvs_c_mane=mane_hgvs_c,
                        )
                    ),
                )
                ensembl_gene_id = (
                    data.get("gene_id") if isinstance(data, dict) else None
                )
                tasks["opentargets_evidence"] = _track(
                    "opentargets_evidence",
                    asyncio.create_task(
                        fetch_opentargets(ensembl_gene_id, req.hpo, gene_symbol=gene)
                    ),
                )
                post_vep_started = True
                vep_done = True
                _maybe_start_domain_plp()
            elif src == "uniprot":
                uniprot_done = True
                _maybe_start_domain_plp()

    evidence["phenotype_audit"] = audit_phenotype_terms(
        req.hpo or "", _gene_disease_categories(evidence),
    )

    t_uniprot_proc = perf_counter()
    up_for_proc = evidence.get("uniprot") or {}
    vep_for_proc = evidence.get("vep") or {}
    if up_for_proc.get("ok") and up_for_proc.get("found") and vep_for_proc.get("ok"):
        protein_position_for_uniprot = (
            _mane_protein_position_from_vep(vep_for_proc)
            or _protein_position_from_vep(vep_for_proc)
        )
        _proband_alt_aa_uniprot = _proband_aa_from_vep(vep_for_proc)[1]
        evidence["uniprot_variants_processed"] = {
            "ok": True,
            **process_natural_variants(
                up_for_proc.get("natural_variants") or [],
                protein_position_for_uniprot,
                up_for_proc.get("domains") or [],
                window=15,
                proband_variant_aa=_proband_alt_aa_uniprot,
            ),
        }
    timing["uniprot_process"] = perf_counter() - t_uniprot_proc

    state["variant_id"] = variant_id
    state["task_elapsed"] = task_elapsed
    state["task_true_elapsed"] = dict(true_elapsed)
    state["gather_elapsed"] = perf_counter() - t_gather_start
    timing["gather_wall"] = state["gather_elapsed"]
    timing["gather_per_source"] = dict(task_elapsed)
    timing["gather_per_source_true"] = dict(true_elapsed)

    def _fmt_per_source(d: dict[str, float]) -> str:
        return ", ".join(
            f"{k}={v:.2f}s" for k, v in sorted(d.items(), key=lambda x: -x[1])
        )

    log.info(
        "[timing] %s %s — DB gather %.2fs (true per-source: %s)",
        gene, hgvs_c, state["gather_elapsed"], _fmt_per_source(true_elapsed),
    )
    log.info(
        "[timing] %s %s — DB gather %.2fs (harvest per-source: %s)",
        gene, hgvs_c, state["gather_elapsed"], _fmt_per_source(task_elapsed),
    )

    spans = []
    for src, dur in true_elapsed.items():
        started_at = task_started.get(src)
        if started_at is None:
            continue
        start_off = started_at - t_gather_start
        spans.append((start_off + dur, start_off, dur, src))
    spans.sort(reverse=True)
    state["critical_path"] = [
        {"source": src, "start_s": round(st, 3), "end_s": round(en, 3),
         "duration_s": round(du, 3)}
        for en, st, du, src in spans
    ]
    timing["critical_path"] = state["critical_path"]
    if spans:
        log.info(
            "[timing] %s %s — CRITICAL PATH of the %.2fs gather (start->end, "
            "slowest-finishing first): %s",
            gene, hgvs_c, state["gather_elapsed"],
            ", ".join(f"{src} {st:.2f}->{en:.2f}s({du:.2f}s)"
                      for en, st, du, src in spans[:8]),
        )
