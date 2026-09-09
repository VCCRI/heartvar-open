from __future__ import annotations

import asyncio
import json
import logging
import os
import re
from functools import lru_cache
from pathlib import Path

import httpx

from ._http_retry import request_with_retry
from ._offline import offline_strict
from .hpo_client import fetch_hpo_descendants

log = logging.getLogger(__name__)

CARDIAC_PARENT_TERMS: frozenset[str] = frozenset({
    "HP:0001627",
    "HP:0030680",
    "HP:0001638",
    "HP:0011675",
    "HP:0001671",
})

_PARENT_TO_CATEGORIES: dict[str, tuple[str, ...]] = {
    "HP:0001627": ("Congenital heart disease", "Other cardiac"),
    "HP:0030680": ("Aortic / vascular", "Congenital heart disease", "Other cardiac"),
    "HP:0001638": ("hcm", "dcm", "cardiomyopathy_other", "Cardiomyopathy"),
    "HP:0011675": ("channelopathy", "conduction", "af", "Arrhythmia / channelopathy"),
    "HP:0001671": ("Congenital heart disease",),
    "HP:0031653": ("Congenital heart disease", "Other cardiac"),
    "HP:0031546": ("conduction", "Arrhythmia / channelopathy"),
    "HP:0031547": ("channelopathy", "Arrhythmia / channelopathy"),
    "HP:0004890": ("Pulmonary hypertension",),
}

CARDIAC_HPO_DESCENDANTS: dict[str, set[str]] = {}
_DESCENDANT_CACHE_LOCK = asyncio.Lock()

_DESCENDANTS_JSON_PATH = (
    Path(__file__).resolve().parents[1] / "data" / "cardiac_hpo_descendants.json"
)

PANELAPP_API = "https://panelapp-aus.org/api/v1"

# data carries no explicit redistribution licence (see the localization audit).
_PANELAPP_SNAPSHOT_DEFAULT = (
    Path(__file__).resolve().parents[2] / "data" / "panelapp_aus_snapshot.json"
)
_panelapp_snapshot: dict | None = None
_panelapp_snapshot_loaded = False

CONFIDENCE_LABELS = {"3": "green", "2": "amber", "1": "red"}

CARDIOVASCULAR_KEYWORDS: tuple[str, ...] = (
    "congenital heart",
    "cardiomyopath",
    "arrhythmia",
    "aort",
    "channelopath",
    "cardiac",
    "vascular",
    "hypertrophic",
    "dilated cardiomyopath",
    "long qt",
    "brugada",
    "arvc",
    "marfan",
    "connective tissue",
    "chd",
    "septal",
    "avsd",
    "tetralogy",
    "hcm",
    "dcm",
    "restrictive",
    "cpvt",
    "sick sinus",
    "aneurysm",
    "atrial fibrillation",
    "conduction",
    "heart block",
    "coronary",
    "non-compaction",
    "noncompaction",
    "non compaction",
    "hypercholesterol",
    "familial hypercholesterol",
    "hyperlipid",
    "dyslipid",
    "pulmonary hypertension",
    "pulmonary arterial",
    "catecholaminergic",
    "arrhythmogenic",
    "ventricular dysplasia",
    "jervell",
    "lange-nielsen",
    "lange nielsen",
    "short qt",
    "timothy syndrome",
    "naxos",
    "carvajal",
    "fabry",
    "danon",
    "emery-dreifuss",
    "emery dreifuss",
)

CATEGORY_KEYWORDS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("hcm",
     ("hypertrophic cardiomyopath", "hypertrophic", " hcm", "hcm ", "(hcm)")),
    ("dcm",
     ("dilated cardiomyopath", " dcm", "dcm ", "(dcm)")),
    ("cardiomyopathy_other",
     ("arvc", "arrhythmogenic right ventricular", "arrhythmogenic",
      "ventricular dysplasia", "naxos", "carvajal", "danon",
      "restrictive cardiomyopath",
      "non-compaction", "noncompaction", "non compaction",
      "left ventricular non-compaction")),
    ("channelopathy",
     ("long qt", "short qt", "brugada", "channelopath", "cpvt",
      "catecholaminergic polymorphic", "catecholaminergic",
      "jervell", "lange-nielsen", "lange nielsen", "timothy syndrome")),
    ("conduction",
     ("conduction", "heart block", "sick sinus")),
    ("af",
     ("atrial fibrillation",)),
    ("scad",
     ("spontaneous coronary",)),
    ("cad",
     ("coronary artery disease", "coronary disease", "coronary")),
    ("Familial hypercholesterolaemia",
     ("familial hypercholesterol", "hypercholesterol", "hyperlipid", "dyslipid")),
    ("Pulmonary hypertension",
     ("pulmonary hypertension", "pulmonary arterial")),
    ("Cardiomyopathy",
     ("cardiomyopath",)),
    ("Congenital heart disease",
     ("congenital heart", "chd ", " chd", "septal", "avsd", "tetralogy")),
    ("Arrhythmia / channelopathy",
     ("arrhythmia",)),
    ("Aortic / vascular",
     ("aort", "marfan", "connective tissue", "vascular", "aneurysm")),
)

CATEGORY_DISPLAY_LABELS: dict[str, str] = {
    "hcm": "Hypertrophic cardiomyopathy (HCM)",
    "dcm": "Dilated cardiomyopathy (DCM)",
    "cardiomyopathy_other": "Other cardiomyopathy (ARVC / restrictive / non-compaction)",
    "channelopathy": "Channelopathy (long QT / Brugada / CPVT)",
    "conduction": "Conduction disease (heart block / sick sinus)",
    "af": "Atrial fibrillation",
    "cad": "Coronary artery disease",
    "scad": "Spontaneous coronary artery dissection",
}

PHENOTYPE_PHRASE_TO_CATEGORIES: dict[str, tuple[str, ...]] = {
    "hypertrophic cardiomyopathy": ("hcm",),
    "hcm": ("hcm",),
    "dilated cardiomyopathy": ("dcm",),
    "dcm": ("dcm",),
    "arrhythmogenic right ventricular cardiomyopathy": ("cardiomyopathy_other",),
    "arrhythmogenic cardiomyopathy": ("cardiomyopathy_other",),
    "arvc": ("cardiomyopathy_other",),
    "restrictive cardiomyopathy": ("cardiomyopathy_other",),
    "left ventricular non-compaction": ("cardiomyopathy_other",),
    "left ventricular noncompaction": ("cardiomyopathy_other",),
    "non-compaction cardiomyopathy": ("cardiomyopathy_other",),
    "noncompaction cardiomyopathy": ("cardiomyopathy_other",),
    "lvnc": ("cardiomyopathy_other",),
    "cardiomyopathy": ("hcm", "dcm", "cardiomyopathy_other", "Cardiomyopathy"),
    "long qt": ("channelopathy",),
    "long qt syndrome": ("channelopathy",),
    "lqts": ("channelopathy",),
    "brugada": ("channelopathy",),
    "brugada syndrome": ("channelopathy",),
    "cpvt": ("channelopathy",),
    "catecholaminergic polymorphic ventricular tachycardia": ("channelopathy",),
    "channelopathy": ("channelopathy",),
    "atrial fibrillation": ("af", "channelopathy"),
    "af": ("af", "channelopathy"),
    "conduction disease": ("conduction",),
    "heart block": ("conduction",),
    "atrioventricular block": ("conduction",),
    "av block": ("conduction",),
    "sick sinus syndrome": ("conduction",),
    "sinus node dysfunction": ("conduction",),
    "coronary artery disease": ("cad",),
    "cad": ("cad",),
    "spontaneous coronary artery dissection": ("scad",),
    "scad": ("scad",),
    "marfan syndrome": ("Aortic / vascular",),
    "marfan": ("Aortic / vascular",),
    "loeys-dietz syndrome": ("Aortic / vascular",),
    "ehlers-danlos syndrome": ("Aortic / vascular",),
    "aortic aneurysm": ("Aortic / vascular",),
    "aortic dissection": ("Aortic / vascular",),
    "thoracic aortic aneurysm": ("Aortic / vascular",),
    "congenital heart disease": ("Congenital heart disease",),
    "chd": ("Congenital heart disease",),
    "tetralogy of fallot": ("Congenital heart disease",),
    "atrial septal defect": ("Congenital heart disease",),
    "ventricular septal defect": ("Congenital heart disease",),
    "atrioventricular septal defect": ("Congenital heart disease",),
    "atrioventricular canal": ("Congenital heart disease",),
    "endocardial cushion defect": ("Congenital heart disease",),
    "asd": ("Congenital heart disease",),
    "vsd": ("Congenital heart disease",),
    "avsd": ("Congenital heart disease",),
    "coarctation of the aorta": ("Congenital heart disease",),
    "hypoplastic left heart syndrome": ("Congenital heart disease",),
    "hlhs": ("Congenital heart disease",),
    "hypoplastic right heart": ("Congenital heart disease",),
    "right heart hypoplasia": ("Congenital heart disease",),
    "single ventricle": ("Congenital heart disease",),
    "truncus arteriosus": ("Congenital heart disease",),
    "persistent truncus arteriosus": ("Congenital heart disease",),
    "interrupted aortic arch": ("Congenital heart disease",),
    "double outlet right ventricle": ("Congenital heart disease",),
    "dorv": ("Congenital heart disease",),
    "transposition of the great arteries": ("Congenital heart disease",),
    "tga": ("Congenital heart disease",),
    "dextro-transposition": ("Congenital heart disease",),
    "ebstein anomaly": ("Congenital heart disease",),
    "ebstein's anomaly": ("Congenital heart disease",),
    "patent ductus arteriosus": ("Congenital heart disease",),
    "pda": ("Congenital heart disease",),
    "pulmonary stenosis": ("Congenital heart disease",),
    "pulmonary atresia": ("Congenital heart disease",),
    "anomalous pulmonary venous return": ("Congenital heart disease",),
    "tapvr": ("Congenital heart disease",),
    "total anomalous pulmonary venous": ("Congenital heart disease",),
    "partial anomalous pulmonary venous": ("Congenital heart disease",),
    "cor triatriatum": ("Congenital heart disease",),
    "dextrocardia": ("Congenital heart disease",),
    "situs inversus": ("Congenital heart disease",),
    "heterotaxy": ("Congenital heart disease",),
    "arrhythmia": ("channelopathy", "conduction", "af", "Arrhythmia / channelopathy"),
}

_HPO_MAP_CACHE: dict | None = None


def _load_hpo_map() -> dict:
    """Lazy-load the category→HPO-subtree mapping from disk; cached on first call."""
    global _HPO_MAP_CACHE
    if _HPO_MAP_CACHE is None:
        path = Path(__file__).resolve().parents[1] / "data" / "panelapp_hpo_map.json"
        with path.open() as f:
            _HPO_MAP_CACHE = json.load(f)
    return _HPO_MAP_CACHE


def _label(level) -> str:
    if level is None:
        return "unknown"
    return CONFIDENCE_LABELS.get(str(level), str(level))


def _categorize_panel(name: str) -> str | None:
    """Bucket a panel by name into one of the cardiovascular categories.

    Specific sub-types (hcm, dcm, channelopathy, conduction, af, cad,
    scad, cardiomyopathy_other) win over the broad fallback buckets
    ("Cardiomyopathy", "Arrhythmia / channelopathy") — see
    CATEGORY_KEYWORDS for the ordering. Returns None when the panel
    isn't cardiovascular at all."""
    n = (name or "").lower()
    if not any(kw in n for kw in CARDIOVASCULAR_KEYWORDS):
        return None
    for category, kws in CATEGORY_KEYWORDS:
        if any(kw in n for kw in kws):
            return category
    return "Other cardiac"


def _normalise_phrase(s: str) -> str:
    """Aggressive phrase-key normalisation so curator typing variations
    all collapse to the same lookup key:
      • lowercase
      • hyphens treated as spaces (so "non-compaction" matches
        "noncompaction" matches "non compaction")
      • runs of whitespace collapsed to a single space
      • trailing punctuation (`.`, `,`, `;`, `:`, `!`) stripped
    Applied symmetrically to both the dict keys (at module load) and
    each free-text chunk at lookup time."""
    s = s.lower()
    s = s.replace("-", " ")
    s = re.sub(r"\s+", " ", s)
    s = s.strip(" .,;:!?")
    return s


_PHENOTYPE_PHRASE_TO_CATEGORIES_NORMALISED: dict[str, tuple[str, ...]] = {
    _normalise_phrase(k): v for k, v in PHENOTYPE_PHRASE_TO_CATEGORIES.items()
}

CATEGORY_PARENTS: dict[str, tuple[str, ...]] = {
    "hcm": ("Cardiomyopathy",),
    "dcm": ("Cardiomyopathy",),
    "cardiomyopathy_other": ("Cardiomyopathy",),
    "channelopathy": ("Arrhythmia / channelopathy",),
    "conduction": ("Arrhythmia / channelopathy",),
    "af": ("Arrhythmia / channelopathy",),
}


_SUBSTRING_MIN_KEY_LEN = 6

_HPO_ID_RE = re.compile(r"^HP:\d{7}$")

CURATED_PHENOTYPE_ACRONYMS: dict[str, str] = {
    "TOF": "HP:0001636",
    "COA": "HP:0001680",
    "BAV": "HP:0001647",
    "IAA": "HP:0011611",
    "AS": "HP:0001650",
    "PS": "HP:0001642",
    "MS": "HP:0001718",
    "MR": "HP:0001653",
    "AR": "HP:0001659",
    "TR": "HP:0005180",
    "LVNC": "HP:0030682",
    "RCM": "HP:0001723",
    "LVH": "HP:0001712",
    "AVB": "HP:0001678",
    "SSS": "HP:0011704",
    "WPW": "HP:0001716",
    "VT": "HP:0004756",
    "VF": "HP:0001663",
    "AF": "HP:0005110",
    "SCD": "HP:0001645",
    "TAA": "HP:0012727",
    "PAH": "HP:0002092",
}
_CURATED_ACRONYMS_NORMALISED: dict[str, str] = {
    _normalise_phrase(k): v for k, v in CURATED_PHENOTYPE_ACRONYMS.items()
}

_COMPOUND_SPLIT_RE = re.compile(r"\s*[-/+&]\s*")
_QUALIFIER_SPLIT_RE = re.compile(r"\s+with\s+", re.IGNORECASE)
_MIN_COMPOUND_PART_LEN = 2

_SYNONYMS_JSON_PATH = (
    Path(__file__).resolve().parents[1] / "data" / "cardiac_hpo_synonyms.json"
)


@lru_cache(maxsize=1)
def _load_synonym_index() -> dict[str, str]:
    """Lazy-load the normalised-phrase → HP-id index; {} when not yet built."""
    try:
        with _SYNONYMS_JSON_PATH.open("r", encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        log.warning("cardiac HPO synonym index unavailable at %s — free-text "
                    "phenotype recognition falls back to the phrase table",
                    _SYNONYMS_JSON_PATH)
        return {}
    return {k: v for k, v in data.items() if isinstance(v, str)}


def _resolve_term(chunk: str) -> dict | None:
    """Resolve ONE phenotype chunk to a phenotype term, or None.

    Resolution order, most explicit first:
      1. ``HP:NNNNNNN`` — the curator named the term outright
      2. curated acronym map — human-signed-off clinical shorthand
      3. generated synonym index — HPO's own names + EXACT synonyms
      4. hand-written phrase table (exact, then >= 6-char substring)

    Returns ``{"hpo_id", "categories", "source"}``; `hpo_id` is None for a
    phrase-table hit, which carries categories but no ontology term.
    Compound splitting is handled by the caller, so this never recurses.
    """
    upper = chunk.upper()
    if _HPO_ID_RE.match(upper):
        return {"hpo_id": upper, "categories": (), "source": "hpo_id",
                "text": chunk}
    key = _normalise_phrase(chunk)
    if not key:
        return None
    hpo_id = _CURATED_ACRONYMS_NORMALISED.get(key)
    if hpo_id:
        return {"hpo_id": hpo_id, "categories": (), "source": "acronym",
                "text": chunk}
    cats = _PHENOTYPE_PHRASE_TO_CATEGORIES_NORMALISED.get(key)
    if cats:
        return {"hpo_id": None, "categories": cats, "source": "phrase",
                "text": chunk, "display_id": _load_synonym_index().get(key)}
    hpo_id = _load_synonym_index().get(key)
    if hpo_id:
        return {"hpo_id": hpo_id, "categories": (), "source": "synonym",
                "text": chunk}
    cats = _substring_phrase_match(key)
    if cats:
        return {"hpo_id": None, "categories": cats, "source": "phrase_substring",
                "text": chunk}
    return None


def _resolve_chunk_or_parts(chunk: str) -> tuple[list[dict], list[str]]:
    """Resolve a chunk, falling back to splitting it as a compound.

    Returns ``(resolved_terms, unresolved_fragments)``. The whole chunk is
    tried intact first; only if that fails is it split on
    ``_COMPOUND_SPLIT_RE``. This is what makes "DORV-TGA" work: both halves
    were always in the vocabulary, but the joined form matched no key and was
    too short for the substring fallback.
    """
    whole = _resolve_term(chunk)
    if whole is not None:
        return [whole], []
    parts = [p.strip() for p in _COMPOUND_SPLIT_RE.split(chunk) if p.strip()]
    resolved: list[dict] = []
    unresolved: list[str] = []
    if len(parts) >= 2:
        for part in parts:
            if len(part) < _MIN_COMPOUND_PART_LEN:
                unresolved.append(part)
                continue
            got = _resolve_term(part)
            (resolved.append(got) if got else unresolved.append(part))
        if resolved:
            return resolved, unresolved
    head = _QUALIFIER_SPLIT_RE.split(chunk)[0].strip()
    if head and head != chunk:
        got = _resolve_term(head)
        if got:
            return [got], []
    return [], [chunk]


def _substring_phrase_match(normalised_chunk: str) -> tuple[str, ...] | None:
    """Substring fallback for the phrase table. Used when the exact
    ``_normalise_phrase`` lookup misses — handles inputs like
    "hypertrophic cardiomyopathy 1" or "HCM with apical involvement"
    that contain a known phrase but aren't an exact match.

    Returns the categories of the LONGEST matching key (so
    "hypertrophic cardiomyopathy" wins over plain "cardiomyopathy"),
    or None if no key of sufficient length is a substring.
    """
    if not normalised_chunk:
        return None
    best_key: str | None = None
    best_len = 0
    for key in _PHENOTYPE_PHRASE_TO_CATEGORIES_NORMALISED:
        if len(key) < _SUBSTRING_MIN_KEY_LEN:
            continue
        if key in normalised_chunk and len(key) > best_len:
            best_key = key
            best_len = len(key)
    if best_key is None:
        return None
    return _PHENOTYPE_PHRASE_TO_CATEGORIES_NORMALISED[best_key]


def _parse_phenotype_input(
    raw,
) -> tuple[list[str], set[str], list[str], list[str]]:
    """Normalise the phenotype input field to:
      (hpo_ids, phrase_matched_categories, recognised_tokens, unrecognised_tokens)

    The input may mix HPO IDs ("HP:0001644"), clinical acronyms ("DORV-TGA",
    "AS"), and free-text phenotype phrases ("dilated cardiomyopathy"),
    separated by commas or semicolons. Each chunk is resolved by
    :func:`_resolve_chunk_or_parts` — see :func:`_resolve_term` for the
    resolution order and the compound-splitting rule.

    A chunk resolving to an ontology term contributes its HP id to
    `hpo_ids`, from where the descendant closure supplies categories in
    :func:`check_hpo_relevance`; one resolving through the hand-written
    phrase table contributes categories directly to
    `phrase_matched_categories`. Anything unresolved lands in
    `unrecognised_tokens` so callers can surface it to the curator rather
    than silently dropping it. `recognised_tokens` preserves the curator's
    own casing and spacing for display.
    """
    if raw is None:
        return [], set(), [], []
    if isinstance(raw, str):
        chunks = [c.strip() for c in re.split(r"[,;]+", raw) if c.strip()]
    else:
        chunks = [str(c).strip() for c in raw if str(c).strip()]
    hpo_ids: list[str] = []
    phrase_cats: set[str] = set()
    recognised: list[str] = []
    unrecognised: list[str] = []
    for chunk in chunks:
        terms, unresolved = _resolve_chunk_or_parts(chunk)
        unrecognised.extend(unresolved)
        if not terms:
            continue
        for term in terms:
            if term["hpo_id"]:
                hpo_ids.append(term["hpo_id"])
            for c in term["categories"]:
                phrase_cats.add(c)
                phrase_cats.update(CATEGORY_PARENTS.get(c, ()))
        if len(terms) == 1 and terms[0]["source"] == "hpo_id":
            recognised.append(terms[0]["hpo_id"])
        else:
            recognised.append(chunk)
    return hpo_ids, phrase_cats, recognised, unrecognised


def resolve_phenotype_terms(raw) -> list[dict]:
    """Per-chunk report of HOW each submitted phenotype term was interpreted.

    A 2–3 letter acronym is ambiguous by nature — "MS" is mitral stenosis in
    a cardiology tool and multiple sclerosis elsewhere — and the phenotype
    gates PP4, so the reading must be visible to the curator instead of
    silently applied. Drives the interpretation line on the validity card.

    Returns one entry per chunk, in the curator's order:
      ``{token, hpo_id, label, source, expanded, parts}``
    where `expanded` marks a term whose meaning was inferred (acronym or
    synonym) rather than typed outright, and `parts` lists the resolved
    labels of a compound like "DORV-TGA". `hpo_id`/`label` are None for a
    phrase-table hit, which carries categories but no single ontology term.
    """
    from .hpo_labels import resolve_hpo_label

    if raw is None:
        return []
    if isinstance(raw, str):
        chunks = [c.strip() for c in re.split(r"[,;]+", raw) if c.strip()]
    else:
        chunks = [str(c).strip() for c in raw if str(c).strip()]
    def _label_of(term: dict) -> str:
        """Best human label for a resolved term: its ontology name where one
        exists, else the fragment the curator actually typed."""
        for candidate in (term.get("hpo_id"), term.get("display_id")):
            if candidate:
                label = resolve_hpo_label(candidate)
                if label:
                    return label
        return term.get("text") or ""

    out: list[dict] = []
    for chunk in chunks:
        terms, _unresolved = _resolve_chunk_or_parts(chunk)
        parts = [_label_of(t) for t in terms]
        first = terms[0] if terms else None
        hpo_id = first["hpo_id"] if first else None
        out.append({
            "token": chunk,
            "hpo_id": hpo_id,
            "label": (resolve_hpo_label(hpo_id) if hpo_id else None),
            "source": (first["source"] if first else None),
            "expanded": bool(
                first and (first["source"] in ("acronym", "synonym") or len(terms) > 1)
            ),
            "parts": parts,
        })
    return out


def _parse_hpo_terms(hpo_raw) -> list[str]:
    """Legacy alias — returns just the HPO IDs from a mixed input.
    Kept for any external caller; new code should use
    _parse_phenotype_input() so free-text phrases aren't dropped."""
    hpo_ids, _, _, _ = _parse_phenotype_input(hpo_raw)
    return hpo_ids


def _load_local_descendants() -> bool:
    """Populate ``CARDIAC_HPO_DESCENDANTS`` from the precomputed offline
    closure (``backend/data/cardiac_hpo_descendants.json``). Returns True on
    a non-empty load, False when the file is absent/empty/unreadable so the
    caller can fall back to the live JAX fetch."""
    if not _DESCENDANTS_JSON_PATH.is_file():
        return False
    try:
        with _DESCENDANTS_JSON_PATH.open() as f:
            data = json.load(f)
    except (OSError, ValueError) as e:
        log.warning("PanelApp local HPO closure unreadable (%r); using live JAX", e)
        return False
    for hpo_id, cats in data.items():
        if cats:
            CARDIAC_HPO_DESCENDANTS[hpo_id] = set(cats)
    return bool(CARDIAC_HPO_DESCENDANTS)


async def _ensure_descendant_cache() -> None:
    """Populate ``CARDIAC_HPO_DESCENDANTS`` on first use.

    Prefers the precomputed offline closure (no network, no cold-start);
    only when that file is missing does it fall back to fetching the full
    descendant set for each of the 5 ``CARDIAC_PARENT_TERMS`` from JAX, each
    descendant mapped to the union of its parents' categories per
    ``_PARENT_TO_CATEGORIES``.

    Failure-silent: if neither source is available, the cache stays empty
    and ``check_hpo_relevance`` falls back to the curated 61-ID seed map.
    The live fallback costs ~8-10 s on the first call (one /descendants
    request per root); subsequent calls are zero-cost.
    """
    if CARDIAC_HPO_DESCENDANTS:
        return
    async with _DESCENDANT_CACHE_LOCK:
        if CARDIAC_HPO_DESCENDANTS:
            return
        if _load_local_descendants():
            log.info(
                "PanelApp HPO descendant cache loaded from local closure: "
                "%d unique IDs", len(CARDIAC_HPO_DESCENDANTS),
            )
            return
        if offline_strict():
            return
        for parent, cats in _PARENT_TO_CATEGORIES.items():
            descendants = await fetch_hpo_descendants(parent)
            for d in descendants:
                CARDIAC_HPO_DESCENDANTS.setdefault(d, set()).update(cats)
        log.info(
            "PanelApp HPO descendant cache built from live JAX: "
            "%d unique IDs across %d parents",
            len(CARDIAC_HPO_DESCENDANTS), len(_PARENT_TO_CATEGORIES),
        )


def check_hpo_relevance(submitted_input, panel_category: str) -> bool:
    """Whether the proband's submitted phenotype is relevant to a panel's
    disease category.

    Match logic (in order):
      1. Descendant-derived HPO match — submitted HPO ID is in the
         descendant cache built from the 5 cardiac ancestor roots, and
         the panel_category is in the descendant's inherited category
         set. Covers any HPO term within the cardiac subtree of the
         ontology, not just the 61 IDs in the curated seed map.
      2. Curated seed fallback — only consulted when the descendant
         cache is empty (e.g. JAX API was unreachable on this run).
         Uses the 61-ID panelapp_hpo_map.json that was the previous
         hand-curated source.
      3. Parent HPO fallback — any submitted HPO ID in
         CARDIAC_PARENT_TERMS counts as a match against EVERY
         cardiovascular category. Intentional broad-match for curators
         who submit a top-level term.
      4. Free-text phrase match — recognised English phrase (e.g.
         "dilated cardiomyopathy") maps to a category via
         PHENOTYPE_PHRASE_TO_CATEGORIES.

    Returns False only when nothing was parsed, or when none of the
    above matches apply to this specific category.
    """
    hpo, phrase_cats, recognised, _unrecognised = _parse_phenotype_input(submitted_input)
    if not hpo and not phrase_cats:
        log.debug("PanelApp HPO check: no terms recognised for %s", panel_category)
        return False
    if CARDIAC_HPO_DESCENDANTS:
        exact = [t for t in hpo
                 if panel_category in CARDIAC_HPO_DESCENDANTS.get(t, set())]
        match_source = "descendants"
    else:
        relevant = set((_load_hpo_map().get(panel_category) or []))
        exact = [t for t in hpo if t in relevant]
        match_source = "seed_map"
    parent = [t for t in hpo if t in CARDIAC_PARENT_TERMS]
    phrase = panel_category in phrase_cats
    matched = bool(exact) or bool(parent) or phrase
    log.info(
        "PanelApp phenotype check — category=%r submitted=%s "
        "source=%s exact=%s parent=%s phrase=%s matched=%s",
        panel_category, recognised, match_source, exact, parent, phrase, matched,
    )
    return matched


def audit_phenotype_terms(raw, gene_categories) -> dict:
    """Per-term verdict on the curator's phenotype input against the disease
    categories this GENE is actually curated for.

    The per-row ✓ ticks in the validity card are computed from the WHOLE
    submitted input, so a term that matches nothing leaves no trace — it
    looks identical to having submitted no phenotype at all. This walks the
    input term by term and names the ones that matched nothing, so the card
    can say so explicitly.

    `gene_categories` is the union of the categories carried by the gene's
    ClinGen curations, GenCC submissions and PanelApp panels, plus
    'Congenital heart disease' when the gene is a CHDgene entry (CHDgene
    itself carries no per-disease categorisation, so a listing contributes
    the congenital-heart bucket and nothing more).

    Match resolution defers to :func:`check_hpo_relevance` per category —
    the same predicate that drives the ticks — so this can never contradict
    a ✓ shown above it.

    Returns ``{matched, unmatched, unrecognised, suppressed}`` where matched
    / unmatched are ``[{token, label}]`` in the curator's own order, with
    duplicates collapsed. `unrecognised` mirrors
    :func:`_parse_phenotype_input` — tokens that parsed as neither an HPO ID
    nor a known phrase. Those are NOT reported as unmatched: an unparseable
    token is an input problem (already flagged on the PanelApp row), not a
    gene-phenotype mismatch.

    `suppressed` is True when there is nothing to say: no phenotype
    submitted, or the gene has no known categories at all. In the latter
    case every term would look unmatched, which blames the curator's
    phenotype for what is really absent gene coverage — so stay silent.
    """
    empty = {"matched": [], "unmatched": [], "unrecognised": [],
             "interpretations": [], "suppressed": True}
    if raw is None:
        return dict(empty)
    if isinstance(raw, str):
        chunks = [c.strip() for c in re.split(r"[,;]+", raw) if c.strip()]
    else:
        chunks = [str(c).strip() for c in raw if str(c).strip()]
    if not chunks:
        return dict(empty)

    cats = {c for c in (gene_categories or ()) if c}
    _, _, _, unrecognised = _parse_phenotype_input(chunks)
    interpretations = [t for t in resolve_phenotype_terms(chunks) if t["expanded"]]
    if not cats:
        return {"matched": [], "unmatched": [], "unrecognised": unrecognised,
                "interpretations": interpretations, "suppressed": True}

    unrecognised_set = set(unrecognised)
    label_by_token = {}
    for t in resolve_phenotype_terms(chunks):
        label_by_token[t["token"]] = (
            " + ".join(p for p in t["parts"] if p) or t["token"]
        )
    matched: list[dict] = []
    unmatched: list[dict] = []
    seen: set[str] = set()
    for chunk in chunks:
        if chunk in unrecognised_set:
            continue
        key = chunk.upper()
        if key in seen:
            continue
        seen.add(key)
        entry = {"token": chunk, "label": label_by_token.get(chunk) or chunk}
        if any(check_hpo_relevance(chunk, c) for c in cats):
            matched.append(entry)
        else:
            unmatched.append(entry)
    return {"matched": matched, "unmatched": unmatched,
            "unrecognised": unrecognised, "interpretations": interpretations,
            "suppressed": not unmatched}


def _load_panelapp_snapshot() -> dict | None:
    """Lazily load the local PanelApp snapshot. Returns the
    {GENE_UPPER: [record, …]} map, or None when no snapshot is built
    (so the caller falls back to the live API). Cached for the process."""
    global _panelapp_snapshot, _panelapp_snapshot_loaded
    if _panelapp_snapshot_loaded:
        return _panelapp_snapshot
    override = os.environ.get("PANELAPP_SNAPSHOT_PATH", "").strip()
    path = Path(override) if override else _PANELAPP_SNAPSHOT_DEFAULT
    if not path.is_file():
        _panelapp_snapshot = None
        _panelapp_snapshot_loaded = True
        return None
    try:
        with path.open() as f:
            data = json.load(f)
        _panelapp_snapshot = data.get("genes", data) if isinstance(data, dict) else None
    except (OSError, ValueError) as e:
        log.warning("PanelApp snapshot unreadable (%r); using live API", e)
        _panelapp_snapshot = None
    _panelapp_snapshot_loaded = True
    if _panelapp_snapshot is not None:
        log.info("PanelApp snapshot loaded: %d genes from %s",
                 len(_panelapp_snapshot), path)
    return _panelapp_snapshot


async def _panelapp_results(gene: str) -> tuple[list | None, str | None]:
    """Return ``(results, error)`` — the raw gene-on-panel records
    ``fetch_panelapp`` post-processes, sourced from the local snapshot when
    present and the live PanelApp Australia API otherwise.

    Each record matches the live ``/genes/`` item shape:
    ``{"panel": {"id", "name", "version"}, "confidence_level",
    "mode_of_inheritance", "phenotypes", "entity_status"}``. ``error`` is a
    string only when the live call fails (the snapshot path never errors —
    a gene absent from the snapshot yields an empty list)."""
    snapshot = _load_panelapp_snapshot()
    if snapshot is not None:
        return list(snapshot.get(gene.strip().upper(), [])), None
    if offline_strict():
        return [], None
    async with httpx.AsyncClient(
        headers={"User-Agent": "HeartVar/1.0 (mailto:heartvar@victorchang.edu.au)"}
    ) as c:
        r = await request_with_retry(
            c, "GET", f"{PANELAPP_API}/genes/",
            params={"entity_name": gene}, timeout=20.0, name="PanelApp",
        )
        if r is None:
            return None, "PanelApp transient failure after retries"
        if r.status_code != 200:
            return None, f"{r.status_code}: {r.text[:200]}"
        return r.json().get("results", []) or [], None


async def fetch_panelapp(gene: str, hpo: str | None = None) -> dict:
    """Find every cardiovascular PanelApp Australia panel the gene appears on.

    One API call to `/genes/?entity_name={gene}` returns every panel the gene
    is on; we filter to cardiovascular panels client-side via
    `_categorize_panel`, then annotate each with HPO relevance for the
    submitted proband phenotype.

    The result shape exposes:
      - panels_found: per-panel records (id, name, confidence, MOI, category,
        hpo_match, contributes_to_pp4)
      - green/amber/red counts
      - per-category rollup (panel_count, any_green, hpo_match, any_pp4)
      - on_chd_panel: legacy back-compat flag (true if any CHD-category panel)
    """
    (
        submitted_hpo_ids,
        submitted_phrase_cats,
        recognised_inputs,
        unrecognised_inputs,
    ) = _parse_phenotype_input(hpo)
    await _ensure_descendant_cache()
    results, error = await _panelapp_results(gene)
    if error is not None:
        return {"ok": False, "error": error}

    seen_ids: set = set()
    panels_found: list[dict] = []
    for g in results:
        panel = g.get("panel") or {}
        pid = panel.get("id")
        name = panel.get("name", "") or ""
        category = _categorize_panel(name)
        if category is None:
            continue
        if pid in seen_ids:
            continue
        seen_ids.add(pid)
        level = g.get("confidence_level")
        label = _label(level)
        hpo_match = check_hpo_relevance(hpo, category)
        panels_found.append({
            "panel_id": pid,
            "panel_name": name,
            "panel_version": panel.get("version"),
            "panel_url": f"https://panelapp-aus.org/panels/{pid}/",
            "confidence_level": level,
            "confidence": label,
            "moi": g.get("mode_of_inheritance"),
            "category": category,
            "phenotypes": g.get("phenotypes") or [],
            "entity_status": g.get("entity_status"),
            "hpo_match": hpo_match,
            "contributes_to_pp4": (label == "green") and hpo_match,
        })

    order = {"green": 0, "amber": 1, "red": 2, "unknown": 3}
    panels_found.sort(key=lambda p: (
        order.get(p["confidence"], 4),
        p["category"],
        p["panel_name"] or "",
    ))

    green = sum(1 for p in panels_found if p["confidence"] == "green")
    amber = sum(1 for p in panels_found if p["confidence"] == "amber")
    red   = sum(1 for p in panels_found if p["confidence"] == "red")

    categories: dict[str, dict] = {}
    for p in panels_found:
        cat = p["category"]
        slot = categories.setdefault(
            cat, {"panel_count": 0, "any_green": False, "hpo_match": False, "any_pp4": False}
        )
        slot["panel_count"] += 1
        if p["confidence"] == "green":
            slot["any_green"] = True
        if p["hpo_match"]:
            slot["hpo_match"] = True
        if p["contributes_to_pp4"]:
            slot["any_pp4"] = True

    any_green_with_hpo = any(p["contributes_to_pp4"] for p in panels_found)
    any_green_panel = any(p["confidence"] == "green" for p in panels_found)

    return {
        "ok": True,
        "gene": gene,
        "submitted_hpo": recognised_inputs,
        "submitted_hpo_ids": submitted_hpo_ids,
        "submitted_phrase_categories": sorted(submitted_phrase_cats),
        "unrecognised_hpo_tokens": unrecognised_inputs,
        "on_cardiovascular_panel": bool(panels_found),
        "on_chd_panel": any(p["category"] == "Congenital heart disease" for p in panels_found),
        "panels_found": panels_found,
        "green_panels": green,
        "amber_panels": amber,
        "red_panels": red,
        "total_panels": len(panels_found),
        "total_panels_overall": len(results),
        "categories": categories,
        "any_green_with_hpo_match": any_green_with_hpo,
        "any_green_no_hpo_match": any_green_panel and not any_green_with_hpo,
        "best_confidence_label": panels_found[0]["confidence"] if panels_found else None,
        "chd_panel_matches": panels_found,
    }
