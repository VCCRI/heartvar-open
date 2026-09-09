"""PMC Open Access full-text client.

Augments the variant-specific PubMed papers with Results/Methods excerpts
from the PMC Open Access subset, when a free full-text record exists. This
deepens every literature-dependent criterion (PS3 functional assays, PS4
case/cohort counts, PP1 segregation) by giving the model the quantitative
body text — case counts, odds ratios, proband numbers — that abstracts omit.

Pipeline:
  PMID --(NCBI ID converter)--> PMCID --(PMC OA API)--> is it open access?
       --(BioC JSON)--> Results + Methods section text.

Variant-specific papers ONLY (the caller passes at most 3 PMIDs); gene-level
papers are excluded — full text is expensive and only the variant-naming
papers carry the case-level evidence we want.

Non-fatal: every failure mode (no PMCID, not in OA, network/parse error)
returns ``{"available": False}`` for that PMID and never raises.
"""

from __future__ import annotations

import json
import logging

import httpx

from ._http_retry import _with_connect_cap
from ._cache import EXTERNAL_CACHE, TTL_PMC
from ._ncbi_throttle import NCBI_CONTACT_EMAIL, NCBI_TOOL, NCBI_USER_AGENT, web_throttle


def _loads(text: str):
    """Parse JSON tolerating the unescaped control characters NCBI
    occasionally emits in E-utilities bodies (see _ncbi_throttle docstring)."""
    return json.loads(text, strict=False)

log = logging.getLogger("heartvar.pmcoa")

IDCONV = "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"
BIOC = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{}/unicode"

MAX_PMIDS = 8
MAX_TEXT_CHARS = 8000
FALLBACK_WORDS = 3000

_TRUNCATION_NOTE = " […truncated]"


async def _pmids_to_pmcids(
    client: httpx.AsyncClient, pmids: list[str]
) -> dict[str, str | None]:
    """Resolve PMIDs to their OWN PMCIDs via the NCBI ID-converter in ONE call.

    idconv accepts a comma-separated ``ids`` list and echoes the input ``pmid``
    in each returned record, so we map back by that ECHOED id — NOT by list
    position: NCBI does not guarantee record order and may drop unresolvable
    ids. ``pmid`` comes back as an integer, so both sides are coerced to str.

    We use the authoritative converter (records carry a ``pmcid`` field, or
    ``status:"error"`` when the article is not in PMC) and do NOT fall back to
    eutils elink: its ``pubmed_pmc_refs`` linkname returns the PMCIDs of
    *cited/citing* articles, not the article's own record, which would silently
    fetch the wrong paper's full text. Missing / error ids map to ``None``."""
    out: dict[str, str | None] = {p: None for p in pmids}
    if not pmids:
        return out
    await web_throttle()
    try:
        r = await client.get(
            IDCONV,
            params={"ids": ",".join(pmids), "format": "json",
                    "tool": NCBI_TOOL, "email": NCBI_CONTACT_EMAIL},
            timeout=_with_connect_cap(20.0),
        )
        if r.status_code == 200 and r.text.lstrip().startswith("{"):
            for rec in (_loads(r.text).get("records") or []):
                if rec.get("status") == "error":
                    continue
                pid = str(rec.get("pmid") or "")
                if pid in out and rec.get("pmcid"):
                    out[pid] = rec["pmcid"]
    except (httpx.HTTPError, ValueError) as e:
        log.warning("idconv batch failed for %s: %r", pmids, e)
    return out


_BIOC_ERROR_PREFIX = "[Error]"


def _bioc_has_article(body: str) -> bool:
    """Whether a BioC response body is an article rather than a not-found.

    ⚠ THIS REPLACED AN oa.fcgi PRE-CHECK THAT HAD SILENTLY KILLED THE WHOLE
    CHANNEL. ``_is_open_access`` asked
    ``https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=<PMCID>`` first and
    treated any non-200 as "not open access". NCBI retired that path, so it now
    returns **404 for every PMCID** — verified 2026-09-07 against PMC5663278,
    PMC4063552 and PMC9492974, all of which the BioC endpoint serves happily
    (75 KB for the first). Every call therefore short-circuited to
    ``available: False`` and the prompt always printed "No open-access full text
    available for variant-specific papers", starving PS3 functional-assay text
    and PS4 proband tallies of the body text they are gathered for.

    BioC is the authority on whether it can serve an article, so it is asked
    directly and no pre-flight is made. Non-OA and not-found are indistinguish-
    able here, and both mean the same thing to the caller: no body text.
    """
    return bool(body) and not body.lstrip().startswith(_BIOC_ERROR_PREFIX)


def _extract_sections(bioc: list | dict) -> tuple[str, str, str, dict]:
    """Pull (results_text, methods_text, title, meta) from a BioC JSON doc.

    BioC structure: ``[{documents:[{passages:[{infons:{section_type,type},
    text}]}]}]``. Section assignment uses ``infons.section_type`` (RESULTS /
    METHODS); when no passage carries a section_type, the caller falls back
    to a body excerpt.
    """
    docs = []
    if isinstance(bioc, list):
        for item in bioc:
            docs.extend((item or {}).get("documents") or [])
    elif isinstance(bioc, dict):
        docs = bioc.get("documents") or []

    _SKIP_FALLBACK = {
        "REF", "FIG", "ABSTRACT", "TITLE", "SUPPL",
        "COMP_INT", "AUTH_CONT", "ACK_FUND",
    }
    results: list[str] = []
    methods: list[str] = []
    tables: list[str] = []
    body: list[str] = []
    title = ""
    meta: dict = {}
    for doc in docs:
        for k in ("year", "journal"):
            v = (doc.get("infons") or {}).get(k)
            if v and k not in meta:
                meta[k] = v
        for psg in doc.get("passages") or []:
            infons = psg.get("infons") or {}
            text = (psg.get("text") or "").strip()
            if not text:
                continue
            sect = (infons.get("section_type") or "").upper()
            ptype = (infons.get("type") or "").lower()
            if sect == "TITLE" or ptype == "front":
                if not title:
                    title = text
            if not meta.get("year") and infons.get("year"):
                meta["year"] = infons["year"]
            if not meta.get("journal") and infons.get("journal"):
                meta["journal"] = infons["journal"]
            if sect == "RESULTS":
                results.append(text)
            elif sect == "METHODS":
                methods.append(text)
            elif sect == "TABLE":
                tables.append(text)
            elif sect not in _SKIP_FALLBACK:
                body.append(text)

    results_text = " ".join(results + tables).strip()
    methods_text = " ".join(methods).strip()
    if not results_text and not methods_text:
        words = " ".join(body).split()
        results_text = " ".join(words[:FALLBACK_WORDS])
    return results_text, methods_text, title, meta


def _cap_combined(results_text: str, methods_text: str) -> tuple[str, str]:
    """Cap results+methods to MAX_TEXT_CHARS combined, prioritising Results
    (where case counts / ORs typically live), with a truncation note."""
    if len(results_text) >= MAX_TEXT_CHARS:
        return results_text[:MAX_TEXT_CHARS - len(_TRUNCATION_NOTE)] + _TRUNCATION_NOTE, ""
    remaining = MAX_TEXT_CHARS - len(results_text)
    if len(methods_text) > remaining:
        methods_text = methods_text[:remaining - len(_TRUNCATION_NOTE)] + _TRUNCATION_NOTE
    return results_text, methods_text


async def _fetch_one(client: httpx.AsyncClient, pmid: str, pmcid: str | None) -> dict:
    pmid = str(pmid).strip()
    if not pmid:
        return {"pmid": pmid, "available": False, "error": "empty pmid"}
    if not pmcid:
        return {"pmid": pmid, "available": False}
    try:
        await web_throttle()
        r = await client.get(BIOC.format(pmcid), timeout=_with_connect_cap(40.0))
        if (r.status_code != 200
                or not _bioc_has_article(r.text)
                or not r.text.lstrip().startswith(("[", "{"))):
            return {"pmid": pmid, "pmcid": pmcid, "available": False}
        results_text, methods_text, title, meta = _extract_sections(r.json())
        if not results_text and not methods_text:
            return {"pmid": pmid, "pmcid": pmcid, "available": False}
        results_text, methods_text = _cap_combined(results_text, methods_text)
        return {
            "pmid": pmid,
            "pmcid": pmcid,
            "title": title,
            "year": str(meta.get("year") or ""),
            "journal": meta.get("journal") or "",
            "results_text": results_text,
            "methods_text": methods_text,
            "available": True,
        }
    except (httpx.HTTPError, ValueError) as e:
        log.warning("PMC OA fetch failed for PMID %s: %r", pmid, e)
        return {"pmid": pmid, "available": False, "error": repr(e)}


def _pmc_cacheable(result: dict) -> bool:
    """Cache a PMC result only when it's a clean, settled answer: non-empty and
    free of any transient ``error`` marker. A per-article ``available: False``
    (not in PMC / not Open Access) is a legitimate permanent answer and IS
    cached; an outer/per-article ``error`` (network/parse hiccup) is not, so it
    re-fetches next time."""
    return bool(result) and not any(
        isinstance(rec, dict) and rec.get("error") for rec in result.values()
    )


async def fetch_pmc_fulltext(pmids: list[str]) -> dict[str, dict]:
    """Cached wrapper around :func:`_fetch_pmc_fulltext_uncached`.

    PMC full text of a published article never changes, so results are cached
    for ``TTL_PMC`` keyed on the (sorted) PMID set, with single-flight
    coalescing. Empty input short-circuits without touching the cache."""
    norm = [str(p).strip() for p in (pmids or []) if str(p).strip()][:MAX_PMIDS]
    if not norm:
        return {}
    key = ("pmc", tuple(sorted(norm)))
    return await EXTERNAL_CACHE.get_or_set(
        key,
        lambda: _fetch_pmc_fulltext_uncached(norm),
        ttl=TTL_PMC,
        should_cache=_pmc_cacheable,
    )


async def _fetch_pmc_fulltext_uncached(pmids: list[str]) -> dict[str, dict]:
    """Fetch PMC OA Results/Methods excerpts for up to ``MAX_PMIDS``
    variant-specific PMIDs. Returns a dict keyed by PMID. Never raises."""
    pmids = [str(p).strip() for p in (pmids or []) if str(p).strip()][:MAX_PMIDS]
    if not pmids:
        return {}
    out: dict[str, dict] = {}
    try:
        async with httpx.AsyncClient(
            headers={"User-Agent": NCBI_USER_AGENT}, follow_redirects=True,
        ) as client:
            pmcid_map = await _pmids_to_pmcids(client, pmids)
            for pmid in pmids:
                out[pmid] = await _fetch_one(client, pmid, pmcid_map.get(pmid))
    except Exception as e:  # noqa: BLE001 — defensive: never propagate
        log.warning("fetch_pmc_fulltext outer failure: %r", e)
        for pmid in pmids:
            out.setdefault(pmid, {"pmid": pmid, "available": False, "error": repr(e)})
    return out
