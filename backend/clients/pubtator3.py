"""NCBI PubTator3 variant→literature client (PS4 retrieval).

Closes the PS4 "FIND gap" diagnosed 2026-06-05: the variant-specific PubMed
search (``pubmed.py``) is a strict ``[tiab]`` string match on the bare HGVS-c or
the three-letter amino-acid change, so it misses one-letter notation (R495W),
table-only mentions, and legacy/IVS splice forms — returning ZERO papers for
recurrent variants whose case-enrichment (PS4) evidence lives in the literature.

PubTator3 (NCBI; supersedes LitVar2) normalizes
every surface form — ``V600E`` / ``p.V600E`` / ``Val600Glu`` / ``c.1799T>A`` —
onto one tmVar3 entity, so an entity-keyed search recovers the variant's papers
regardless of nomenclature (measured ~600x more hits than an rsID/HGVS string).

Flow:
  1. ``entity/autocomplete``      — resolve (gene, AA, HGVS-c, rsID) → @VARIANT id
  2. ``search``                   — entity id → PMIDs (relevance-scored, w/ pmcid)
  3. ``publications/export``      — BioC-JSON full text → mention-anchored excerpts

This is ADDITIVE to ``pubmed.py`` / ``pmcoa.py`` — it never replaces them. The
excerpts feed the prompt as PS4 *candidate* text only; the PS4 gate is unchanged.
Non-fatal: every failure returns a result dict with empty lists, never raises.
PubTator3 is a separate NCBI service (not E-utilities), so all GETs go through
the conservative, key-independent ``web_get_json`` bucket (~3 req/s) — the
E-utilities ``api_key`` does NOT raise PubTator3's limits, and driving it at
the keyed eutils rate is what triggers a temporary block.
"""

from __future__ import annotations

import re

import httpx

from ._cache import EXTERNAL_CACHE, TTL_LITERATURE, _ok
from ._http_retry import make_async_client
from ._ncbi_throttle import NCBIError, web_get_json

PUBTATOR3_BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"

MAX_PAPERS = 20
MAX_EXPORT_PMIDS = 12
MAX_EXCERPTS = 16
MAX_RESOLVE_TERMS = 4
EXCERPT_HALF_WINDOW = 400

_AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}
_AA3_RE = re.compile(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def _one_letter_change(aa: str) -> str | None:
    """Convert a three-letter missense token ("Arg495Trp") to one-letter
    ("R495W"). Returns None for non-missense forms (Ter/fs/del/dup/splice)
    or anything that doesn't parse — those have no clean one-letter mutation
    string to search on."""
    m = _AA3_RE.match((aa or "").strip())
    if not m:
        return None
    ref, pos, alt = m.group(1), m.group(2), m.group(3)
    r1, a1 = _AA3_TO_1.get(ref), _AA3_TO_1.get(alt)
    if not r1 or not a1:
        return None
    return f"{r1}{pos}{a1}"


def _bare_c(hgvs_c: str) -> str:
    raw = (hgvs_c or "").strip()
    return raw.split(":")[-1] if ":" in raw else raw


def _candidate_terms(gene: str, hgvs_c: str, amino_acid: str | None,
                     rsid: str | None) -> list[str]:
    """Ordered autocomplete query candidates, most-specific-and-recoverable
    first. One-letter gene+AA is the single highest-value term (what the tiab
    query structurally cannot match)."""
    gene = (gene or "").strip()
    aa = (amino_acid or "").strip()
    one = _one_letter_change(aa)
    bare = _bare_c(hgvs_c)
    terms: list[str] = []
    if rsid:
        terms.append(rsid.strip())
    if gene and one:
        terms.append(f"{gene} {one}")
    if gene and aa:
        terms.append(f"{gene} {aa}")
    if gene and bare:
        terms.append(f"{gene} {bare}")
    seen: set[str] = set()
    out: list[str] = []
    for t in terms:
        if t and t not in seen:
            seen.add(t)
            out.append(t)
    return out


def _gene_matches(gene_up: str, row: dict) -> bool:
    """True only if the @VARIANT entity is for THIS gene. The id embeds the gene
    (@VARIANT_p.R495W_MYBPC3_human) and the description carries "GENE (human)".
    Matching on a word boundary avoids substring false hits (e.g. RAF1 in BRAF1)."""
    ent_id = (row.get("_id") or "").upper()
    desc = (row.get("description") or "").upper()
    return bool(re.search(rf"(^|[_ (]){re.escape(gene_up)}([_ )]|$)", ent_id)
                or re.search(rf"(^|[ (]){re.escape(gene_up)}([ )]|$)", desc))


async def _resolve_entity(client: httpx.AsyncClient, gene: str,
                          terms: list[str]) -> dict | None:
    """Try each candidate term against entity/autocomplete; return the first
    variant entity that MATCHES THIS GENE ({"_id", "db_id", "term"}).

    Requiring a gene match is essential for precision: a bare/ambiguous mutation
    token (e.g. I284F) otherwise resolves to a high-volume wrong-gene entity
    (a SARS-CoV-2 ORF1ab variant). If nothing gene-matches, return None — the
    pipeline falls back to PubMed; better no excerpts than another variant's."""
    gene_up = (gene or "").upper()
    if not gene_up:
        return None
    for term in terms[:MAX_RESOLVE_TERMS]:
        url = f"{PUBTATOR3_BASE}/entity/autocomplete/"
        params = {"query": term, "limit": 10}
        try:
            payload = await web_get_json(client, url, params, timeout=20.0)
        except NCBIError:
            continue
        rows = payload if isinstance(payload, list) else (payload.get("results") or [])
        for row in rows:
            if not isinstance(row, dict) or row.get("biotype") != "variant":
                continue
            if row.get("_id") and _gene_matches(gene_up, row):
                return {"_id": row["_id"], "db_id": row.get("db_id") or "", "term": term}
    return None


async def _search_pmids(client: httpx.AsyncClient, entity_id: str) -> list[dict]:
    """Entity-keyed literature search → relevance-sorted paper stubs."""
    url = f"{PUBTATOR3_BASE}/search/"
    params = {"text": entity_id}
    try:
        payload = await web_get_json(client, url, params, timeout=25.0)
    except NCBIError:
        return []
    results = payload.get("results") if isinstance(payload, dict) else None
    if not results:
        return []
    papers: list[dict] = []
    for r in results[:MAX_PAPERS]:
        if not isinstance(r, dict) or not r.get("pmid"):
            continue
        papers.append({
            "pmid": str(r.get("pmid")),
            "pmcid": r.get("pmcid") or None,
            "title": (r.get("title") or "").strip(),
            "journal": r.get("journal") or "",
            "year": str(r.get("date") or "")[:4],
            "score": r.get("score"),
        })
    return papers


def _iter_documents(payload):
    """Yield BioC documents from the export payload, tolerating the several
    shapes PubTator3 returns (``{"PubTator3":[...]}``, ``{"documents":[...]}``,
    a bare list, or a single document dict)."""
    if isinstance(payload, dict):
        for key in ("PubTator3", "documents", "docs"):
            if isinstance(payload.get(key), list):
                yield from payload[key]
                return
        if payload.get("passages"):
            yield payload
            return
    elif isinstance(payload, list):
        yield from payload


async def _export_excerpts(client: httpx.AsyncClient, pmids: list[str],
                           entity_id: str) -> list[dict]:
    """Pull BioC-JSON full text for ``pmids`` and window out the sentence around
    each variant mention matching ``entity_id`` — the text where proband counts
    / odds ratios live. Only entity-matched mentions are kept (precision)."""
    pmids = [p for p in pmids if p][:MAX_EXPORT_PMIDS]
    if not pmids:
        return []
    url = f"{PUBTATOR3_BASE}/publications/export/biocjson"
    params = {"pmids": ",".join(pmids), "full": "true"}
    try:
        payload = await web_get_json(client, url, params, timeout=40.0)
    except NCBIError:
        return []
    excerpts: list[dict] = []
    for doc in _iter_documents(payload):
        if not isinstance(doc, dict):
            continue
        pmid = str(doc.get("pmid") or doc.get("id") or "")
        pmcid = doc.get("pmcid")
        for passage in (doc.get("passages") or []):
            if not isinstance(passage, dict):
                continue
            text = passage.get("text") or ""
            if not text:
                continue
            p_off = passage.get("offset") or 0
            section = ((passage.get("infons") or {}).get("section_type")
                       or (passage.get("infons") or {}).get("type") or "")
            for ann in (passage.get("annotations") or []):
                infons = ann.get("infons") or {}
                if infons.get("accession") != entity_id:
                    continue
                locs = ann.get("locations") or []
                rel = (locs[0].get("offset", p_off) - p_off) if locs else 0
                start = max(0, rel - EXCERPT_HALF_WINDOW)
                end = min(len(text), rel + EXCERPT_HALF_WINDOW)
                snippet = text[start:end].strip()
                if start > 0:
                    snippet = "…" + snippet
                if end < len(text):
                    snippet = snippet + "…"
                excerpts.append({
                    "pmid": pmid, "pmcid": pmcid, "section": section,
                    "mention": ann.get("text") or "", "excerpt": snippet,
                })
                if len(excerpts) >= MAX_EXCERPTS:
                    return excerpts
                break
    return excerpts


async def fetch_pubtator3(gene: str, hgvs_c: str,
                          amino_acid: str | None = None,
                          rsid: str | None = None) -> dict:
    """Cached wrapper around :func:`_fetch_pubtator3_uncached` — keyed on
    (gene, HGVSc, amino-acid token, rsID), cached for ``TTL_LITERATURE`` with
    single-flight coalescing; only successful lookups are retained."""
    key = ("pubtator3", (gene or "").upper(), hgvs_c or "",
           amino_acid or "", rsid or "")
    return await EXTERNAL_CACHE.get_or_set(
        key,
        lambda: _fetch_pubtator3_uncached(gene, hgvs_c, amino_acid, rsid),
        ttl=TTL_LITERATURE,
        should_cache=_ok,
    )


async def _fetch_pubtator3_uncached(gene: str, hgvs_c: str,
                                    amino_acid: str | None = None,
                                    rsid: str | None = None) -> dict:
    """Resolve the variant to a PubTator3 entity, find its papers, and return
    mention-anchored full-text excerpts as PS4 candidate evidence.

    Never raises. On any failure (no entity resolved, HTTP/parse error) returns
    ``{"ok": True/False, "papers": [], "mention_excerpts": [], ...}`` so the
    pipeline proceeds exactly as before — PubMed stays the floor. Purely additive
    recall: it can never reduce what PubMed already finds.
    """
    result: dict = {
        "ok": True, "gene": gene, "hgvs_c": hgvs_c,
        "amino_acid": amino_acid or None,
        "entity_id": None, "db_id": None, "resolved_term": None,
        "query_terms_tried": [],
        "papers": [], "mention_excerpts": [],
    }
    terms = _candidate_terms(gene, hgvs_c, amino_acid, rsid)
    result["query_terms_tried"] = terms
    if not terms:
        return result
    try:
        async with make_async_client() as client:
            entity = await _resolve_entity(client, gene, terms)
            if not entity:
                result["error"] = "no PubTator3 variant entity resolved"
                return result
            result["entity_id"] = entity["_id"]
            result["db_id"] = entity["db_id"]
            result["resolved_term"] = entity["term"]
            result["papers"] = await _search_pmids(client, entity["_id"])
            export_pmids = [p["pmid"] for p in result["papers"]
                            if p.get("pmcid")] or [p["pmid"] for p in result["papers"]]
            result["mention_excerpts"] = await _export_excerpts(
                client, export_pmids, entity["_id"]
            )
    except (httpx.HTTPError, httpx.TimeoutException) as e:
        result["ok"] = False
        result["error"] = repr(e)
    return result
