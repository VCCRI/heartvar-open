"""PubMed literature search for HeartVar evidence gathering.

Runs two searches via NCBI E-utilities:
  Search A — variant-specific (gene + HGVS-c and/or amino-acid change)
  Search B — gene + cardiac-disease MeSH + functional-evidence keywords

For each PMID returned, fetches PMID, title, first author, year, journal,
abstract (truncated to 300 words), and DOI.

Notes
-----
* The MeSH terms in the user-supplied query template ("congenital heart
  disease", "cardiomyopathy") are not valid MeSH descriptors. The correct
  forms are "Heart Defects, Congenital" and "Cardiomyopathies". "Channelopathy"
  is dropped — it's not a routine CHD MeSH and was returning 0 results.
* All failures are non-fatal: any HTTP error or parse failure returns empty
  result lists, never raises.
* All E-utilities requests go through the process-wide throttle in
  ``_ncbi_throttle`` so the variant-paper, gene-literature, and MedGen
  callers share a single 3 req/s budget (10 req/s with ``NCBI_API_KEY``).
  Without that shared limiter the first-wave esearches fired concurrently
  by ``app.py`` would hit NCBI's 429 rate-cap and surface as silent empty
  results.
"""

from __future__ import annotations

from xml.etree import ElementTree as ET

import httpx

from ._cache import EXTERNAL_CACHE, TTL_LITERATURE, _ok
from ._ncbi_throttle import (
    NCBI_API_KEY,
    NCBI_CONTACT_EMAIL,
    NCBI_TOOL,
    NCBI_USER_AGENT,
    NCBIError,
    eutils_get_json,
    eutils_get_text,
)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

ABSTRACT_WORD_CAP = 350
GENE_LIT_SNIPPET_CHARS = 300
GENE_LIT_RETMAX = 10
PER_SEARCH_RETMAX = 10
VARIANT_TOTAL_CAP = 10


def _eutils_params(**extra) -> dict:
    """Build E-utilities query params with the NCBI-policy fields attached
    (tool + contact email — see _ncbi_throttle), plus the API key when set.
    Covers both the variant search/efetch and the gene-literature paths since
    all of them route through here."""
    params = dict(extra)
    params["tool"] = NCBI_TOOL
    params["email"] = NCBI_CONTACT_EMAIL
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    return params

_DEFAULT_CARDIAC_KEYWORDS: list[str] = [
    "congenital heart disease",
    "cardiomyopathy",
    "cardiac",
]


def _abstract_truncate(text: str, words: int = ABSTRACT_WORD_CAP) -> str:
    toks = (text or "").split()
    if len(toks) <= words:
        return text or ""
    return " ".join(toks[:words]) + " …"


def _parse_pubmed_xml(xml_text: str) -> list[dict]:
    """Parse efetch's PubmedArticleSet XML into a list of paper dicts.

    Tolerant of missing fields and structured AbstractText (Background /
    Methods / Results / Conclusion labels)."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []

    out: list[dict] = []
    for art in root.findall(".//PubmedArticle"):
        pmid = art.findtext(".//PMID") or ""
        title = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""

        first_author = ""
        first = art.find(".//AuthorList/Author")
        if first is not None:
            last = first.findtext("LastName") or ""
            init = first.findtext("Initials") or ""
            collective = first.findtext("CollectiveName") or ""
            first_author = f"{last} {init}".strip() if last else collective

        year = (art.findtext(".//ArticleDate/Year")
                or art.findtext(".//PubDate/Year")
                or art.findtext(".//PubDate/MedlineDate") or "")[:4]

        journal = (art.findtext(".//Journal/ISOAbbreviation")
                   or art.findtext(".//Journal/Title") or "")

        abstract_parts: list[str] = []
        for at in art.findall(".//Abstract/AbstractText"):
            chunk = "".join(at.itertext()).strip()
            if not chunk:
                continue
            label = at.attrib.get("Label", "")
            abstract_parts.append(f"{label}: {chunk}" if label else chunk)
        abstract = _abstract_truncate(" ".join(abstract_parts))

        doi = ""
        for aid in art.findall(".//ArticleId"):
            if aid.attrib.get("IdType") == "doi" and aid.text:
                doi = aid.text.strip()
                break

        out.append({
            "pmid": pmid,
            "title": title.strip(),
            "first_author": first_author,
            "year": year,
            "journal": journal.strip(),
            "abstract": abstract,
            "doi": doi,
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
        })
    return out


async def _esearch(client: httpx.AsyncClient, query: str) -> tuple[int, list[str]]:
    """Run an esearch and return ``(total, idlist)``.

    ``total`` is the true PubMed hit count (``esearchresult.count``), which
    can exceed ``len(idlist)`` because ``idlist`` is capped at
    ``PER_SEARCH_RETMAX``. The frontend needs the true total for the
    displayed paper count, while the (capped) idlist drives the efetch.
    On any NCBI error returns ``(0, [])``.
    """
    params = _eutils_params(
        db="pubmed",
        term=query,
        retmode="json",
        retmax=str(PER_SEARCH_RETMAX),
        sort="relevance",
    )
    try:
        payload = await eutils_get_json(
            client, f"{EUTILS}/esearch.fcgi", params, timeout=20.0
        )
    except NCBIError:
        return (0, [])
    esearchresult = (payload.get("esearchresult") or {})
    idlist = esearchresult.get("idlist") or []
    try:
        total = int(esearchresult.get("count"))
    except (TypeError, ValueError):
        total = len(idlist)
    return (total, idlist)


async def _efetch(client: httpx.AsyncClient, pmids: list[str]) -> list[dict]:
    if not pmids:
        return []
    params = _eutils_params(
        db="pubmed",
        id=",".join(pmids),
        rettype="abstract",
        retmode="xml",
    )
    try:
        xml_text = await eutils_get_text(
            client, f"{EUTILS}/efetch.fcgi", params, timeout=25.0
        )
    except NCBIError:
        return []
    return _parse_pubmed_xml(xml_text)


def _strip_transcript_prefix(hgvs_c: str) -> str:
    """Drop ``"ENST…:"`` / ``"NM_…:"`` so the bare ``c.…`` form is what
    PubMed sees. Papers cite ``c.703C>T`` — never the transcript-prefixed
    form — so leaving the prefix on returned zero hits.
    """
    raw = (hgvs_c or "").strip()
    if ":" in raw:
        return raw.split(":")[-1]
    return raw


def _build_variant_query(gene: str, hgvs_c: str, amino_acid: str | None) -> str:
    """One canonical variant-specific PubMed query.

    Form: ``{gene}[tiab] AND ("{hgvsc}"[tiab] OR "{aa}"[tiab])`` (the
    ``OR "{aa}"`` half is omitted when no protein change is known). No
    pathogenic/case-report gate and no residue-position fallback —
    decisions reached with the curator are that the displayed list and
    the click-through link must match exactly, so the query stays
    narrow and surfaces only papers that literally mention either the
    HGVS-c or the amino-acid change.

    The ``[gene]`` field tag was previously used here. PubMed has no
    ``[gene]`` index for the article database (it is silently stripped,
    triggering the "unknown field was ignored" banner). ``[tiab]``
    searches title + abstract, which is the right surrogate for "papers
    that discuss this gene".
    """
    aa = (amino_acid or "").strip()
    hgvs_bare = _strip_transcript_prefix(hgvs_c)
    parts: list[str] = []
    if hgvs_bare:
        parts.append(f'"{hgvs_bare}"[tiab]')
    if aa:
        parts.append(f'"{aa}"[tiab]')
    if not parts:
        return ""
    or_clause = parts[0] if len(parts) == 1 else "(" + " OR ".join(parts) + ")"
    return f'{gene}[tiab] AND {or_clause}'


def _build_gene_disease_query(gene: str, disease_phrases: list[str]) -> str:
    """One canonical gene-disease PubMed query.

    Form: ``{gene}[tiab] AND ("phrase1"[tiab] OR "phrase2"[tiab] …)``
    with ``disease_phrases`` already resolved to plain-text phenotype
    labels by the caller (HPO IDs translated, free-text passed through,
    deduped). When no phrases are available the cardiac fallback
    keywords are used so the search remains anchored on cardiac
    literature.

    No MeSH clause — decided with the curator that free-text title/
    abstract phrases are easier to reason about than MeSH, and lenient
    enough to catch papers that haven't been fully MeSH-indexed yet.
    """
    phrases = [p.strip() for p in (disease_phrases or []) if p and p.strip()]
    if not phrases:
        phrases = list(_DEFAULT_CARDIAC_KEYWORDS)
    or_clause = " OR ".join(f'"{p}"[tiab]' for p in phrases)
    return f'{gene}[tiab] AND ({or_clause})'


def build_pubmed_search_url(query: str) -> str:
    """Return the human-facing PubMed search URL for ``query``, sorted
    by "Best Match" (relevance). Backend ships this string to the
    frontend so the click-through link is *exactly* the same search
    the displayed paper list came from — no syntax drift between the
    pre-fetched results and what the curator sees on click."""
    import urllib.parse as _u
    if not query:
        return "https://pubmed.ncbi.nlm.nih.gov/"
    return (
        "https://pubmed.ncbi.nlm.nih.gov/?term="
        + _u.quote(query, safe="")
        + "&sort=relevance"
    )


async def get_gene_disease_literature(
    gene_symbol: str,
    phenotype_keywords: list[str] | None = None,
) -> dict:
    """Cached wrapper around :func:`_get_gene_disease_literature_uncached`.

    Keyed on (gene, phenotype phrases) and cached for ``TTL_LITERATURE`` with
    single-flight coalescing, so concurrent or repeat curations of the same
    gene+phenotype share one set of NCBI E-utilities calls. Only successful
    lookups are retained."""
    key = ("gene_lit", (gene_symbol or "").upper(),
           tuple(phenotype_keywords or ()))
    return await EXTERNAL_CACHE.get_or_set(
        key,
        lambda: _get_gene_disease_literature_uncached(gene_symbol, phenotype_keywords),
        ttl=TTL_LITERATURE,
        should_cache=_ok,
    )


async def _get_gene_disease_literature_uncached(
    gene_symbol: str,
    phenotype_keywords: list[str] | None = None,
) -> dict:
    """Search PubMed for gene-level literature anchored on the proband's
    phenotype (free-text title/abstract phrases — no MeSH).

    Single canonical query:
      ``{gene}[tiab] AND ("phrase1"[tiab] OR "phrase2"[tiab] ...)``

    ``phenotype_keywords`` should already be the human-readable disease
    phrases (resolved HPO labels and/or curator-typed free text). Empty
    input falls back to the broad cardiac vocabulary so the search stays
    anchored on cardiac literature.

    Results are sorted by PubMed "Best Match" relevance, capped at
    ``GENE_LIT_RETMAX``. The ``search_url`` field on the result is the
    exact PubMed URL the frontend should use for the click-through link,
    so paper #1 in the panel matches paper #1 on PubMed.
    """
    gene = (gene_symbol or "").strip()
    if not gene:
        return {
            "ok": False,
            "gene": gene,
            "phenotype_keywords": [],
            "query": "",
            "search_url": "",
            "papers": [],
            "error": "gene_symbol is required",
        }

    keywords = [k for k in (phenotype_keywords or []) if k and k.strip()]
    if not keywords:
        keywords = list(_DEFAULT_CARDIAC_KEYWORDS)
    query = _build_gene_disease_query(gene, keywords)
    search_url = build_pubmed_search_url(query)

    result: dict = {
        "ok": True,
        "gene": gene,
        "phenotype_keywords": keywords,
        "query": query,
        "search_url": search_url,
        "papers": [],
    }

    try:
        async with httpx.AsyncClient(headers={"User-Agent": NCBI_USER_AGENT}) as client:
            params = _eutils_params(
                db="pubmed",
                term=query,
                retmode="json",
                retmax=str(GENE_LIT_RETMAX),
                sort="relevance",
            )
            try:
                payload = await eutils_get_json(
                    client, f"{EUTILS}/esearch.fcgi", params, timeout=20.0
                )
            except NCBIError as e:
                result["ok"] = False
                result["error"] = f"esearch {e}"
                return result
            pmids = (
                (payload.get("esearchresult") or {}).get("idlist") or []
            )
            if not pmids:
                return result
            papers = await _efetch(client, pmids)
    except (httpx.HTTPError, httpx.TimeoutException) as e:
        return {**result, "ok": False, "error": repr(e), "papers": []}

    for p in papers:
        abstract = p.get("abstract") or ""
        snippet = abstract[:GENE_LIT_SNIPPET_CHARS]
        if len(abstract) > GENE_LIT_SNIPPET_CHARS:
            snippet = snippet.rstrip() + "…"
        result["papers"].append({**p, "abstract_snippet": snippet})

    return result


async def fetch_pubmed(gene: str, hgvs_c: str, amino_acid: str | None = None) -> dict:
    """Cached wrapper around :func:`_fetch_pubmed_uncached` — keyed on
    (gene, HGVSc, amino-acid token), cached for ``TTL_LITERATURE`` with
    single-flight coalescing; only successful lookups are retained."""
    key = ("pubmed", (gene or "").upper(), hgvs_c or "", amino_acid or "")
    return await EXTERNAL_CACHE.get_or_set(
        key,
        lambda: _fetch_pubmed_uncached(gene, hgvs_c, amino_acid),
        ttl=TTL_LITERATURE,
        should_cache=_ok,
    )


async def _fetch_pubmed_uncached(gene: str, hgvs_c: str, amino_acid: str | None = None) -> dict:
    """Run the canonical variant-specific PubMed search and return the
    top-N results (sorted by PubMed "Best Match" relevance).

    ``variant_papers`` is the capped top-N list, but ``variant_total`` is
    the TRUE esearch hit count and may exceed ``len(variant_papers)``. The
    frontend should display ``variant_total`` as the paper count so it
    matches the click-through ``variant_search_url`` (which is uncapped),
    while still rendering only the top-N papers in the list.

    The returned ``variant_search_url`` is the exact human-facing URL
    PubMed served the pre-fetched papers from — the frontend uses it
    verbatim for the "Search PubMed for this variant" click-through,
    guaranteeing that paper #1 in the displayed list is paper #1 on
    the click-through page. No more drift between what the curator
    sees in the card and what they see on PubMed.

    Gene-level papers are no longer pulled here — the Gene-tab "Gene
    literature" panel runs its own phenotype-anchored search via
    :func:`get_gene_disease_literature`.

    Never raises. On any HTTP, network, or parse failure, returns the
    result dict with an empty paper list and an `error` field populated.
    """
    variant_query = _build_variant_query(gene, hgvs_c, amino_acid)
    variant_search_url = build_pubmed_search_url(variant_query)

    result: dict = {
        "ok": True,
        "gene": gene,
        "hgvs_c": hgvs_c,
        "amino_acid": amino_acid or None,
        "variant_query": variant_query,
        "variant_search_url": variant_search_url,
        "variant_papers": [],
        "variant_total": 0,
    }

    if not variant_query:
        return result

    try:
        async with httpx.AsyncClient(headers={"User-Agent": NCBI_USER_AGENT}) as client:
            total, idlist = await _esearch(client, variant_query)
            result["variant_total"] = total
            pmids = idlist[:VARIANT_TOTAL_CAP]
            if pmids:
                result["variant_papers"] = await _efetch(client, pmids)
    except (httpx.HTTPError, httpx.TimeoutException) as e:
        result["ok"] = False
        result["error"] = repr(e)

    return result
