"""ClinGen Evidence Repository (eRepo) client — READ-ONLY CURATOR REFERENCE.

Surfaces the VCEP's published classification + applied ACMG criteria for a
variant when ClinGen has already curated it. This is shown to the human
curator as an independent reference — it is **never** passed to the LLM as
evidence and never influences HeartVar's own criteria.

CRITICAL — it gives the literal answer (classification + criteria codes), so
it is a reference for the curator, not an input to automated reasoning: eRepo
data is never passed to the LLM or used by HeartVar's own deterministic
criteria.

LOCAL-FIRST lookup: the deploy downloads the FULL eRepo dump to
``backend/data/erepo_all.tsv`` (``scripts/build_erepo_dump.py``, refreshed
~monthly). That dump carries the same fields the live API returns
(classification, met-criteria, expert panel, link) for ALL ~12.6k records
across every VCEP — which is *strictly more complete* than the live
``?gene=`` endpoint, which caps at ~25 records/gene and ignores offset
paging (so a live ``found: False`` can be a false negative). We therefore
read the local TSV when present and only fall back to the live API when the
TSV is absent. This removes the per-curation ``erepo.clinicalgenome.org``
call entirely on a normal deploy. Path overridable via ``HEARTVAR_EREPO_TSV``.

API notes (from the 2026-06 source investigation):
  * Public, no auth: GET .../api/classifications?gene={gene} (JSON).
  * The endpoint returns at most ~25 records per gene and does NOT honour
    offset/paging params — so live coverage is best-effort.

Non-fatal: any failure returns {"ok": False, ...} and never raises.
"""

from __future__ import annotations

import asyncio
import csv
import logging
import os
import re
from pathlib import Path

import httpx
from ._http_retry import _with_connect_cap
from ..localio import run_local

log = logging.getLogger("heartvar.erepo")

EREPO_API = "https://erepo.clinicalgenome.org/evrepo/api/classifications"

_DEFAULT_TSV = Path(__file__).resolve().parent.parent / "data" / "erepo_all.tsv"


def _tsv_path() -> Path:
    override = os.environ.get("HEARTVAR_EREPO_TSV", "").strip()
    return Path(override) if override else _DEFAULT_TSV


_MIN_INTERVAL = 1.0
_lock = asyncio.Lock()
_last_call = 0.0

_cache: dict[tuple[str, str], dict] = {}

_tsv_index: dict[tuple[str, str], dict] | None = None
_tsv_index_loaded = False
_tsv_lock = asyncio.Lock()

_TX_PREFIX_RE = re.compile(r"^[A-Za-z]+_?\d+(?:\.\d+)?(?:\([^)]+\))?:")


def _strip_tx(hgvs: str) -> str:
    """Reduce an HGVS string to its bare ``c.…`` form, lowercased, for
    transcript-agnostic matching (``NM_005343.3(HRAS):c.34G>A`` → ``c.34g>a``)."""
    s = (hgvs or "").strip()
    s = _TX_PREFIX_RE.sub("", s)
    s = s.split(" ", 1)[0]
    return s.lower()


def _is_retracted(value: str) -> bool:
    return (value or "").strip().lower() in {"true", "1", "yes"}


def _result_from_tsv_row(row: dict) -> dict:
    """Build the same shape ``_extract`` returns, from one TSV row."""
    classification = (row.get("Assertion") or "").strip()
    vcep = (row.get("Expert Panel") or "").strip()
    url = (row.get("Evidence Repo Link") or "").strip()
    raw_criteria = (row.get("Applied Evidence Codes (Met)") or "").strip()
    criteria: list[str] = []
    seen: set[str] = set()
    for part in raw_criteria.split(","):
        code = part.strip()
        if code and code not in seen:
            seen.add(code)
            criteria.append(code)
    return {
        "ok": True,
        "found": True,
        "classification": classification,
        "criteria": criteria,
        "vcep": vcep,
        "url": url,
    }


_VARIATION_GENE_RE = re.compile(r"\(([A-Za-z0-9_.\-]+)\)\s*:")


def _variation_gene(row: dict) -> str:
    """Gene symbol from the row's own HGVS, independent of its label column."""
    m = _VARIATION_GENE_RE.search(row.get("#Variation") or "")
    return (m.group(1).strip().upper() if m else "")


def _build_tsv_index(path: Path) -> dict[tuple[str, str], dict]:
    """Parse the full eRepo dump into a (gene_upper, stripped_hgvs) index.

    Every HGVS expression a record carries is indexed (transcript-agnostic),
    mirroring the live client's ``any(_strip_tx(h) == target …)`` match.
    Retracted rows are skipped. The FIRST record to claim a given key wins
    (the dump is ordered newest-first per gene, matching the live "first
    matching record" behaviour).

    INDEXED UNDER BOTH THE LABEL AND THE HGVS's OWN GENE, because the two
    disagree on 66 of the export's 12,702 rows and the label is not always the
    one a caller will ask with. The case that matters here: eRepo labels 11
    records ``LRRC56`` whose HGVS is written on ``NM_005343.3(HRAS)`` — the two
    genes are adjacent at 11p15.5 — and all 11 are RASopathy-VCEP curations of
    HRAS variants. Keyed on the label alone, a curation of HRAS never finds
    them, and they also fall through to the generic frequency path: BA1 at 0.05
    instead of the RASopathy spec's 0.0005, a hundred times looser, plus a PM2
    that fires on present variants where that spec is absence-only.

    Both keys rather than a replacement, deliberately. Most of the other 55
    disagreements are mitochondrial or overlapping-gene rows (MT-TL1/MT-CYB,
    KLLN/PTEN, CDKL5/RS1) where the LABEL is plausibly the curated choice and
    the HGVS's gene is just whichever transcript the expression was written on.
    Adding an alias cannot lose a lookup that works today; replacing the symbol
    could."""
    index: dict[tuple[str, str], dict] = {}
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if _is_retracted(row.get("Retracted", "")):
                continue
            genes = {
                g for g in (
                    (row.get("HGNC Gene Symbol") or "").strip().upper(),
                    _variation_gene(row),
                ) if g and g != "N/A"
            }
            if not genes:
                continue
            result = _result_from_tsv_row(row)
            for expr in (row.get("HGVS Expressions") or "").split(","):
                stripped = _strip_tx(expr)
                if not stripped:
                    continue
                for gene in genes:
                    index.setdefault((gene, stripped), result)
    return index


async def _ensure_tsv_index() -> dict[tuple[str, str], dict] | None:
    """Load the local TSV index once. Returns None when the dump is absent
    (so the caller falls back to the live API)."""
    global _tsv_index, _tsv_index_loaded
    if _tsv_index_loaded:
        return _tsv_index
    async with _tsv_lock:
        if _tsv_index_loaded:
            return _tsv_index
        path = _tsv_path()
        if not path.is_file():
            _tsv_index = None
            _tsv_index_loaded = True
            log.info("eRepo local dump not found at %s; using live API fallback", path)
            return None
        try:
            _tsv_index = await run_local(_build_tsv_index, path)
            log.info("eRepo local dump loaded: %d (gene,hgvs) keys from %s",
                     len(_tsv_index), path)
        except (OSError, csv.Error, ValueError) as e:
            log.warning("eRepo local dump unreadable (%r); falling back to live API", e)
            _tsv_index = None
        _tsv_index_loaded = True
        return _tsv_index


async def _throttle() -> None:
    global _last_call
    async with _lock:
        import time
        wait = (_last_call + _MIN_INTERVAL) - time.monotonic()
        if wait > 0:
            await asyncio.sleep(wait)
        _last_call = time.monotonic()


def _extract(vi: dict) -> dict:
    """Pull classification / met-criteria / VCEP / url from one
    variantInterpretation record (LIVE API shape)."""
    classification = ""
    vcep = ""
    criteria: list[str] = []
    for g in vi.get("guidelines") or []:
        g_outcome = (g.get("outcome") or {}).get("label")
        if g_outcome and not classification:
            classification = g_outcome
        for a in g.get("agents") or []:
            a_outcome = (a.get("outcome") or {}).get("label")
            if a_outcome and not classification:
                classification = a_outcome
            if a.get("affiliation") and not vcep:
                vcep = a["affiliation"]
            for ec in a.get("evidenceCodes") or []:
                if ec.get("status") == "Met" and ec.get("label"):
                    criteria.append(ec["label"])
    at_id = vi.get("@id") or ""
    url = at_id.replace("/api/interpretation/", "/ui/classification/") if at_id else ""
    seen: set[str] = set()
    uniq = [c for c in criteria if not (c in seen or seen.add(c))]
    return {
        "ok": True,
        "found": True,
        "classification": classification,
        "criteria": uniq,
        "vcep": vcep,
        "url": url,
    }


async def _fetch_erepo_live(gene: str, target: str, hgvs_c: str) -> dict:
    """Live ClinGen eRepo API lookup (fallback when the local dump is absent).

    Used only when ``backend/data/erepo_all.tsv`` is not present; on a normal
    deploy the local dump is authoritative and this never runs."""
    try:
        seen_caids: set[str] = set()
        match: dict | None = None
        async with httpx.AsyncClient(
            follow_redirects=True,
            headers={"Accept": "application/json",
                     "User-Agent": "HeartVar/1.0 (mailto:heartvar@victorchang.edu.au)"},
        ) as client:
            for page in range(3):
                await _throttle()
                r = await client.get(
                    EREPO_API,
                    params={"gene": gene, "_offset": page * 25},
                    timeout=_with_connect_cap(25.0),
                )
                if r.status_code != 200 or not r.text.lstrip().startswith("{"):
                    break
                vis = (r.json().get("variantInterpretations") or [])
                new = 0
                for vi in vis:
                    caid = vi.get("caid") or vi.get("@id") or ""
                    if caid in seen_caids:
                        continue
                    seen_caids.add(caid)
                    new += 1
                    if any(_strip_tx(h) == target for h in (vi.get("hgvs") or [])):
                        match = vi
                        break
                if match or new == 0:
                    break
        return _extract(match) if match else {"ok": True, "found": False}
    except (httpx.HTTPError, ValueError) as e:
        log.warning("eRepo live lookup failed for %s %s: %r", gene, hgvs_c, e)
        return {"ok": False, "found": False, "error": repr(e)}


async def fetch_erepo(gene: str, hgvs_c: str) -> dict:
    """Return the ClinGen eRepo VCEP verdict for ``gene``/``hgvs_c``.

    Reads the local full-dump TSV (``backend/data/erepo_all.tsv``) when
    present — authoritative and more complete than the live per-gene API —
    and only falls back to the live ``erepo.clinicalgenome.org`` API when the
    dump is absent. READ-ONLY reference — never feed the result to the LLM.
    Returns ``{"ok": True, "found": False}`` when no match,
    ``{"ok": False, "error": ...}`` on failure. Never raises."""
    gene = (gene or "").strip()
    target = _strip_tx(hgvs_c)
    if not gene or not target:
        return {"ok": False, "found": False, "error": "gene and hgvs_c required"}

    key = (gene.upper(), target)
    if key in _cache:
        return _cache[key]

    index = await _ensure_tsv_index()
    if index is not None:
        result = index.get(key, {"ok": True, "found": False})
    else:
        result = await _fetch_erepo_live(gene.upper(), target, hgvs_c)

    _cache[key] = result
    return result
