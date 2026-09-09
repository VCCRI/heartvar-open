"""GenCC (Gene Curation Coalition) gene–disease classifications.

GenCC does not publish a public REST/JSON API — only a bulk submissions
TSV at https://thegencc.org/download. The snapshot is built by the MONTHLY data
job (scripts/build_gencc_db.py, which calls `_refresh_cache` here) and read from
`backend/data/gencc_submissions.json`, or from GENCC_SNAPSHOT_PATH on the mount.

⚠ LOOKUPS ARE LOCAL-ONLY, ALWAYS. The read path never downloads. It used to
refresh whenever the snapshot passed 30 days, which — for a monthly source —
meant a curation wore the whole TSV download every month. Staleness is reported
on the meta (`stale`, `age_days`) and fixed by the build, not by the request.

Each entry returns disease title, GenCC classification (Definitive /
Strong / Moderate / Limited / Disputed / Refuted / Animal Model Only),
submitting organisation, and mode of inheritance.
"""

from __future__ import annotations

import asyncio
import csv
import io
import json
import logging
import os
import time
from datetime import datetime, timezone
from pathlib import Path

from ._http_retry import _with_connect_cap, make_async_client
from .panelapp import _categorize_panel, check_hpo_relevance

DATA_FILE = Path(
    os.environ.get("GENCC_SNAPSHOT_PATH")
    or (Path(__file__).resolve().parent.parent / "data" / "gencc_submissions.json")
)
log = logging.getLogger("heartvar.gencc")

SOURCE_TSV = "https://thegencc.org/download/action/submissions-export-tsv?format=new"
CACHE_TTL_SECONDS = 30 * 24 * 3600

_INDEX: dict[str, list[dict]] | None = None
_META: dict | None = None
_LOAD_LOCK = asyncio.Lock()


def _parse_tsv(tsv_text: str) -> tuple[dict[str, list[dict]], dict]:
    """Group TSV rows by uppercased gene_symbol. Each submission carries
    the row's ``gene_curie`` ("HGNC:NNNN…") so the frontend can build a
    direct gene-page URL — search.thegencc.org rejects symbol-keyed URLs
    and only resolves the curie form.
    """
    reader = csv.DictReader(io.StringIO(tsv_text), delimiter="\t")
    by_gene: dict[str, list[dict]] = {}
    rows_in = 0
    for row in reader:
        sym = (row.get("gene_symbol") or "").strip().upper()
        if not sym:
            continue
        rows_in += 1
        by_gene.setdefault(sym, []).append({
            "disease": row.get("disease_title") or row.get("submitted_as_disease_name") or "",
            "disease_curie": row.get("disease_curie") or "",
            "classification": row.get("classification_title") or "",
            "moi": row.get("moi_title") or "",
            "submitter": row.get("submitter_title") or "",
            "pmids": row.get("submitted_as_pmids") or "",
            "submission_date": row.get("submitted_as_date") or "",
            "hgnc_id": (row.get("gene_curie") or "").strip(),
        })
    meta = {
        "source_url": SOURCE_TSV,
        "fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "row_count": rows_in,
        "gene_count": len(by_gene),
    }
    return by_gene, meta


async def _refresh_cache() -> tuple[dict[str, list[dict]], dict]:
    async with make_async_client() as c:
        r = await c.get(SOURCE_TSV, timeout=_with_connect_cap(60.0), follow_redirects=True)
        r.raise_for_status()
        index, meta = _parse_tsv(r.text)
    DATA_FILE.parent.mkdir(parents=True, exist_ok=True)
    DATA_FILE.write_text(json.dumps({"meta": meta, "genes": index}, indent=2) + "\n")
    return index, meta


async def _ensure_loaded() -> tuple[dict[str, list[dict]] | None, dict]:
    """The snapshot, or ``(None, {"error": ...})``. NEVER a live download.

    ⚠ THIS USED TO FETCH thegencc.org MID-CURATION. The old rule was "serve the
    local snapshot if it is under CACHE_TTL_SECONDS (30 days) old, otherwise
    refresh". GenCC is a MONTHLY-cadence source — build_all.sh:488 builds it and
    :492 copies it to the mount — so the snapshot crosses that boundary EVERY
    MONTH by construction, and the next curation after it did paid for the whole
    TSV download inside a user's request.

    A stale snapshot is a BUILD problem, and fixing it here is the wrong place
    three times over: it makes one unlucky curation wear a multi-second
    download; it does so at the moment the data job has already failed, i.e.
    when the network is the least trustworthy thing in the system; and it hides
    the build failure that caused it. HEARTVAR_OFFLINE_STRICT would have
    prevented it, but that flag is not set in the deploy and the guarantee
    should not depend on an operator remembering it.

    So: local only, always, and staleness is REPORTED (``stale``, ``age_days``
    on the meta) rather than repaired. ``_refresh_cache`` remains the single
    parser — scripts/build_gencc_db.py calls it directly, which is what keeps
    the snapshot current.
    """
    global _INDEX, _META
    if _INDEX is not None and _META is not None:
        return _INDEX, _META
    async with _LOAD_LOCK:
        if _INDEX is not None and _META is not None:
            return _INDEX, _META
        if not DATA_FILE.exists():
            msg = (
                f"GenCC snapshot missing at {DATA_FILE} — the monthly data "
                "build (scripts/build_gencc_db.py) has not produced it. GenCC "
                "evidence is unavailable for this curation; it is NOT fetched "
                "live."
            )
            log.warning("%s", msg)
            _META = {"error": msg}
            _INDEX = None
            return None, _META
        try:
            payload = json.loads(DATA_FILE.read_text())
        except (OSError, ValueError) as exc:
            msg = f"GenCC snapshot at {DATA_FILE} is unreadable: {exc!r}"
            log.warning("%s", msg)
            _META = {"error": msg}
            _INDEX = None
            return None, _META
        age_days = (time.time() - DATA_FILE.stat().st_mtime) / 86400.0
        _INDEX = payload.get("genes") or {}
        _META = {**(payload.get("meta") or {}), "age_days": round(age_days, 1)}
        if age_days * 86400 >= CACHE_TTL_SECONDS:
            _META["stale"] = True
            log.warning(
                "GenCC snapshot is %.0f days old (>%d days) — serving it anyway; "
                "the monthly data build should refresh it.",
                age_days, CACHE_TTL_SECONDS // 86400,
            )
        return _INDEX, _META


_DISPUTED_TIERS = ("Disputed Evidence", "Refuted Evidence")


def _aggregate(subs: list[dict]) -> tuple[str | None, bool]:
    """Return (best_classification, has_disputed_or_refuted) for a submission list.
    Assumes `subs` is already sorted by tier rank ascending."""
    best = subs[0]["classification"] if subs else None
    disputed = any(s.get("classification") in _DISPUTED_TIERS for s in subs)
    return best, disputed


async def fetch_gencc(gene: str, submitted_hpo=None) -> dict:
    """Look up GenCC gene-disease curations, tagging each submission with
    a cardiovascular category derived from its disease title and a
    proband-HPO match flag. The phenotype-matched aggregates are what the
    UI and prompt should read for the headline — the gene-wide aggregates
    are kept for completeness but can be misleading on pleiotropic genes
    (e.g. RAF1: Definitive Noonan + Disputed CFC + Strong DCM-1NN — the
    gene-wide "Disputed" flag is wrong for a DCM proband)."""
    index, meta = await _ensure_loaded()
    if index is None:
        return {"ok": False, "error": meta.get("error", "GenCC unavailable")}
    submissions = index.get(gene.strip().upper()) or []
    if not submissions:
        return {
            "ok": True,
            "gene": gene,
            "found": False,
            "submissions": [],
            "fetched_at": meta.get("fetched_at"),
        }
    rank = [
        "Definitive", "Strong", "Moderate", "Limited",
        "Animal Model Only", "Disputed Evidence", "Refuted Evidence",
    ]
    rank_index = {r: i for i, r in enumerate(rank)}
    sorted_subs = sorted(
        submissions,
        key=lambda s: rank_index.get(s.get("classification", ""), 99),
    )
    for s in sorted_subs:
        cat = _categorize_panel(s.get("disease") or "")
        s["category"] = cat
        s["hpo_match"] = bool(cat) and check_hpo_relevance(submitted_hpo, cat) \
            if submitted_hpo else False
    best, has_disputed = _aggregate(sorted_subs)
    pheno_subs = [s for s in sorted_subs if s.get("hpo_match")]
    pheno_best, pheno_has_disputed = _aggregate(pheno_subs)
    clingen_subs = [
        s for s in sorted_subs
        if (s.get("submitter") or "").strip().lower() == "clingen"
    ]
    best_clingen, clingen_has_disputed = _aggregate(clingen_subs)
    clingen_pheno = [s for s in clingen_subs if s.get("hpo_match")]
    clingen_pheno_best, clingen_pheno_has_disputed = _aggregate(clingen_pheno)
    hgnc_id = next(
        (s.get("hgnc_id") for s in sorted_subs if s.get("hgnc_id")),
        "",
    )
    url = (
        f"https://thegencc.org/genes/{hgnc_id}"
        if hgnc_id
        else f"https://search.thegencc.org/genes/HGNC?q={gene}"
    )
    return {
        "ok": True,
        "gene": gene,
        "hgnc_id": hgnc_id or None,
        "found": True,
        "best_classification": best,
        "has_disputed_or_refuted": has_disputed,
        "submission_count": len(submissions),
        "submissions": sorted_subs,
        "phenotype_matched_count": len(pheno_subs),
        "phenotype_matched_best_classification": pheno_best,
        "phenotype_matched_has_disputed_or_refuted": pheno_has_disputed,
        "submitted_hpo": submitted_hpo,
        "clingen_submissions": clingen_subs,
        "clingen_submission_count": len(clingen_subs),
        "best_clingen_classification": best_clingen,
        "clingen_has_disputed_or_refuted": clingen_has_disputed,
        "clingen_phenotype_matched_count": len(clingen_pheno),
        "clingen_phenotype_matched_best_classification": clingen_pheno_best,
        "clingen_phenotype_matched_has_disputed_or_refuted": clingen_pheno_has_disputed,
        "url": url,
        "fetched_at": meta.get("fetched_at"),
    }
