"""Local-DB-backed ClinVar lookup.

Reads from ``data/clinvar.db`` built by ``build_clinvar_db.py``. Replaces
the previous NCBI E-utilities client: no network calls, no rate-limit
retries. Return shape is the same as the legacy client so the curate
pipeline picks up the swap without any downstream changes.

Rebuild the database monthly:

    python3 scripts/build_clinvar_db.py
"""

from __future__ import annotations

import logging
import re
import sqlite3
from time import perf_counter

from ..localio import connect_ro, log_safe, run_local
from ._paths import db_path

log = logging.getLogger("heartvar.clinvar")


_DB_MISSING = {
    "ok": False,
    "error": (
        "The local ClinVar database could not be opened (the data mount may be "
        "unavailable). No ClinVar evidence was retrieved."
    ),
}


def _conn(*, rows: bool):
    """A REUSED read-only connection, with ``row_factory`` set explicitly.

    WHY REUSE. Every query here used to open its own ``sqlite3.connect``, six
    sites of them, so each one started with an EMPTY page cache and re-read the
    index and table pages it needed from the Azure Files SMB mount. clinvar.db
    is the largest database we ship and its gene queries carry
    ``ORDER BY number_submitters DESC``, which means reading every row for the
    gene before ``LIMIT 100`` can apply. Measured 2026-08-31 with the fixed
    per-source instrument: the `clinvar` task took 14.63 s on CHD7, 10.22 s on
    MYH7 and 7.75 s on MYBPC3, against 0.06-0.30 s for every other local
    source. It scales with the gene's record count, which is the signature of
    re-reading pages rather than of query complexity.

    ``localio.connect_ro`` caches one connection per thread per database, so the
    page cache survives between queries and between curations. It was added in
    the PR #53-#58 sweep and clinvar — the database that needed it most — was
    never migrated.

    WHY ``rows`` IS EXPLICIT AND NOT DEFAULTED. This is the exact hazard that
    got connection reuse deferred the first time, recorded at the time as:
    "callers set row_factory per call, so a shared cached connection leaks it
    between them (clinvar expects Row, biogrid expects tuples) = silently wrong
    data." It is a real hazard and it is real HERE: five sites in this module
    want ``sqlite3.Row`` and ``_gene_phenotype_strings_sync`` indexes its result
    as a TUPLE. On a shared connection whichever ran first would decide, and the
    loser would read the wrong column with no error.

    So the setting is not inherited, it is asserted on every acquisition. A
    keyword-only argument with no default means a new call site cannot forget to
    choose, and test_clinvar_connection_reuse enforces that every site here goes
    through this function.
    """
    conn = connect_ro(DB_PATH)
    if conn is not None:
        conn.row_factory = sqlite3.Row if rows else None
    return conn

DB_PATH = db_path("clinvar.db", "CLINVAR_DB_PATH")

MAX_RECORDS = 5


MAX_LANDSCAPE_POINTS = 6000

_AA_POSITION_RE = re.compile(r"p\.[A-Z][a-z]{2}(\d+)[A-Z][a-z]{2}")

_AA_CHANGE_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")
_AA_RESIDUE_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)(.*?)(?:\)|$)")

_HGVS_C_RE = re.compile(r"(c\.[^ )]+)")

_N_CONDITIONS_PLACEHOLDER_RE = re.compile(r"^\d+\s+conditions?$", re.IGNORECASE)


def _stars_for(review_status: str | None) -> int:
    """Map a ClinVar `review_status` string to the canonical 0–4 ★ scale.

    Reference: https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
      4 ★ — practice guideline
      3 ★ — reviewed by expert panel
      2 ★ — criteria provided, multiple submitters, no conflicts
      1 ★ — criteria provided, single submitter / conflicting interpretations
      0 ★ — no assertion criteria provided / no assertion provided
    """
    rs = (review_status or "").lower()
    if "practice guideline" in rs:
        return 4
    if "expert panel" in rs:
        return 3
    if "multiple submitters" in rs and "no conflict" in rs:
        return 2
    if "single submitter" in rs or "conflicting" in rs:
        return 1
    return 0


def _vcv(variation_id: int | None) -> str | None:
    """Build the canonical VCV accession from the integer variation_id —
    the variant_summary file doesn't carry the formatted accession."""
    if variation_id is None:
        return None
    return f"VCV{variation_id:09d}"


def _split_phenotypes(phenotype_list: str | None) -> list[str]:
    """variant_summary stores PhenotypeList as ``A|B|C`` (pipe-delimited).
    Drop empties and ClinVar's "not provided" / "not specified" stubs."""
    if not phenotype_list:
        return []
    out = []
    seen = set()
    for raw in phenotype_list.split("|"):
        p = raw.strip()
        if not p:
            continue
        low = p.lower()
        if low in {"not provided", "not specified"}:
            continue
        if _N_CONDITIONS_PLACEHOLDER_RE.match(p):
            continue
        if p in seen:
            continue
        seen.add(p)
        out.append(p)
    return out


def _hgvs_match(name: str, hgvs_c: str) -> bool:
    """Confirm the HGVS appears as a complete token in the ClinVar Name.

    ClinVar's ``Name`` looks like ``NM_000257.4(MYH7):c.1234C>T (p.Pro412Leu)``
    or ``GENE:c.886-2A>G``. The SQL LIKE pre-filter is a substring match,
    which catches false positives (e.g. ``c.123A>T`` matching
    ``c.1234A>T``). We tighten with a word-boundary check here.

    Matching is CASE-INSENSITIVE: HGVS nucleotide case is not semantic
    (``c.886g>a`` is the same variant as ``c.886G>A``), so a lowercased
    curator input must still match the uppercase ClinVar Name. This cannot
    create a false match — the nucleotide identity is still required — it only
    recovers a correct match that case would otherwise drop (which would
    silently suppress the exact-variant PP5/BP6 / PM5 surfaces).
    """
    if not name or not hgvs_c:
        return False
    pattern = r"(?:^|[^A-Za-z0-9])" + re.escape(hgvs_c) + r"(?:$|[^A-Za-z0-9])"
    return bool(re.search(pattern, name, re.IGNORECASE))


def _query_sync(
    gene: str,
    hgvs_c: str,
    hgvs_c_mane: str | None = None,
    chrom: str | None = None,
    pos: int | None = None,
    end: int | None = None,
) -> dict:
    """Synchronous DB lookup. Called from the async wrapper via to_thread.

    THREE MATCHING LAYERS, in decreasing confidence.

    1. ``hgvs_c`` — the curator's own coding change, matched as a whole token
       inside ClinVar's ``Name``. Allele-exact. Unchanged behaviour.
    2. ``hgvs_c_mane`` — the SAME variant's coding change on the MANE Select
       transcript, from VEP. Also allele-exact, and the layer that fixes the
       reported gap: ClinVar writes ``Name`` on its preferred (MANE) transcript,
       so a curator working on GATA4 ``NM_002052.5`` searched for ``c.886G>A``
       while ClinVar holds the same variant as ``NM_001308093.3:c.889G>A``
       (VCV000009030, Pathogenic/Likely pathogenic, 2★) and HeartVar reported
       it absent — silently removing PP5 from consideration.
    3. ``chrom``/``pos``/``end`` — other ClinVar records at the SAME NUCLEOTIDE.
       These are returned separately as ``same_position_records`` and are
       explicitly NOT counted as finding the variant, because the local mirror
       CANNOT confirm the allele: ``scripts/build_clinvar_db.py`` populates
       reference_allele/alternate_allele from variant_summary's
       ``ReferenceAllele``/``AlternateAllele`` columns, which are the literal
       string ``'na'`` in 4,456,618 of 4,456,656 rows (measured 2026-08-24).
       The real alleles live in the ``ReferenceAlleleVCF``/``AlternateAlleleVCF``
       columns, which the build does not read. Until that build changes, a
       position hit means "another submitted allele at this base", not "this
       variant". Each such record carries ``allele_confirmed: False`` so no
       downstream consumer can mistake it for an exact match.

    INDEL GUARD. Position matching is restricted to single-base loci
    (``pos == end`` on the query side, ``obj_type == 'single nucleotide
    variant'`` on the row side). The mirror stores ClinVar's ``Start``, which
    equals ``PositionVCF`` for 499,997 of 500,000 sampled GRCh38 SNV rows
    (99.9994%) but disagrees for 65.1% of sampled non-SNV rows — so matching
    indels on ``start`` would be wrong far more often than right. Indels fall
    through to layers 1-2 only.
    """
    if not DB_PATH.exists():
        return {
            "ok": False,
            "error": (
                f"ClinVar local DB not found at {DB_PATH}. "
                "Run `python3 scripts/build_clinvar_db.py` to build it."
            ),
        }

    _t_conn = perf_counter()
    conn = _conn(rows=True)
    if conn is None:
        return dict(_DB_MISSING)
    _d_conn = perf_counter() - _t_conn
    _d_exec = _d_fetch = 0.0
    try:
        cur = conn.cursor()
        _snv_pos = pos if (pos is not None and end is not None and pos == end) else None
        if hgvs_c or hgvs_c_mane or _snv_pos is not None:
            clauses, params = [], [gene]
            for tok in (hgvs_c, hgvs_c_mane):
                if tok:
                    clauses.append("name LIKE ?")
                    params.append(f"%{tok}%")
            if _snv_pos is not None:
                clauses.append(
                    "(start = ? AND stop = ? "
                    "AND obj_type = 'single nucleotide variant')"
                )
                params.extend([_snv_pos, _snv_pos])
            _t_exec = perf_counter()
            cur.execute(
                f"""
                SELECT variation_id, obj_type, name, gene_symbol,
                       clinical_significance, review_status,
                       number_submitters, phenotype_list, last_evaluated,
                       start, stop
                FROM variants
                WHERE gene_symbol = ?
                  AND ({" OR ".join(clauses)})
                ORDER BY number_submitters DESC, variation_id ASC
                LIMIT 100
                """,
                params,
            )
        else:
            _t_exec = perf_counter()
            cur.execute(
                """
                SELECT variation_id, obj_type, name, gene_symbol,
                       clinical_significance, review_status,
                       number_submitters, phenotype_list, last_evaluated,
                       start, stop
                FROM variants
                WHERE gene_symbol = ?
                ORDER BY number_submitters DESC, variation_id ASC
                LIMIT 100
                """,
                (gene,),
            )
        _d_exec = perf_counter() - _t_exec
        _t_fetch = perf_counter()
        rows = cur.fetchall()
        _d_fetch = perf_counter() - _t_fetch
    except sqlite3.Error as e:
        log.exception("ClinVar local DB query failed for %s %s", gene, hgvs_c)
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    finally:
        pass

    _t_proc = perf_counter()
    matched = []
    same_position: list[dict] = []
    total_submissions = 0
    all_conditions_seen: set[str] = set()
    all_conditions: list[str] = []
    _snv_pos = pos if (pos is not None and end is not None and pos == end) else None
    for row in rows:
        exact = False
        matched_via = None
        if hgvs_c and _hgvs_match(row["name"], hgvs_c):
            exact, matched_via = True, "hgvs_c"
        elif hgvs_c_mane and _hgvs_match(row["name"], hgvs_c_mane):
            exact, matched_via = True, "hgvs_c_mane"
        elif not hgvs_c and not hgvs_c_mane:
            exact, matched_via = True, "gene_only"
        if not exact:
            if _snv_pos is not None and row["start"] == _snv_pos:
                same_position.append({
                    "uid": str(row["variation_id"]) if row["variation_id"] is not None else "",
                    "accession": _vcv(row["variation_id"]),
                    "title": row["name"],
                    "variation_id": row["variation_id"],
                    "clinical_significance": row["clinical_significance"],
                    "review_status": row["review_status"],
                    "stars": _stars_for(row["review_status"]),
                    "last_evaluated": row["last_evaluated"],
                    "position": row["start"],
                    "consequence": _classify_consequence(row["name"]),
                    "allele_confirmed": False,
                })
            continue
        review_status = row["review_status"]
        n_sub = row["number_submitters"] or 1
        row_conditions = _split_phenotypes(row["phenotype_list"])
        total_submissions += n_sub
        for cond in row_conditions:
            if cond not in all_conditions_seen:
                all_conditions_seen.add(cond)
                all_conditions.append(cond)
        matched.append({
            "uid": str(row["variation_id"]) if row["variation_id"] is not None else "",
            "accession": _vcv(row["variation_id"]),
            "title": row["name"],
            "variation_id": row["variation_id"],
            "obj_type": row["obj_type"],
            "clinical_significance": row["clinical_significance"],
            "review_status": review_status,
            "stars": _stars_for(review_status),
            "last_evaluated": row["last_evaluated"],
            "conditions": row_conditions,
            "number_submitters": n_sub,
            "matched_via": matched_via,
        })

    display_records = matched[:MAX_RECORDS]
    log.info(
        "[timing] clinvar %s: conn=%.3fs exec=%.3fs fetch=%.3fs proc=%.3fs "
        "rows=%d matched=%d",
        gene, _d_conn, _d_exec, _d_fetch, perf_counter() - _t_proc,
        len(rows), len(matched),
    )
    return {
        "ok": True,
        "found": bool(matched),
        "uids": [r["uid"] for r in display_records],
        "records": display_records,
        "total_records": len(matched),
        "total_submissions": total_submissions,
        "phenotype_matched_submissions": total_submissions,
        "all_conditions": all_conditions,
        "same_position_records": sorted(
            same_position, key=lambda r: (-(r["stars"] or 0), str(r["title"])),
        )[:MAX_RECORDS],
        "same_position_count": len(same_position),
        "same_position_allele_confirmed": False,
        "matched_on_mane": any(
            r.get("matched_via") == "hgvs_c_mane" for r in matched
        ),
        "position_matching_used": _snv_pos is not None,
    }


async def fetch_clinvar(
    gene: str,
    hgvs_c: str,
    hgvs_c_mane: str | None = None,
    chrom: str | None = None,
    pos: int | None = None,
    end: int | None = None,
) -> dict:
    """Look up a variant in the local ClinVar SQLite cache.

    Returns the same dict shape as the previous E-utilities client:
    ``{"ok": bool, "found": bool, "uids": [...], "records": [...]}``, plus
    ``same_position_records`` / ``same_position_count`` when a single-base
    genomic position is supplied. SQLite calls are wrapped in
    ``run_local`` so the parallel DB gather doesn't block the loop.

    ``hgvs_c_mane`` is the same coding change expressed on the MANE Select
    transcript (from VEP). ClinVar writes its ``Name`` on its preferred —
    i.e. MANE Select — transcript, so passing this is what lets a curator
    working on a non-MANE transcript find the variant's own ClinVar record
    at all. All new arguments are optional and default to the previous
    behaviour exactly, so existing callers are unaffected.
    """
    gene = (gene or "").strip()
    hgvs_c = (hgvs_c or "").strip()
    hgvs_c_mane = (hgvs_c_mane or "").strip() or None
    if not gene:
        return {"ok": False, "error": "gene is required"}
    if hgvs_c_mane and hgvs_c_mane.lower() == hgvs_c.lower():
        hgvs_c_mane = None
    _t_dispatch = perf_counter()
    result = await run_local(
        _query_sync, gene, hgvs_c, hgvs_c_mane, chrom, pos, end
    )
    _d_total = perf_counter() - _t_dispatch
    if isinstance(result, dict):
        log.info("[timing] clinvar %s: run_local total=%.3fs "
                 "(pool queue = this minus the phases above)",
                 log_safe(gene), _d_total)
    return result


def _classify_tier(clinical_significance: str | None) -> str | None:
    """Map a ClinVar ``clinical_significance`` string to one of
    ``P / LP / VUS / LB / B``, or None when the entry is uninterpretable
    (``-``, ``not provided``, ``no classification from unflagged…``).

    Compound assertions ("Pathogenic/Likely pathogenic",
    "Benign/Likely benign") collapse into the stronger tier so the
    landscape totals don't double-count a single record. "Conflicting
    classifications of pathogenicity" maps to VUS — the curator should
    treat conflicts as an undetermined call, not a benign signal.
    """
    if not clinical_significance:
        return None
    cs = clinical_significance.strip().lower()
    if cs in ("", "-", "not provided", "no classifications from unflagged records"):
        return None
    if "conflicting" in cs:
        return "VUS"
    if cs.startswith("pathogenic"):
        return "P"
    if "likely pathogenic" in cs:
        return "LP"
    if cs.startswith("benign"):
        return "B"
    if "likely benign" in cs:
        return "LB"
    if "uncertain significance" in cs:
        return "VUS"
    return None


def _phenotype_matches(phenotype_list: str | None, keywords: list[str]) -> bool:
    """True if any keyword appears as a case-insensitive substring of
    ``phenotype_list``. Used by the landscape filter to narrow the
    variant set to disease-relevant submissions (e.g. ``["cardiomyopathy",
    "hypertrophic"]`` for an HCM proband)."""
    if not keywords:
        return True
    if not phenotype_list:
        return False
    haystack = phenotype_list.lower()
    return any(k.lower() in haystack for k in keywords if k)


_CSQ_P_RE = re.compile(r"p\.\(?([A-Za-z0-9*=]+)\)?")
_CSQ_C_RE = re.compile(r"c\.([^ )]+)")
_MISSENSE_P_RE = re.compile(r"^[A-Za-z]{3}\d+[A-Za-z]{3}$")


def _classify_consequence(name: str | None) -> str:
    """Coarse consequence bucket for a ClinVar Name: one of ``missense`` /
    ``truncating`` / ``splice`` / ``inframe`` / ``synonymous`` / ``utr`` /
    ``other``. Heuristic, from the HGVS protein/coding tokens only."""
    s = name or ""
    pm = _CSQ_P_RE.search(s)
    if pm:
        p = pm.group(1)
        low = p.lower()
        if "ter" in low or "*" in p or "fs" in low:
            return "truncating"
        if p.endswith("="):
            return "synonymous"
        if "del" in low or "dup" in low or "ins" in low:
            return "inframe"
        if _MISSENSE_P_RE.match(p):
            return "missense"
        return "other"
    cm = _CSQ_C_RE.search(s)
    if cm:
        c = cm.group(1)
        if c.startswith("-") or c.startswith("*"):
            return "utr"
        if re.search(r"\d[+-]\d", c):
            return "splice"
    return "other"


def _gene_variant_landscape_sync(
    gene: str, condition_keywords: list[str] | None
) -> dict:
    """Synchronous DB pull for ``get_gene_variant_landscape``."""
    if not DB_PATH.exists():
        return {
            "ok": False,
            "error": (
                f"ClinVar local DB not found at {DB_PATH}. "
                "Run `python3 scripts/build_clinvar_db.py` to build it."
            ),
        }

    conn = _conn(rows=True)
    if conn is None:
        return dict(_DB_MISSING)
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT variation_id, name, clinical_significance, review_status,
                   number_submitters, phenotype_list, start
            FROM variants
            WHERE gene_symbol = ?
              AND (
                   clinical_significance LIKE '%pathogenic%'
                OR clinical_significance LIKE '%benign%'
                OR clinical_significance LIKE '%uncertain significance%'
              )
            ORDER BY number_submitters DESC, variation_id ASC
            """,
            (gene,),
        )
        rows = cur.fetchall()
    except sqlite3.Error as e:
        log.exception("ClinVar landscape query failed for %s", gene)
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    finally:
        pass

    keywords = [k for k in (condition_keywords or []) if k and k.strip()]
    tier_counts = {"P": 0, "LP": 0, "VUS": 0, "LB": 0, "B": 0}
    tier_counts_by_csq: dict[str, dict[str, int]] = {}
    has_two_star_plp = False
    filtered_total = 0
    skipped_by_keywords = 0
    # Two-pass tally so each P/LP record is attributed to exactly ONE
    plp_record_phens: list[list[str]] = []
    plp_record_display: list[dict[str, str]] = []
    plp_phen_overall: dict[str, int] = {}
    pos_agg: dict[tuple[int, str, str], dict] = {}
    for row in rows:
        if keywords and not _phenotype_matches(row["phenotype_list"], keywords):
            skipped_by_keywords += 1
            continue
        tier = _classify_tier(row["clinical_significance"])
        if tier is None:
            continue
        tier_counts[tier] += 1
        filtered_total += 1
        _csq = _classify_consequence(row["name"])
        tier_counts_by_csq.setdefault(
            _csq, {"P": 0, "LP": 0, "VUS": 0, "LB": 0, "B": 0}
        )[tier] += 1
        _gpos = row["start"]
        if _gpos is not None:
            _stars = _stars_for(row["review_status"])
            _aa_m = _AA_POSITION_RE.search(row["name"] or "")
            _aa = int(_aa_m.group(1)) if _aa_m else None
            _key = (_gpos, tier, _csq)
            _agg = pos_agg.get(_key)
            if _agg is None:
                pos_agg[_key] = {
                    "gpos": _gpos, "aa": _aa, "tier": tier, "csq": _csq,
                    "stars": _stars, "count": 1,
                }
            else:
                _agg["count"] += 1
                if _stars > _agg["stars"]:
                    _agg["stars"] = _stars
                if _agg["aa"] is None and _aa is not None:
                    _agg["aa"] = _aa
        if tier in ("P", "LP"):
            if _stars_for(row["review_status"]) >= 2:
                has_two_star_plp = True
            record_keys: list[str] = []
            record_display: dict[str, str] = {}
            seen_in_row: set[str] = set()
            for part in _PHENOTYPE_DELIM_RE.split(row["phenotype_list"] or ""):
                cleaned = part.strip()
                key = cleaned.lower()
                if (
                    not key
                    or key in _NON_INFORMATIVE_PHENOTYPES
                    or _N_CONDITIONS_PLACEHOLDER_RE.match(key)
                    or key in seen_in_row
                ):
                    continue
                seen_in_row.add(key)
                record_keys.append(key)
                record_display[key] = cleaned
                plp_phen_overall[key] = plp_phen_overall.get(key, 0) + 1
            plp_record_phens.append(record_keys)
            plp_record_display.append(record_display)

    plp_phen_counts: dict[str, int] = {}
    plp_phen_display: dict[str, str] = {}
    for record_keys, record_display in zip(
        plp_record_phens, plp_record_display, strict=True
    ):
        if not record_keys:
            continue
        primary = record_keys[0]
        for k in record_keys[1:]:
            if plp_phen_overall[k] > plp_phen_overall[primary]:
                primary = k
        plp_phen_counts[primary] = plp_phen_counts.get(primary, 0) + 1
        if primary not in plp_phen_display:
            plp_phen_display[primary] = record_display[primary]

    plp = tier_counts["P"] + tier_counts["LP"]
    plp_fraction_pct = (
        round(100.0 * plp / filtered_total, 1) if filtered_total else None
    )
    top_phen = sorted(
        plp_phen_counts.items(), key=lambda kv: (-kv[1], plp_phen_display[kv[0]])
    )[:4]
    plp_top_phenotypes = [
        {"phenotype": plp_phen_display[k], "count": c} for k, c in top_phen
    ]

    positions = list(pos_agg.values())
    positions_truncated = False
    if len(positions) > MAX_LANDSCAPE_POINTS:
        _tier_pri = {"P": 0, "LP": 1, "B": 2, "LB": 3, "VUS": 4}
        positions.sort(
            key=lambda d: (_tier_pri.get(d["tier"], 9), -d["stars"], -d["count"])
        )
        positions = positions[:MAX_LANDSCAPE_POINTS]
        positions_truncated = True
    positions.sort(key=lambda d: (d["gpos"], d["tier"]))
    g_min = min((p["gpos"] for p in positions), default=None)
    g_max = max((p["gpos"] for p in positions), default=None)

    return {
        "ok": True,
        "gene": gene,
        "total_classified": filtered_total,
        "tier_counts": tier_counts,
        "tier_counts_by_csq": tier_counts_by_csq,
        "positions": positions,
        "positions_truncated": positions_truncated,
        "g_min": g_min,
        "g_max": g_max,
        "plp_fraction_pct": plp_fraction_pct,
        "has_two_star_plp": has_two_star_plp,
        "condition_keywords": keywords or None,
        "records_examined": len(rows),
        "records_filtered_out": skipped_by_keywords,
        "plp_top_phenotypes": plp_top_phenotypes,
    }


_PHENOTYPE_DELIM_RE = re.compile(r"[|;]+")

_NON_INFORMATIVE_PHENOTYPES: frozenset[str] = frozenset({
    "", "-", "not provided", "not specified",
    "no classifications from unflagged records",
})


def _gene_phenotype_strings_sync(gene: str) -> list[str]:
    """Synchronous DB pull for ``get_gene_phenotype_strings``."""
    if not DB_PATH.exists():
        return []
    conn = _conn(rows=False)
    if conn is None:
        return []
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT DISTINCT phenotype_list
            FROM variants
            WHERE gene_symbol = ?
              AND phenotype_list IS NOT NULL
              AND phenotype_list != ''
            """,
            (gene,),
        )
        rows = cur.fetchall()
    except sqlite3.Error:
        log.exception("ClinVar phenotype-strings query failed for %s", gene)
        return []
    finally:
        pass

    seen: set[str] = set()
    out: list[str] = []
    for (raw,) in rows:
        for part in _PHENOTYPE_DELIM_RE.split(raw or ""):
            p = part.strip().lower()
            if (
                p in _NON_INFORMATIVE_PHENOTYPES
                or _N_CONDITIONS_PLACEHOLDER_RE.match(p)
                or p in seen
            ):
                continue
            seen.add(p)
            out.append(p)
    return out


async def get_gene_phenotype_strings(gene_symbol: str) -> list[str]:
    """Return every distinct ClinVar PhenotypeList entry for ``gene_symbol``
    as a lowercased, deduplicated list.

    PhenotypeList values are split on both ``|`` (ClinVar's group
    separator) and ``;`` (within-group separator) so the caller gets
    individual disease titles rather than compound multi-disease
    strings. Stubs like "not provided" / "not specified" are filtered
    out — they don't carry phenotype information, and including them
    as keywords would broaden the substring filter to every variant
    with an empty phenotype.

    Fast local SQLite query — no network call. Used by the landscape
    keyword-builder in app.py to derive the filter directly from the
    gene's actual ClinVar phenotype vocabulary, which avoids the
    earlier failure mode of trying to anticipate every possible
    phenotype string via a curated pattern list.
    """
    gene = (gene_symbol or "").strip()
    if not gene:
        return []
    return await run_local(_gene_phenotype_strings_sync, gene)


async def get_gene_variant_landscape(
    gene_symbol: str, condition_keywords: list[str] | None = None
) -> dict:
    """Summarise the ClinVar variant landscape for ``gene_symbol``.

    Returns counts per ACMG tier (P / LP / VUS / LB / B), the P+LP
    fraction as a percentage of classified records, and a flag for any
    P/LP records with ≥2-star review status. When ``condition_keywords``
    is provided, only records whose ``PhenotypeList`` mentions one of the
    keywords (case-insensitive) are counted — useful for narrowing to
    disease-specific submissions.

    Implementation note: the task description called for an eUtils
    esearch query, but this project has migrated to a local SQLite cache
    of ``variant_summary.txt.gz`` (see ``backend/clients/clinvar.py``
    header). The local pull is faster, has no rate limits, and surfaces
    the same fields. Rebuild the DB monthly with
    ``python3 scripts/build_clinvar_db.py``.
    """
    gene = (gene_symbol or "").strip()
    if not gene:
        return {"ok": False, "error": "gene_symbol is required"}
    return await run_local(
        _gene_variant_landscape_sync, gene, condition_keywords
    )


def _normalize_hgvs_c(hgvs_c: str | None) -> str | None:
    """Lowercase + strip whitespace from an HGVS-c token for equality
    comparison. ClinVar stores ``c.1785T>G``; the proband token from the
    curate request may arrive with surrounding whitespace or a transcript
    prefix (``NM_004333.6:c.1785T>G``). Reduce both to the bare ``c.…``
    coding change so the proband's OWN ClinVar record is matched exactly.
    Returns None when no ``c.`` token is present."""
    if not hgvs_c:
        return None
    m = _HGVS_C_RE.search(hgvs_c)
    if not m:
        return None
    return m.group(1).strip().lower()


def _pm5_evidence_sync(
    gene: str,
    protein_position: int,
    proband_alt_aa: str | None = None,
    proband_hgvs_c: str | None = None,
    mane_protein_position: int | None = None,
    proband_hgvs_c_mane: str | None = None,
) -> dict:
    """Synchronous DB pull for ``get_pm5_evidence``.

    Scans the gene's P/LP variants and keeps those whose protein-change
    token sits at exactly ``protein_position``. Skips records without a
    parseable ``p.XxxNNNYyy`` token (synonymous / non-coding / del-ins
    notations that don't follow the 3-letter format).

    Proband self-exclusion (added to break ClinVar self-confirmation
    circularity — eRepo VCEP variants now sit in ClinVar at ≥2★, so the
    proband's own record would otherwise satisfy its own PM5):
      - ``proband_alt_aa``: 3-letter alt AA of the proband (e.g. "Leu").
        Candidates encoding the SAME amino-acid substitution as the
        proband are dropped from the PM5 ``candidates`` list — PM5 is a
        *different* missense change at the residue, never the same one.
        Those same-AA records are surfaced separately as ``ps1_candidates``
        (ACMG-2015 PS1 = same amino-acid change, different nucleotide),
        after removing the proband's own record by ``proband_hgvs_c``.
      - ``proband_hgvs_c``: the proband's coding change (``c.1785T>G``).
        Used only to drop the proband's OWN record from the PS1 surface so
        a variant can't confirm itself; a DIFFERENT nucleotide encoding the
        same AA is a valid PS1 comparison variant and is kept.
      - ``proband_hgvs_c_mane``: the SAME coding change written on the MANE
        Select transcript. Required, not optional-in-practice, for the same
        reason the residue number is: ClinVar's Name carries the MANE coding
        change, so comparing it against the curator's token alone makes the
        proband's OWN record look like a different nucleotide. GATA4
        NM_002052.5:c.886G>A is ClinVar's c.889G>A — self-exclusion by
        "c.886g>a" != "c.889g>a" would let the proband's own 2★
        Pathogenic/Likely-pathogenic record satisfy its own PS1 at Strong
        (+4). Both tokens are checked; a match on either means "self".

    The task also mentioned a fallback by genomic region ±3 bp around
    the variant's codon — we omit that here because the method signature
    doesn't carry genomic coordinates, and PM5 is defined at the residue
    level ("different missense change at the same amino acid residue").
    Protein-position match is the canonical PM5 query.

    TRANSCRIPT ANCHORING (``mane_protein_position``). ClinVar's ``Name`` is
    written on ClinVar's preferred transcript, which is the MANE Select one.
    ``protein_position`` comes from VEP's hgvsp on the transcript the CURATOR
    supplied, which is frequently NOT MANE. Where the two isoforms are
    numbered differently the raw integer compare below is between two
    DIFFERENT residues — it both missed real same-residue evidence and
    matched unrelated residues (GATA4: NM_002052.5 residue 303 is the same
    residue as MANE NM_001308093.3 residue 304, and a query at 296 surfaced
    the pathogenic p.Cys296Ser, which is residue 295 of NM_002052.5, as a
    PS1 comparator). So we MATCH on ``mane_protein_position`` when it is
    available and keep ``protein_position`` for DISPLAY. When no MANE
    annotation exists we fall back to ``protein_position`` — i.e. unchanged
    pre-fix behaviour — rather than guessing an offset.
    """
    if not DB_PATH.exists():
        return {
            "ok": False,
            "error": (
                f"ClinVar local DB not found at {DB_PATH}. "
                "Run `python3 scripts/build_clinvar_db.py` to build it."
            ),
        }

    conn = _conn(rows=True)
    if conn is None:
        return dict(_DB_MISSING)
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT variation_id, name, clinical_significance, review_status,
                   number_submitters, phenotype_list
            FROM variants
            WHERE gene_symbol = ?
              AND name LIKE '%p.%'
              AND (
                   clinical_significance LIKE '%pathogenic%'
                OR clinical_significance LIKE '%likely pathogenic%'
              )
              AND clinical_significance NOT LIKE '%benign%'
              AND clinical_significance NOT LIKE '%conflicting%'
            ORDER BY number_submitters DESC, variation_id ASC
            """,
            (gene,),
        )
        rows = cur.fetchall()
    except sqlite3.Error as e:
        log.exception("ClinVar PM5 query failed for %s pos=%s", gene, protein_position)
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    finally:
        pass

    candidates: list[dict] = []
    ps1_candidates: list[dict] = []
    non_missense: list[dict] = []
    self_excluded: list[dict] = []
    match_position = (
        mane_protein_position if mane_protein_position is not None
        else protein_position
    )
    target = str(match_position)
    proband_alt = (proband_alt_aa or "").strip() or None
    proband_cs = {
        c for c in (
            _normalize_hgvs_c(proband_hgvs_c),
            _normalize_hgvs_c(proband_hgvs_c_mane),
        ) if c
    }
    for row in rows:
        m = _AA_CHANGE_RE.search(row["name"] or "")
        if not m or m.group(2) != target:
            continue
        tier = _classify_tier(row["clinical_significance"])
        if tier not in ("P", "LP"):
            continue
        candidate_alt = m.group(3)
        rec = {
            "variation_id": row["variation_id"],
            "accession": _vcv(row["variation_id"]),
            "name": row["name"],
            "clinical_significance": row["clinical_significance"],
            "tier": tier,
            "review_status": row["review_status"],
            "stars": _stars_for(row["review_status"]),
        }
        if (
            _classify_consequence(row["name"]) != "missense"
            or candidate_alt == m.group(1)
        ):
            non_missense.append(rec)
            continue
        is_same_aa = proband_alt is not None and candidate_alt == proband_alt
        if is_same_aa:
            cand_c = _normalize_hgvs_c(row["name"])
            is_self = (not proband_cs) or (cand_c in proband_cs)
            if not is_self:
                ps1_candidates.append(rec)
            elif proband_cs:
                self_excluded.append(rec)
            continue
        candidates.append(rec)

    two_star = [c for c in candidates if (c.get("stars") or 0) >= 2]
    ps1_two_star = [c for c in ps1_candidates if (c.get("stars") or 0) >= 2]
    return {
        "ok": True,
        "gene": gene,
        "protein_position": protein_position,
        "matched_protein_position": match_position,
        "numbering_differs": match_position != protein_position,
        "proband_alt_aa": proband_alt,
        "candidates": candidates,
        "count": len(candidates),
        "count_two_star": len(two_star),
        "non_missense_excluded": non_missense,
        "non_missense_excluded_count": len(non_missense),
        "non_missense_excluded_two_star_count": len(
            [c for c in non_missense if (c.get("stars") or 0) >= 2]
        ),
        "ps1_candidates": ps1_candidates,
        "ps1_count": len(ps1_candidates),
        "ps1_count_two_star": len(ps1_two_star),
        "ps1_self_excluded": self_excluded,
        "ps1_self_excluded_count": len(self_excluded),
    }


def _same_residue_sync(
    gene: str,
    match_position: int,
    display_position: int,
    proband_alt_aa: str | None,
    proband_hgvs_c: str | None,
    proband_hgvs_c_mane: str | None = None,
) -> dict:
    """EVERY classified ClinVar record at one residue, each labelled with
    whether it can support PM5 and, when it cannot, why not.

    Distinct from ``_pm5_evidence_sync``, which returns only the P/LP
    missense records that actually drive the criterion. This is the curator-
    facing context surface: it deliberately includes VUS, benign and
    SYNONYMOUS records, because the whole point is to show what has been
    submitted at the residue without letting any of it be mistaken for PM5
    support. The prompt carries a hard "MUST NOT infer PS1 or PM5 from
    nearby variants" rule; every row here therefore carries an explicit
    ``pm5_eligible`` flag and a reason, so the surface cannot undermine it.

    PM5 (Richards 2015) needs a DIFFERENT MISSENSE change at the residue that
    is itself established Pathogenic/Likely pathogenic. So a record is
    PM5-ineligible when it is synonymous or truncating (not a missense
    change), when it encodes the SAME amino-acid substitution as the proband
    (that is PS1 territory, or the proband's own record), or when it is not
    classified P/LP. The GATA4 report is exactly the case this guards:
    ``p.Gly304=`` is Gly→Gly and cannot support PM5 under any numbering,
    while ``p.Gly304Arg`` is a different missense at the residue but is a VUS,
    not an established pathogenic comparator.

    ``match_position`` is the MANE residue number (ClinVar's frame);
    ``display_position`` is the curator's. Both are echoed back so the UI can
    label the panel with both.
    """
    if not DB_PATH.exists():
        return {"ok": False, "error": f"ClinVar local DB not found at {DB_PATH}."}
    conn = _conn(rows=True)
    if conn is None:
        return dict(_DB_MISSING)
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT variation_id, name, clinical_significance, review_status,
                   number_submitters, phenotype_list, start
            FROM variants
            WHERE gene_symbol = ?
              AND name LIKE '%p.%'
            ORDER BY number_submitters DESC, variation_id ASC
            """,
            (gene,),
        )
        rows = cur.fetchall()
    except sqlite3.Error as e:
        log.exception("ClinVar same-residue query failed for %s %s", gene, match_position)
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    finally:
        pass

    target = str(match_position)
    proband_alt = (proband_alt_aa or "").strip() or None
    proband_cs = {
        c for c in (
            _normalize_hgvs_c(proband_hgvs_c),
            _normalize_hgvs_c(proband_hgvs_c_mane),
        ) if c
    }
    records: list[dict] = []
    for row in rows:
        m = _AA_RESIDUE_RE.search(row["name"] or "")
        if not m or m.group(2) != target:
            continue
        tier = _classify_tier(row["clinical_significance"])
        csq = _classify_consequence(row["name"])
        alt_m = _AA_CHANGE_RE.search(row["name"] or "")
        alt_aa = alt_m.group(3) if alt_m else None
        is_self = bool(
            proband_cs and _normalize_hgvs_c(row["name"]) in proband_cs
        )
        same_aa = bool(proband_alt and alt_aa and alt_aa == proband_alt)
        if csq == "synonymous":
            eligible, reason = False, (
                "synonymous (no amino-acid change) — PM5 requires a DIFFERENT "
                "MISSENSE change at the residue"
            )
        elif csq != "missense":
            eligible, reason = False, (
                f"{csq} — PM5 requires a missense change at the residue"
            )
        elif is_self:
            eligible, reason = False, "the proband's own record — circular"
        elif same_aa:
            eligible, reason = False, (
                "same amino-acid change as the proband — PS1 territory, not PM5"
            )
        elif tier not in ("P", "LP"):
            eligible, reason = False, (
                f"classified {row['clinical_significance']!r} — PM5 requires an "
                "established Pathogenic/Likely-pathogenic comparator"
            )
        elif _stars_for(row["review_status"]) < 2:
            eligible, reason = False, (
                "below the ≥2★ review bar HeartVar requires for PM5"
            )
        else:
            eligible, reason = True, None
        records.append({
            "variation_id": row["variation_id"],
            "accession": _vcv(row["variation_id"]),
            "name": row["name"],
            "clinical_significance": row["clinical_significance"],
            "tier": tier,
            "review_status": row["review_status"],
            "stars": _stars_for(row["review_status"]),
            "consequence": csq,
            "alt_aa": alt_aa,
            "position": row["start"],
            "is_proband_own_record": is_self,
            "pm5_eligible": eligible,
            "pm5_ineligible_reason": reason,
        })
    records.sort(key=lambda r: (not r["pm5_eligible"], -(r["stars"] or 0), r["name"]))
    return {
        "ok": True,
        "gene": gene,
        "display_position": display_position,
        "matched_position": match_position,
        "numbering_differs": match_position != display_position,
        "records": records,
        "count": len(records),
        "pm5_eligible_count": sum(1 for r in records if r["pm5_eligible"]),
    }


async def get_same_residue_records(
    gene_symbol: str,
    display_position: int | None,
    match_position: int | None = None,
    proband_alt_aa: str | None = None,
    proband_hgvs_c: str | None = None,
    proband_hgvs_c_mane: str | None = None,
) -> dict:
    """Async wrapper for :func:`_same_residue_sync`. ``match_position``
    defaults to ``display_position`` when no MANE annotation was available."""
    gene = (gene_symbol or "").strip()
    if not gene:
        return {"ok": False, "error": "gene_symbol is required"}
    if display_position is None:
        return {"ok": True, "not_applicable": True, "records": [], "count": 0,
                "reason": "no protein residue for this variant"}
    try:
        disp = int(display_position)
        match = int(match_position) if match_position is not None else disp
    except (TypeError, ValueError):
        return {"ok": False, "error": "positions must be int-coercible"}
    return await run_local(
        _same_residue_sync, gene, match, disp, proband_alt_aa, proband_hgvs_c,
        proband_hgvs_c_mane,
    )


def _domain_plp_sync(
    gene: str,
    aa_start: int,
    aa_end: int,
    exclude_position: int | None,
) -> dict:
    """Synchronous DB pull for ``get_domain_plp_evidence``.

    Same protein-change filter as the PM5 query (``name LIKE '%p.%'``
    plus the (likely)pathogenic significance filter); the residue
    selector is applied in Python via ``_AA_POSITION_RE`` because the
    SQLite column doesn't carry a structured aa_position field.
    """
    if not DB_PATH.exists():
        return {
            "ok": False,
            "error": (
                f"ClinVar local DB not found at {DB_PATH}. "
                "Run `python3 scripts/build_clinvar_db.py` to build it."
            ),
        }

    conn = _conn(rows=True)
    if conn is None:
        return dict(_DB_MISSING)
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT variation_id, name, clinical_significance, review_status,
                   number_submitters, phenotype_list
            FROM variants
            WHERE gene_symbol = ?
              AND name LIKE '%p.%'
              AND (
                   clinical_significance LIKE '%pathogenic%'
                OR clinical_significance LIKE '%likely pathogenic%'
              )
              AND clinical_significance NOT LIKE '%benign%'
              AND clinical_significance NOT LIKE '%conflicting%'
            ORDER BY number_submitters DESC, variation_id ASC
            """,
            (gene,),
        )
        rows = cur.fetchall()
    except sqlite3.Error as e:
        log.exception(
            "ClinVar domain P/LP query failed for %s %s-%s", gene, aa_start, aa_end,
        )
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    finally:
        pass

    hits: list[dict] = []
    count_p = 0
    count_lp = 0
    has_two_star_plp = False
    exclude_str = str(exclude_position) if exclude_position is not None else None
    for row in rows:
        m = _AA_POSITION_RE.search(row["name"] or "")
        if not m:
            continue
        pos = int(m.group(1))
        if not (aa_start <= pos <= aa_end):
            continue
        if exclude_str is not None and m.group(1) == exclude_str:
            continue
        tier = _classify_tier(row["clinical_significance"])
        if tier not in ("P", "LP"):
            continue
        stars = _stars_for(row["review_status"])
        if stars >= 2:
            has_two_star_plp = True
        if tier == "P":
            count_p += 1
        else:
            count_lp += 1
        hits.append({
            "variation_id": row["variation_id"],
            "accession": _vcv(row["variation_id"]),
            "name": row["name"],
            "clinical_significance": row["clinical_significance"],
            "tier": tier,
            "review_status": row["review_status"],
            "stars": stars,
            "number_submitters": row["number_submitters"],
            "aa_position": pos,
            "conditions": _split_phenotypes(row["phenotype_list"]),
        })

    return {
        "ok": True,
        "gene": gene,
        "domain_start": aa_start,
        "domain_end": aa_end,
        "exclude_position": exclude_position,
        "total_plp": count_p + count_lp,
        "count_P": count_p,
        "count_LP": count_lp,
        "has_two_star_plp": has_two_star_plp,
        "all_hits": hits,
    }


async def get_domain_plp_evidence(
    gene_symbol: str,
    aa_start: int | None,
    aa_end: int | None,
    exclude_position: int | None = None,
    numbering_offset: int = 0,
) -> dict:
    """Count P/LP ClinVar variants within ``[aa_start, aa_end]`` of
    ``gene_symbol`` — the empirical basis for ACMG PM1 ("located in a
    mutational hotspot or critical and well-established functional
    domain").

    ``exclude_position`` lets the caller drop the curated variant's own
    residue from the count so a variant doesn't count itself when its
    residue already carries a P/LP submission.

    ``numbering_offset`` translates ``aa_start``/``aa_end``/
    ``exclude_position`` — which arrive in the numbering of the transcript
    the CURATOR supplied (the domain boundaries come from UniProt and are
    selected with that same position) — into the numbering of ClinVar's
    ``Name`` field, which is written on ClinVar's preferred (MANE Select)
    transcript. It is ``mane_residue - supplied_residue`` for the queried
    variant. Without it, a GATA4 query on NM_002052.5 excluded residue 303
    while ClinVar names the proband's own residue 304, so the proband could
    count itself toward its own PM1, and the window was off by one at both
    ends. ASSUMPTION (documented, not verified per gene): the offset is
    constant across the domain, which holds when the isoform difference lies
    outside the domain and fails when it lies inside it. Defaults to 0 —
    unchanged pre-fix behaviour — whenever no MANE annotation is available.

    Returns ``{"ok": True, "not_applicable": True, "reason": ...}`` when
    the domain boundaries aren't available (variant outside any
    annotated domain, or UniProt lookup failed). ``top_variants`` is
    capped at 5 (sorted by stars DESC then submitter count DESC); the
    aggregate counts cover every P/LP hit in the range.
    """
    gene = (gene_symbol or "").strip()
    if not gene:
        return {"ok": False, "error": "gene_symbol is required"}
    if aa_start is None or aa_end is None:
        return {
            "ok": True,
            "not_applicable": True,
            "reason": "domain boundaries not available",
        }
    try:
        s = int(aa_start)
        e = int(aa_end)
    except (TypeError, ValueError):
        return {
            "ok": False,
            "error": (
                f"aa_start/aa_end must be int-coercible, got "
                f"{aa_start!r}/{aa_end!r}"
            ),
        }
    if s > e:
        s, e = e, s
    excl: int | None = None
    if exclude_position is not None:
        try:
            excl = int(exclude_position)
        except (TypeError, ValueError):
            excl = None
    try:
        offset = int(numbering_offset or 0)
    except (TypeError, ValueError):
        offset = 0
    q_start, q_end = s + offset, e + offset
    q_excl = (excl + offset) if excl is not None else None
    result = await run_local(
        _domain_plp_sync, gene, q_start, q_end, q_excl
    )
    if not result.get("ok"):
        return result
    if offset:
        result["numbering_offset"] = offset
        result["matched_aa_start"] = q_start
        result["matched_aa_end"] = q_end
    hits = result.pop("all_hits", [])
    top = sorted(
        hits, key=lambda h: (-(h.get("stars") or 0), -(h.get("number_submitters") or 0)),
    )[:5]
    result["top_variants"] = top
    return result


async def get_pm5_evidence(
    gene_symbol: str,
    protein_position: int | None,
    proband_alt_aa: str | None = None,
    proband_hgvs_c: str | None = None,
    mane_protein_position: int | None = None,
    proband_hgvs_c_mane: str | None = None,
) -> dict:
    """Find P/LP missense variants at ``protein_position`` in
    ``gene_symbol`` — the canonical PM5 query ("different missense change
    at the same amino acid residue as a known pathogenic missense
    variant"). Returns ``{"ok": True, "candidates": [...]}`` with HGVS
    notation and star rating for each match; empty ``candidates`` list
    when no match is found.

    ``proband_alt_aa`` (3-letter, e.g. "Leu") and ``proband_hgvs_c``
    (``c.1785T>G``) drive proband self-exclusion: candidates encoding the
    SAME amino-acid substitution as the proband are removed from the PM5
    list (PM5 is a *different* change) and instead surfaced under
    ``ps1_candidates`` (same-AA = ACMG-2015 PS1), with the proband's own
    record dropped by coding change. When the alt AA isn't supplied, the
    same-residue behaviour is unchanged and ``ps1_candidates`` is empty.

    ``proband_hgvs_c_mane`` is the proband's own coding change on the MANE
    Select transcript — needed so the PS1 self-exclusion can recognise the
    proband's own ClinVar record, which is named in MANE coordinates.

    ``mane_protein_position`` is the same residue expressed on the MANE
    Select transcript. ClinVar Names are written on ClinVar's preferred
    (MANE Select) transcript, so that is the number the retrieval matches
    on; ``protein_position`` (the curator's transcript) is carried through
    for display. Pass None when the variant has no MANE annotation and the
    supplied-transcript number is used for both.

    The caller is responsible for ensuring the proband variant is a
    missense — when it isn't, pass ``protein_position=None`` and the
    method short-circuits with ``not_applicable``.
    """
    gene = (gene_symbol or "").strip()
    if not gene:
        return {"ok": False, "error": "gene_symbol is required"}
    if protein_position is None:
        return {
            "ok": True,
            "gene": gene,
            "not_applicable": True,
            "reason": "protein_position not supplied (variant is not missense or VEP failed)",
            "candidates": [],
            "count": 0,
            "ps1_candidates": [],
            "ps1_count": 0,
            "ps1_count_two_star": 0,
            "splice_ps1": await get_splice_ps1_evidence(gene, proband_hgvs_c_mane
                                                        or proband_hgvs_c),
        }
    try:
        position = int(protein_position)
    except (TypeError, ValueError):
        return {
            "ok": False,
            "error": f"protein_position must be int-coercible, got {protein_position!r}",
        }
    try:
        mane_position = (
            int(mane_protein_position) if mane_protein_position is not None
            else None
        )
    except (TypeError, ValueError):
        mane_position = None
    return await run_local(
        _pm5_evidence_sync, gene, position, proband_alt_aa, proband_hgvs_c,
        mane_position, proband_hgvs_c_mane,
    )


_CANONICAL_SPLICE_RE = re.compile(
    r"c\.(\d+)([+-])([12])(?![0-9])(.*)$", re.IGNORECASE
)


def parse_canonical_splice_offset(hgvs_c: str | None) -> tuple[str, str] | None:
    """``('732+1', 'G>A')`` for a canonical +/-1,2 splice variant, else None.

    The returned position string is what ClinVar's own ``name`` column spells,
    so it can be matched textually without a coordinate lift-over. The second
    element is the allele change, used only to exclude the proband's own record.
    """
    m = _CANONICAL_SPLICE_RE.search(str(hgvs_c or ""))
    if not m:
        return None
    exon_end, sign, offset, allele = m.groups()
    return f"{exon_end}{sign}{offset}", (allele or "").strip()


def _splice_ps1_sync(gene: str, position: str, proband_allele: str) -> dict:
    """P/LP ClinVar records at the same canonical splice position, >=2*.

    Mirrors the missense PS1 bar: established Pathogenic/Likely pathogenic at
    >=2 stars, reached via a DIFFERENT nucleotide change, with the proband's own
    record excluded.
    """
    if not DB_PATH.exists():
        return {"ok": False, "error": f"ClinVar local DB not found at {DB_PATH}."}
    conn = _conn(rows=True)
    try:
        rows = conn.execute(
            "SELECT name, clinical_significance, review_status FROM variants "
            "WHERE gene_symbol = ? AND name LIKE ? "
            "AND clinical_significance LIKE '%athogenic%'",
            (gene.upper(), f"%c.{position}%"),
        ).fetchall()
    except sqlite3.Error as exc:
        return {"ok": False, "error": f"ClinVar query failed: {exc}"}
    finally:
        conn.close()

    cands, self_excluded = [], []
    for row in rows:
        sig = (row["clinical_significance"] or "").lower()
        if "conflicting" in sig or "benign" in sig:
            continue
        parsed = parse_canonical_splice_offset(row["name"])
        if not parsed or parsed[0] != position:
            continue
        rec = {
            "name": row["name"],
            "significance": row["clinical_significance"],
            "stars": _stars_for(row["review_status"]),
        }
        if proband_allele and parsed[1].upper() == proband_allele.upper():
            self_excluded.append(rec)
            continue
        if rec["stars"] >= 2:
            cands.append(rec)
    cands.sort(key=lambda r: (-(r["stars"] or 0), r["name"]))
    return {
        "ok": True,
        "position": position,
        "candidates": cands,
        "count": len(cands),
        "self_excluded": self_excluded,
    }


async def get_splice_ps1_evidence(
    gene_symbol: str, hgvs_c: str | None
) -> dict:
    """PS1-by-splicing-similarity evidence for a canonical splice variant."""
    gene = (gene_symbol or "").strip()
    parsed = parse_canonical_splice_offset(hgvs_c)
    if not gene or not parsed:
        return {
            "ok": True, "not_applicable": True, "candidates": [], "count": 0,
            "reason": (
                "not a canonical +/-1,2 splice variant — the Walker 2023 PS1 "
                "splice rule is applied only where two changes at the same "
                "position abolish the same site by construction"
            ),
        }
    return await run_local(_splice_ps1_sync, gene, parsed[0], parsed[1])
