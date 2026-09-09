"""Local-first gnomAD lookup — gene constraint + variant frequency.

HeartVar fans out ~16 live lookups per ``/api/curate``; this client moves the
two gnomAD pulls off the live GraphQL API and onto local caches built at
deploy time, keeping the live API as a graceful fallback so we stop
overloading / being blacklisted by the Broad's gnomAD endpoint.

Two local sources, both built by ``scripts/build_gnomad_*_db.py``:

  * gene CONSTRAINT — ``data/gnomad_constraint.db`` (one row per gene, the
    MANE Select / canonical transcript's pLI / LOEUF / oe_mis / mis_z / syn_z).
    Served in place of the live ``gene`` query. Tens of MB; built for real.

  * variant FREQUENCY — ``data/gnomad_freq.db`` (one row per panel variant,
    the EXACT ``variant`` blob the GraphQL ``variant`` query returns, stored as
    a JSON payload). Built by tabix-slicing the gnomAD v4.1 sites VCFs to the
    cardiac panel intervals (multi-TB source → slice at deploy).

ABSENCE SEMANTICS (criteria-critical — drives PM2/BA1/BS1):
  * freq DB present AND variant's chrom-pos is INSIDE a panel interval but NOT
    in the table → variant genuinely absent from gnomAD over the panel slice →
    ``variant_found=False, variant=None, ok=True`` and NO live call (correct
    PM2 "absent" signal).
  * variant OFF-panel, OR freq DB absent/unbuilt → live fallback
    (``_fetch_gnomad_live``), because the local slice can't speak to a variant
    outside its coverage. Fail safe to live rather than report a wrong (absent)
    frequency.

The public ``fetch_gnomad(variant_id, gene, indel_hgvs=None)`` signature and
the return-dict shape (keys, the nested ``variant``/``gene`` GraphQL blobs, and
the optional indel ``resolved_by``/``indel_rsid``/``indel_unresolved``/
``indel_spdi`` fields) are byte-identical to the previous live-only client —
app.py, the prompt builder, and the frontend consume specific keys. Never raises.

The indel rsID fallback is LIVE only when the annotation path is live: it needs
Ensembl's variant_recoder, so it is skipped under ``HEARTVAR_VEP_OFFLINE`` (and
under ``HEARTVAR_OFFLINE_STRICT``), leaving the frequency UNRESOLVED — the same
conservative outcome the call itself produces when Ensembl is down, which on
2026-08-31 it was, at a cost of 29 s per indel. See ``_reconcile_indel``.

Env overrides (resolved at call time so tests/env win):
  * ``GNOMAD_CONSTRAINT_DB_PATH`` → constraint DB (default data/gnomad_constraint.db)
  * ``GNOMAD_FREQ_DB_PATH``       → frequency DB (default data/gnomad_freq.db)

Rebuild the caches on a new gnomAD release:
    python3 scripts/build_gnomad_constraint_db.py
    python3 scripts/build_gnomad_freq_db.py
"""

from __future__ import annotations

import json
import logging
import os
import sqlite3
from pathlib import Path

import httpx

from ..localio import run_local
from ._http_retry import make_async_client, request_with_retry
from ._offline import offline_strict, vep_offline_enabled
from ._paths import PROJECT_ROOT

log = logging.getLogger("heartvar.gnomad")

GNOMAD_API = "https://gnomad.broadinstitute.org/api"

_PROJECT_ROOT = PROJECT_ROOT
_DEFAULT_CONSTRAINT_DB = _PROJECT_ROOT / "data" / "gnomad_constraint.db"
_DEFAULT_FREQ_DB = _PROJECT_ROOT / "data" / "gnomad_freq.db"

_CONSTRAINT_FIELDS = (
    "pLI", "oe_lof", "oe_lof_upper", "oe_mis", "oe_mis_upper", "mis_z", "syn_z",
)


def _constraint_db_path() -> Path:
    """Resolve the constraint DB path at CALL TIME so env / tests win."""
    override = os.environ.get("GNOMAD_CONSTRAINT_DB_PATH")
    return Path(override) if override else _DEFAULT_CONSTRAINT_DB


def _freq_db_path() -> Path:
    """Resolve the frequency DB path at CALL TIME so env / tests win."""
    override = os.environ.get("GNOMAD_FREQ_DB_PATH")
    return Path(override) if override else _DEFAULT_FREQ_DB


_PANEL_BY_CHROM: dict[str, list[tuple[int, int]]] | None | bool = None


def _load_panel_intervals() -> dict[str, list[tuple[int, int]]] | None:
    """Load the cardiac panel intervals once into a {chrom: [(start, end)…]}
    map. Returns None if the panel can't be resolved (then every variant is
    treated as OFF-panel → live fallback, which is the fail-safe direction).

    Reads the cardiac-panel BED directly (``data/cardiac_panel.grch38.bed`` — a
    deploy-built cache, the SAME file SpliceAI reads). A missing/unreadable BED
    never breaks the freq lookup: every variant is then treated as off-panel and
    routed to live (the fail-safe direction). No network and no import from the
    build-only ``scripts/`` package, which is excluded from the deploy image
    (.dockerignore) — so this resolves correctly inside the container.
    """
    global _PANEL_BY_CHROM
    if _PANEL_BY_CHROM is not None and _PANEL_BY_CHROM is not False:
        return _PANEL_BY_CHROM  # type: ignore[return-value]
    if _PANEL_BY_CHROM is False:
        return None
    bed = _PROJECT_ROOT / "data" / "cardiac_panel.grch38.bed"
    if not bed.is_file():
        _PANEL_BY_CHROM = False
        return None
    try:
        by_chrom: dict[str, list[tuple[int, int]]] = {}
        for raw in bed.read_text().splitlines():
            row = raw.strip()
            if not row or row.startswith("#"):
                continue
            cols = row.split("\t")
            if len(cols) < 3:
                continue
            try:
                start, end = int(cols[1]), int(cols[2])
            except ValueError:
                continue
            by_chrom.setdefault(_norm_chrom(cols[0]), []).append((start, end))
        _PANEL_BY_CHROM = by_chrom
        return by_chrom
    except OSError:
        log.warning("gnomAD: cardiac panel BED unreadable at %s; "
                    "freq lookups fall back to live for all variants", bed)
        _PANEL_BY_CHROM = False
        return None


def _norm_chrom(chrom: str) -> str:
    """Strip a leading ``chr`` so panel BED ('7') and variant_id ('7') agree."""
    c = str(chrom).strip()
    return c[3:] if c.lower().startswith("chr") else c


def _parse_variant_id(variant_id: str) -> tuple[str, int] | None:
    """Parse ``chrom-pos-ref-alt`` (GRCh38, no chr prefix) → (chrom, pos).
    Returns None on any malformation (then we fail safe to live)."""
    parts = (variant_id or "").split("-")
    if len(parts) < 4:
        return None
    chrom = _norm_chrom(parts[0])
    try:
        pos = int(parts[1])
    except (ValueError, TypeError):
        return None
    return chrom, pos


def _in_panel(variant_id: str) -> bool | None:
    """True if the variant's chrom-pos is inside a cardiac-panel interval,
    False if it's off-panel, None if membership can't be determined (panel
    unavailable / unparseable id) → caller routes to live (fail safe)."""
    by_chrom = _load_panel_intervals()
    if by_chrom is None:
        return None
    parsed = _parse_variant_id(variant_id)
    if parsed is None:
        return None
    chrom, pos = parsed
    for start, end in by_chrom.get(chrom, ()):
        if start <= pos <= end:
            return True
    return False


def _constraint_row_sync(gene: str) -> dict | None:
    """Read the gene's constraint row from the local DB and build the
    GraphQL ``gene`` block ``{gene_id, symbol, gnomad_constraint{…}}``.

    Returns None when the DB is absent or the gene isn't in it → caller falls
    back to the live ``gene`` query. Never raises.
    """
    db_path = _constraint_db_path()
    if not db_path.exists():
        return None
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    except sqlite3.Error:
        return None
    try:
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        cur.execute(
            "SELECT gene_symbol, gene_id, pLI, oe_lof, oe_lof_upper, "
            "oe_mis, oe_mis_upper, mis_z, syn_z "
            "FROM gnomad_constraint WHERE UPPER(gene_symbol) = UPPER(?) "
            "LIMIT 1",
            (gene,),
        )
        row = cur.fetchone()
    except sqlite3.Error:
        log.exception("gnomAD constraint local DB query failed for %s", gene)
        return None
    finally:
        conn.close()
    if row is None:
        return None
    return {
        "gene_id": row["gene_id"],
        "symbol": row["gene_symbol"],
        "gnomad_constraint": {f: row[f] for f in _CONSTRAINT_FIELDS},
    }


def _freq_payload_sync(variant_id: str) -> dict | None:
    """Read the variant's frequency payload (the exact GraphQL ``variant``
    blob) from the local DB. Returns the parsed dict on a hit, None on a clean
    miss / DB-absent / parse error. Never raises.

    The caller distinguishes "DB absent" from "in-panel miss" via
    :func:`_freq_db_exists` + :func:`_in_panel` BEFORE calling this, so a None
    here only means "row not present in an existing table".
    """
    db_path = _freq_db_path()
    if not db_path.exists():
        return None
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    except sqlite3.Error:
        return None
    try:
        cur = conn.cursor()
        cur.execute(
            "SELECT payload FROM gnomad_freq WHERE variant_id = ? LIMIT 1",
            (variant_id,),
        )
        row = cur.fetchone()
    except sqlite3.Error:
        log.exception("gnomAD freq local DB query failed for %s", variant_id)
        return None
    finally:
        conn.close()
    if row is None or row[0] is None:
        return None
    try:
        payload = json.loads(row[0])
    except (ValueError, TypeError):
        log.warning("gnomAD freq payload for %s is not valid JSON", variant_id)
        return None
    return payload if isinstance(payload, dict) else None


def _same_site_sync(chrom: str, positions: list[int]) -> list[dict]:
    """Every gnomAD allele recorded at each of ``positions`` on ``chrom``.

    Keyed by genomic position and normalised allele, NEVER by rsID: at
    chr8:11750234 the G>T allele carries no rsID in gnomAD while the G>A
    allele at the same base carries rs1205549216, so an rsID-keyed lookup
    silently drops or conflates alleles at multi-allelic sites.

    Uses a RANGE scan, not a prefix ``LIKE``. ``variant_id`` is the table's
    TEXT PRIMARY KEY, so ``>= '8-11750234-' AND < '8-11750234.'`` becomes
    SEARCH … USING COVERING INDEX and returns in ~3 ms for a whole codon;
    the equivalent ``LIKE '8-11750234-%'`` degrades to a full index SCAN of
    9.1 M rows and took 3.1 s (both measured 2026-08-24). '.' is the next
    ASCII character after '-', which bounds the position prefix exactly.
    """
    db_path = _freq_db_path()
    if not db_path.exists():
        return []
    out: list[dict] = []
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    except sqlite3.Error:
        return []
    try:
        cur = conn.cursor()
        for pos in positions:
            lo = f"{chrom}-{pos}-"
            hi = f"{chrom}-{pos}."
            try:
                rows = cur.execute(
                    "SELECT variant_id, payload FROM gnomad_freq "
                    "WHERE variant_id >= ? AND variant_id < ? "
                    "ORDER BY variant_id",
                    (lo, hi),
                ).fetchall()
            except sqlite3.Error:
                log.exception("gnomAD same-site query failed at %s:%s", chrom, pos)
                continue
            for vid, payload in rows:
                try:
                    blob = json.loads(payload) if payload else {}
                except (ValueError, TypeError):
                    continue
                if not isinstance(blob, dict):
                    continue
                ex = blob.get("exome") or {}
                ge = blob.get("genome") or {}
                ac = (ex.get("ac") or 0) + (ge.get("ac") or 0)
                an = (ex.get("an") or 0) + (ge.get("an") or 0)
                faf = [
                    v for v in (
                        (ex.get("faf95") or {}).get("popmax"),
                        (ge.get("faf95") or {}).get("popmax"),
                    ) if v is not None
                ]
                parts = str(vid).split("-")
                out.append({
                    "variant_id": vid,
                    "position": pos,
                    "ref": parts[2] if len(parts) > 3 else None,
                    "alt": parts[3] if len(parts) > 3 else None,
                    "rsid": blob.get("rsid"),
                    "ac": ac,
                    "an": an,
                    "af": (ac / an) if an else None,
                    "faf95_popmax": max(faf) if faf else None,
                })
    finally:
        conn.close()
    return out


async def fetch_same_site_frequencies(
    chrom: str | None, positions: list[int] | None,
) -> dict:
    """gnomAD alleles at the variant's own base and the rest of its codon.

    Local-only by design. The frequency DB covers the cardiac panel
    intervals; off-panel genes degrade to
    ``{"ok": True, "available": False, "reason": "off-panel"}`` rather than
    firing a live per-position gnomAD query, which would be a new network
    round-trip per base for a context panel.
    """
    if not chrom or not positions:
        return {"ok": True, "available": False, "reason": "no codon span"}
    chrom = _norm_chrom(chrom)
    if not _freq_db_exists():
        return {"ok": True, "available": False, "reason": "frequency DB absent"}
    in_panel = _in_panel(f"{chrom}-{positions[0]}-N-N")
    if in_panel is False:
        return {
            "ok": True, "available": False, "reason": "off-panel",
            "note": ("same-site data not available off-panel — the local gnomAD "
                     "frequency DB covers the cardiac panel intervals only"),
        }
    alleles = await run_local(_same_site_sync, chrom, list(positions))
    return {
        "ok": True,
        "available": True,
        "chrom": chrom,
        "positions": list(positions),
        "alleles": alleles,
        "count": len(alleles),
    }


def _freq_db_exists() -> bool:
    return _freq_db_path().exists()


async def fetch_gnomad(
    variant_id: str | None, gene: str, indel_hgvs: str | None = None,
    variant_id_canonical: bool = False,
) -> dict:
    """Local-first gnomAD lookup: variant frequency + gene constraint.

    ``variant_id`` is "chrom-pos-ref-alt" (GRCh38, no chr prefix). If None,
    only constraint is fetched (constraint-only path).

    ``variant_id_canonical`` — True when ``variant_id`` was built from VEP's
    left-aligned ``vcf_string`` (matches gnomAD's own anchored key). It gates
    whether an in-panel indel MISS may be declared genuinely absent: a miss on
    a non-canonical hand-built indel key is meaningless, so absence is only
    trusted for a canonical key (see ``_resolve_indel_freq_by_rsid``).

    Resolution order:
      1. gene CONSTRAINT served from the local constraint DB; live fallback on
         a DB-absent / gene-not-found miss.
      2. variant FREQUENCY served from the local freq DB with strict absence
         semantics (see module docstring): in-panel-missing → genuinely absent
         (no live call); off-panel or DB-absent → live fallback.
      3. The indel rsID fallback stays LIVE (gnomAD canonicalises the indel id
         from the rsID in a way we can't derive locally).

    Return shape is byte-identical to the previous live-only client. Never
    raises.
    """
    if not variant_id:
        return await _constraint_only(gene)

    gene_block = await run_local(_constraint_row_sync, gene)
    constraint_local_hit = gene_block is not None

    freq_db_present = await run_local(_freq_db_exists)
    if freq_db_present:
        payload = await run_local(_freq_payload_sync, variant_id)
        if payload is not None:
            if gene_block is None:
                gene_block = await _fetch_constraint_block_live(gene)
            return {
                "ok": True,
                "variant_id": variant_id,
                "variant_found": True,
                "variant": payload,
                "gene": gene_block,
                "errors": None,
            }
        in_panel = await run_local(_in_panel, variant_id)
        if in_panel is True:
            if gene_block is None:
                gene_block = await _fetch_constraint_block_live(gene)
            result = {
                "ok": True,
                "variant_id": variant_id,
                "variant_found": False,
                "variant": None,
                "gene": gene_block,
                "errors": None,
            }
            return await _resolve_indel_freq_by_rsid(
                result, indel_hgvs, variant_id_canonical,
            )

    result = await _fetch_gnomad_live(
        variant_id, gene, indel_hgvs, variant_id_canonical,
    )
    if constraint_local_hit and gene_block is not None:
        if result.get("ok"):
            result["gene"] = gene_block
    return result


async def _constraint_only(gene: str) -> dict:
    """Constraint-only path (variant_id is None). Local first, live fallback."""
    gene_block = await run_local(_constraint_row_sync, gene)
    if gene_block is not None:
        return {
            "ok": True,
            "variant_id": None,
            "variant_found": False,
            "variant": None,
            "gene": gene_block,
        }
    return await _fetch_constraint_only_live(gene)


VARIANT_QUERY = """
query VariantAndGene($variantId: String!, $gene: String!) {
  variant(variantId: $variantId, dataset: gnomad_r4) {
    variantId
    rsid
    exome {
      ac
      an
      af
      ac_hom
      ac_hemi
      filters
      populations { id ac an ac_hom ac_hemi }
      faf95 { popmax popmax_population }
    }
    genome {
      ac
      an
      af
      ac_hom
      ac_hemi
      filters
      populations { id ac an ac_hom ac_hemi }
      faf95 { popmax popmax_population }
    }
  }
  gene(gene_symbol: $gene, reference_genome: GRCh38) {
    gene_id
    symbol
    gnomad_constraint {
      pLI
      oe_lof
      oe_lof_upper
      oe_mis
      oe_mis_upper
      mis_z
      syn_z
    }
  }
}
"""


RSID_QUERY = """
query VariantByRsid($rsid: String!) {
  variant(rsid: $rsid, dataset: gnomad_r4) {
    variantId
    rsid
    exome { ac an af ac_hom ac_hemi filters populations { id ac an ac_hom ac_hemi } faf95 { popmax popmax_population } }
    genome { ac an af ac_hom ac_hemi filters populations { id ac an ac_hom ac_hemi } faf95 { popmax popmax_population } }
  }
}
"""


async def _query_variant_by_rsid(rsid: str) -> dict | None:
    """Return the gnomAD `variant` blob for an rsID, or None on any failure /
    no match. Never raises."""
    if offline_strict():
        return None
    try:
        async with make_async_client() as client:
            r = await request_with_retry(
                client, "POST", GNOMAD_API,
                json={"query": RSID_QUERY, "variables": {"rsid": rsid}},
                timeout=20.0, name="gnomAD/rsid",
            )
        if r is None or r.status_code != 200:
            return None
        return (r.json().get("data") or {}).get("variant")
    except (httpx.HTTPError, ValueError):
        return None


async def _resolve_indel_freq_by_rsid(
    result: dict, indel_hgvs: str | None, variant_id_canonical: bool = False,
) -> dict:
    """INDEL/dup frequency fallback via dbSNP rsID. A hand-built
    ``chrom-pos-ref-alt`` id is unreliable for indels (gnomAD keys on a
    left-aligned, anchored id), so when a variant lookup found NOTHING and this
    is an indel, resolve the variant's dbSNP rsID (Ensembl variant_recoder) and
    re-query gnomAD by rsID — gnomAD then matches its own canonical id. Recovers
    BA1/BS1 (and stops the spurious "absent → PM2") on common indels.

    Mutates + returns ``result``. When the frequency still can't be confirmed it
    sets ``indel_unresolved=True`` so the criteria layer treats the variant as
    frequency-UNKNOWN rather than confidently absent. No-op when the variant was
    already found or this is not an indel. Used by BOTH the live path and the
    local-first in-panel-miss path. Never raises.

    ``variant_id_canonical`` — when the primary lookup key was VEP's canonical
    left-aligned ``vcf_string`` (which matches gnomAD's anchored key), an
    in-panel miss IS meaningful, so the indel may be declared GENUINELY ABSENT
    (``indel_unresolved`` NOT set → PM2_Supporting fires). Genuine absence is
    only trusted when ALL hold:
      - ``variant_id_canonical`` is True (a non-canonical key could never have
        matched, so its miss proves nothing), AND
      - the recoder call SUCCEEDED (``rec['ok'] is True``) — a FAILED/offline
        recoder call is NOT evidence that no dbSNP entry exists and MUST keep
        ``indel_unresolved=True`` (else a common indel hitting a failed call
        offline would be wrongly declared absent → wrong PM2 + wrong BA1/BS1
        suppression = catastrophic FP), AND
      - the recoder returned NO rsID (``rsid is None``) — a positive "this
        variant has no dbSNP entry" answer.
    Anything else leaves the frequency UNRESOLVED (conservative).
    """
    if result.get("variant_found") or not indel_hgvs:
        return result
    if offline_strict():
        result["indel_unresolved"] = True
        return result
    if vep_offline_enabled():
        result["indel_unresolved"] = True
        return result
    from .ensembl_vep import fetch_variant_recoder_rsid
    rec = await fetch_variant_recoder_rsid(indel_hgvs)
    rsid = rec.get("rsid")
    if rsid:
        rs_variant = await _query_variant_by_rsid(rsid)
        if rs_variant:
            result["variant"] = rs_variant
            result["variant_found"] = True
            result["resolved_by"] = "rsid"
            result["indel_rsid"] = rsid
    if not result["variant_found"]:
        confirmed_absent = (
            variant_id_canonical
            and rec.get("ok") is True
            and rsid is None
        )
        if not confirmed_absent:
            result["indel_unresolved"] = True
        result["indel_spdi"] = rec.get("spdi")
    return result


async def _fetch_gnomad_live(
    variant_id: str | None, gene: str, indel_hgvs: str | None = None,
    variant_id_canonical: bool = False,
) -> dict:
    """Live gnomAD GraphQL query for variant frequencies and gene constraint.
    This is the original ``fetch_gnomad`` body, kept verbatim as the fallback.

    `variant_id` is "chrom-pos-ref-alt" (GRCh38). If None, only constraint is
    fetched.

    `indel_hgvs` (the variant's transcript:c. HGVS) enables an INDEL/dup
    fallback: a hand-built chrom-pos-ref-alt id is unreliable for indels (gnomAD
    keys on a left-aligned, anchored id we can't derive without the reference),
    so when the primary lookup finds nothing AND this is an indel, we resolve the
    variant's dbSNP rsID (Ensembl variant_recoder) and re-query gnomAD by rsID —
    gnomAD then matches its own canonical id. This recovers BA1/BS1 (and stops the
    spurious "absent → PM2") on common indels that previously read as not-found.
    """
    if offline_strict():
        return {"ok": False, "error": "gnomAD unavailable (offline-strict mode; live query disabled)"}

    if not variant_id:
        return await _fetch_constraint_only_live(gene)

    async with make_async_client() as client:
        r = await request_with_retry(
            client, "POST", GNOMAD_API,
            json={"query": VARIANT_QUERY, "variables": {"variantId": variant_id, "gene": gene}},
            timeout=20.0, name="gnomAD/variant",
        )
        if r is None:
            return {"ok": False, "error": "gnomAD transient failure after retries"}
        if r.status_code != 200:
            return {"ok": False, "error": f"{r.status_code}: {r.text[:200]}"}
        payload = r.json()
        data = payload.get("data") or {}
        errors = payload.get("errors")
        variant = data.get("variant")
        gene_data = data.get("gene")
        result = {
            "ok": True,
            "variant_id": variant_id,
            "variant_found": variant is not None,
            "variant": variant,
            "gene": gene_data,
            "errors": errors,
        }

    return await _resolve_indel_freq_by_rsid(
        result, indel_hgvs, variant_id_canonical,
    )


_CONSTRAINT_QUERY = """
query Gene($gene: String!) {
  gene(gene_symbol: $gene, reference_genome: GRCh38) {
    gene_id
    symbol
    gnomad_constraint {
      pLI oe_lof oe_lof_upper oe_mis oe_mis_upper mis_z syn_z
    }
  }
}
"""


async def _fetch_constraint_only_live(gene: str) -> dict:
    """Original ``_constraint_only`` body — the live constraint-only fallback."""
    if offline_strict():
        return {"ok": False, "error": "gnomAD unavailable (offline-strict mode; live query disabled)"}

    async with make_async_client() as client:
        r = await request_with_retry(
            client, "POST", GNOMAD_API,
            json={"query": _CONSTRAINT_QUERY, "variables": {"gene": gene}},
            timeout=30.0, name="gnomAD/constraint",
        )
        if r is None:
            return {"ok": False, "error": "gnomAD transient failure after retries"}
        if r.status_code != 200:
            return {"ok": False, "error": f"{r.status_code}: {r.text[:200]}"}
        data = r.json().get("data") or {}
        return {
            "ok": True,
            "variant_id": None,
            "variant_found": False,
            "variant": None,
            "gene": data.get("gene"),
        }


async def _fetch_constraint_block_live(gene: str) -> dict | None:
    """Fetch just the ``gene`` GraphQL block (for the variant path when the
    local constraint DB missed but the freq cache hit). Returns the gene dict
    or None on failure. Never raises."""
    if offline_strict():
        return None
    try:
        res = await _fetch_constraint_only_live(gene)
    except Exception:  # noqa: BLE001
        return None
    return res.get("gene") if res.get("ok") else None
