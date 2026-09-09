"""SpliceAI delta-score lookup — LOCAL panel slice first, Broad API fallback.

Local-first
-----------
For GRCh38 in-panel SNVs, ``fetch_spliceai`` reads a small bgzipped + tabix-
indexed VCF sliced from Ensembl's precomputed masked-SNV SpliceAI release
(built by ``scripts/build_spliceai_db.py`` →
``data/spliceai_cardiac.masked.grch38.vcf.gz``; path overridable via the
``SPLICEAI_DB_PATH`` env var). The sliced VCF's INFO carries one or more
comma-separated SpliceAI annotations:

    SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL

Each matching gene annotation becomes one ``scores_per_transcript`` entry
(``transcript_id=None`` — the precomputed file is per-gene). This keeps the
public deployment off the Broad's Cloud Run instance for the common case.

Live fallback
-------------
The original Broad SpliceAI Lookup client is preserved verbatim as
``_fetch_spliceai_live`` and is used whenever the local path cannot answer
*authoritatively*:

  * GRCh37 (the slice is GRCh38-only),
  * the slice file is absent,
  * the variant is OFF-PANEL (outside every cardiac-panel gene span), or
  * panel membership cannot be determined (BED absent / unreadable).

An IN-PANEL GRCh38 position with NO SpliceAI record in the slice is a real
answer ("no score at this position"), NOT a fallback — it returns
``{ok:False, lookup_failed:False, reason:"No SpliceAI score at this
position"}`` to match the live "no scores" branch. The Broad endpoints:

  hg38: https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/
  hg19: https://spliceai-37-xwkwwwxdwq-uc.a.run.app/spliceai/

Query parameters:
  variant   chrom-pos-ref-alt (1-based, no "chr" prefix required)
  hg        "38" or "37"
  distance  base-pairs scanned around the variant (default 500)
  mask      "0" (raw) or "1" (mask scores below 0.1)
  raw       "0" or "1" — raw=1 returns the unfiltered numeric scores

Response (success) — `scores` is a list, one entry per overlapping
transcript:

  {
    "variant": "...",
    "scores": [
      {
        "DS_AG": 0.02, "DS_AL": 0.01, "DS_DG": 0.99, "DS_DL": 0.99,
        "DP_AG": -10, "DP_AL": -10, "DP_DG": 1, "DP_DL": 0,
        "gene_name": "MYH7", "transcript_id": "...", ...
      }
    ]
  }

Response (no scores) — the API returns an `error` string when the
variant does not overlap a GENCODE-basic exon/intron, or when the
variant cannot be parsed. We surface that as ok=False with the error
message so downstream PVS1/BP7 logic falls back to insufficient_data.

CRITERIA-CRITICAL: a wrong delta score = a dangerous PP3/BP4/BP7/PVS1
misclassification. On ANY ambiguity in the local path (REF/ALT not matched
unambiguously, malformed INFO, unreadable index) we fail safe to the live
fallback rather than returning a possibly-wrong score.

The client never raises — any HTTPError / timeout / JSON / tabix error is
caught and returned as ok=False / falls back so the DB-gather pipeline
stays robust.
"""

from __future__ import annotations

import asyncio
import logging
import os
from pathlib import Path

import httpx

from ..localio import run_local
from ._http_retry import _with_connect_cap, make_async_client
from ._offline import offline_strict
from ._throttle import global_semaphore, throttle_host
from ._paths import PROJECT_ROOT

log = logging.getLogger("heartvar.spliceai")

_PROJECT_ROOT = PROJECT_ROOT
DB_PATH = _PROJECT_ROOT / "data" / "spliceai_cardiac.masked.grch38.vcf.gz"
PANEL_BED_PATH = _PROJECT_ROOT / "data" / "cardiac_panel.grch38.bed"

_PANEL_INTERVALS: dict[str, list[tuple[int, int]]] | None = None


def _db_path() -> Path:
    """Resolve the slice path at call time so the SPLICEAI_DB_PATH env
    override (or a test monkeypatching the module attribute) always wins."""
    return Path(os.environ.get("SPLICEAI_DB_PATH") or DB_PATH)


def _bed_path() -> Path:
    """Resolve the panel-BED path at call time (env override honoured)."""
    return Path(os.environ.get("SPLICEAI_PANEL_BED_PATH") or PANEL_BED_PATH)


def _reset_panel_cache() -> None:
    """Clear the cached panel intervals — used by tests so a monkeypatched
    BED path is re-read instead of a stale closure."""
    global _PANEL_INTERVALS
    _PANEL_INTERVALS = None

_SPLICEAI_MAX_ATTEMPTS = 3
_SPLICEAI_BACKOFF = (1.0, 2.0)

API_BY_BUILD = {
    "GRCh38": "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/",
    "GRCh37": "https://spliceai-37-xwkwwwxdwq-uc.a.run.app/spliceai/",
}

TIMEOUT = 25.0


def _to_float(v):
    """Coerce a numeric or numeric-string value to float, returning None on
    failure. The Broad API serialises delta scores as strings ("0.99"); older
    docs showed them as numbers. Accept both forms defensively."""
    if v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _max_delta(score_entry: dict) -> float | None:
    """Pick the largest of the four delta scores from one transcript entry."""
    candidates = []
    for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL"):
        f = _to_float(score_entry.get(k))
        if f is not None:
            candidates.append(f)
    return max(candidates) if candidates else None


_NUC_CHARS = set("ACGTN")


def _is_indel_variant_id(variant_id: str) -> tuple[bool, str | None]:
    """Detect indel-format gnomAD variant IDs that the Broad SpliceAI API
    cannot consume. The API requires REF and ALT to be valid nucleotide
    strings; VEP serialises deletions as "-" for the missing allele which
    blows up the API parser. Returns (is_indel, reason)."""
    if not variant_id:
        return False, None
    parts = variant_id.split("-")
    ref = alt = None
    if len(parts) == 4:
        ref, alt = parts[2], parts[3]
    elif len(parts) == 5 and parts[4] == "":
        ref, alt = parts[2], parts[3]
    else:
        return False, None
    if not ref or ref == "-" or any(c.upper() not in _NUC_CHARS for c in ref):
        return True, "Insertion variant" if (not ref or ref == "-") else None
    if not alt or alt == "-" or any(c.upper() not in _NUC_CHARS for c in alt):
        return True, "Deletion variant" if (not alt or alt == "-") else None
    return False, None


async def _fetch_spliceai_live(variant_id: str | None, assembly: str | None) -> dict:
    """Live Broad SpliceAI Lookup API call — the original implementation,
    preserved verbatim as the graceful fallback for the local-first
    ``fetch_spliceai``.

    `variant_id` is "chrom-pos-ref-alt" in gnomAD style (matches what
    Ensembl VEP returns via `gnomad_variant_id`). `assembly` is the VEP
    assembly_name ("GRCh38" or "GRCh37"); when missing we default to
    GRCh38.

    Returns a dict with at minimum {"ok": bool, "lookup_failed": bool}.
    The two-flag scheme lets the UI distinguish three distinct outcomes:

      ok=True, lookup_failed=False
        Real scores returned. Includes `max_delta`,
        `scores_per_transcript`, `model_message`.
      ok=False, lookup_failed=False
        API responded successfully but had no scores at this position
        (variant outside GENCODE-basic exons/introns, empty scores
        list, etc.). Includes `reason`. This is a *real answer*, not a
        failure — the loading-screen row should tick green.
      ok=False, lookup_failed=True
        Transport-level failure: timeout, HTTP non-200, JSON parse
        error. Includes `error`. The loading-screen row should mark
        this as failed since we don't know what the score would be.

    For indel variants the Broad API cannot parse (REF or ALT is "-" or
    empty), returns {"skipped": True, "not_applicable": True, ...} so
    the frontend renders the row as "Not applicable" rather than red.
    """
    if not variant_id:
        return {
            "ok": False,
            "lookup_failed": True,
            "error": "no variant id (VEP did not resolve coords)",
        }

    is_indel, indel_kind = _is_indel_variant_id(variant_id)
    if is_indel:
        reason = (
            f"{indel_kind or 'Indel variant'} — SpliceAI lookup not supported "
            f"for this allele format"
        )
        log.info("SpliceAI %s skipped — %s", variant_id, reason)
        return {
            "skipped": True,
            "not_applicable": True,
            "variant_id": variant_id,
            "reason": reason,
        }

    raw_build = (assembly or "GRCh38").strip()
    build = next(
        (k for k in API_BY_BUILD if k.lower() == raw_build.lower()),
        "GRCh38",
    )
    url = API_BY_BUILD[build]
    params = {
        "variant": variant_id,
        "hg": "38" if build == "GRCh38" else "37",
        "distance": "500",
        "mask": "1",
        "raw": "0",
    }

    log.info(
        "SpliceAI request: %s?variant=%s&hg=%s&distance=%s&mask=%s&raw=%s",
        url, params["variant"], params["hg"], params["distance"], params["mask"], params["raw"],
    )

    payload = None
    last_err = "no attempts"
    async with make_async_client() as client:
        for attempt in range(_SPLICEAI_MAX_ATTEMPTS):
            try:
                await throttle_host(url)
                async with global_semaphore():
                    r = await client.get(url, params=params, timeout=_with_connect_cap(TIMEOUT))
                if r.status_code == 200:
                    payload = r.json()
                    break
                if r.status_code < 500:
                    msg = f"HTTP {r.status_code}: {r.text[:200]}"
                    log.warning("SpliceAI %s — %s", variant_id, msg)
                    return {
                        "ok": False, "lookup_failed": True,
                        "variant_id": variant_id, "error": msg,
                    }
                last_err = f"HTTP {r.status_code}"
            except (httpx.HTTPError, ValueError) as e:
                last_err = f"{type(e).__name__}: {e}"
            log.warning(
                "SpliceAI %s attempt %d/%d: %s",
                variant_id, attempt + 1, _SPLICEAI_MAX_ATTEMPTS, last_err,
            )
            if attempt < _SPLICEAI_MAX_ATTEMPTS - 1:
                await asyncio.sleep(_SPLICEAI_BACKOFF[min(attempt, len(_SPLICEAI_BACKOFF) - 1)])
    if payload is None:
        log.warning("SpliceAI %s — giving up after %d attempts (%s)",
                    variant_id, _SPLICEAI_MAX_ATTEMPTS, last_err)
        return {
            "ok": False, "lookup_failed": True,
            "variant_id": variant_id, "error": last_err,
        }

    if payload.get("error"):
        return {
            "ok": False,
            "lookup_failed": False,
            "variant_id": variant_id,
            "reason": payload["error"],
        }

    scores = payload.get("scores") or []
    if not scores:
        return {
            "ok": False,
            "lookup_failed": False,
            "variant_id": variant_id,
            "reason": "API returned an empty scores list",
        }

    per_transcript = []
    overall_max = 0.0
    for s in scores:
        mx = _max_delta(s)
        if mx is not None and mx > overall_max:
            overall_max = mx
        per_transcript.append({
            "gene": s.get("gene_name") or s.get("SYMBOL"),
            "transcript_id": s.get("transcript_id"),
            "DS_AG": _to_float(s.get("DS_AG")),
            "DS_AL": _to_float(s.get("DS_AL")),
            "DS_DG": _to_float(s.get("DS_DG")),
            "DS_DL": _to_float(s.get("DS_DL")),
            "DP_AG": s.get("DP_AG"),
            "DP_AL": s.get("DP_AL"),
            "DP_DG": s.get("DP_DG"),
            "DP_DL": s.get("DP_DL"),
            "max_delta": mx,
        })

    if overall_max >= 0.8:
        interp = "High confidence splice impact"
    elif overall_max >= 0.5:
        interp = "Moderate splice impact"
    elif overall_max >= 0.2:
        interp = "Low confidence splice impact"
    else:
        interp = "No predicted splice impact"

    return {
        "ok": True,
        "lookup_failed": False,
        "variant_id": variant_id,
        "assembly": build,
        "max_delta": overall_max,
        "model_message": interp,
        "scores_per_transcript": per_transcript,
    }


def _interp_for(overall_max: float) -> str:
    """Plain-English confidence label for an overall max-delta — IDENTICAL
    thresholds to the live path (≥0.8 high, ≥0.5 moderate, ≥0.2 low, else
    none). Kept as a tiny shared helper so the local slice produces the
    byte-identical ``model_message`` string the frontend + prompt expect."""
    if overall_max >= 0.8:
        return "High confidence splice impact"
    if overall_max >= 0.5:
        return "Moderate splice impact"
    if overall_max >= 0.2:
        return "Low confidence splice impact"
    return "No predicted splice impact"


def _parse_variant_id(variant_id: str) -> tuple[str, int, str, str] | None:
    """Split a gnomAD-style ``chrom-pos-ref-alt`` id into
    (chrom, pos:int, ref:UPPER, alt:UPPER). Returns ``None`` when the id is
    not a clean 4-token SNV-shaped string — the caller then fails safe to the
    live API rather than guessing.

    NOTE: indels are already filtered upstream by ``_is_indel_variant_id``;
    this only sees REF/ALT that passed that nucleotide check. We still verify
    the pos is an integer and REF/ALT are single, non-empty tokens."""
    if not variant_id:
        return None
    parts = variant_id.split("-")
    if len(parts) != 4:
        return None
    chrom, pos_s, ref, alt = parts
    if not chrom or not ref or not alt:
        return None
    try:
        pos = int(pos_s)
    except ValueError:
        return None
    return chrom, pos, ref.upper(), alt.upper()


def _load_panel_intervals() -> dict[str, list[tuple[int, int]]] | None:
    """Load the cardiac-panel intervals from the BED into
    ``{contig -> [(start, end), …]}`` (1-based inclusive), cached for the
    process lifetime.

    Returns ``None`` (caller fails safe to live) when the BED is absent or
    unreadable — we must NOT treat "couldn't read the panel" as "off-panel".
    Contigs are stored under BOTH their bare ("3") and chr-prefixed ("chr3")
    spellings so membership matches regardless of how the variant_id names the
    contig."""
    global _PANEL_INTERVALS
    if _PANEL_INTERVALS is not None:
        return _PANEL_INTERVALS
    bed = _bed_path()
    if not bed.is_file():
        return None
    intervals: dict[str, list[tuple[int, int]]] = {}
    try:
        for line in bed.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            chrom = cols[0]
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                continue
            bare = chrom[3:] if chrom.startswith("chr") else chrom
            chr_pref = chrom if chrom.startswith("chr") else f"chr{chrom}"
            for key in {bare, chr_pref}:
                intervals.setdefault(key, []).append((start, end))
    except OSError:
        log.warning("SpliceAI panel BED unreadable at %s — falling back to live", bed)
        return None
    _PANEL_INTERVALS = intervals
    return intervals


def _is_in_panel(chrom: str, pos: int) -> bool | None:
    """True/False if ``(chrom, pos)`` is inside / outside any cardiac-panel
    interval; ``None`` when membership can't be determined (BED absent), so
    the caller fails safe to live."""
    intervals = _load_panel_intervals()
    if intervals is None:
        return None
    key = chrom if chrom in intervals else (
        f"chr{chrom}" if f"chr{chrom}" in intervals else
        (chrom[3:] if chrom.startswith("chr") and chrom[3:] in intervals else None)
    )
    if key is None:
        return False
    return any(start <= pos <= end for (start, end) in intervals[key])


_SPLICEAI_FIELDS = (
    "ALLELE", "SYMBOL",
    "DS_AG", "DS_AL", "DS_DG", "DS_DL",
    "DP_AG", "DP_AL", "DP_DG", "DP_DL",
)


def _extract_spliceai_info(info: str) -> list[str]:
    """Pull the comma-separated SpliceAI annotation strings out of a VCF INFO
    column. Returns the list of pipe-delimited annotation payloads (each the
    ALLELE|SYMBOL|… string), or an empty list when the INFO has no SpliceAI
    key."""
    for field in info.split(";"):
        if field.startswith("SpliceAI="):
            payload = field[len("SpliceAI="):]
            return [a for a in payload.split(",") if a]
    return []


def _dp_to_int(v: str):
    """Coerce a delta-position token to int (the live API returns ints for
    DP_*); leave non-integers as None so a malformed field doesn't poison the
    record."""
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _local_lookup_sync(path: str, chrom: str, pos: int, ref: str, alt: str) -> dict | None:
    """Synchronous tabix read of the local SpliceAI slice. Wrapped in
    run_local by the caller.

    Returns:
      * a success dict (ok=True) when a record at ``pos`` matches REF/ALT and
        the SpliceAI ALLELE,
      * ``{"_no_score": True}`` sentinel when the position is present in the
        slice's contig but no SpliceAI record matches this variant (in-panel
        miss → the caller renders the "No SpliceAI score" no-score answer),
      * ``None`` on ANY ambiguity / read error so the caller fails safe to the
        live API (CRITERIA-CRITICAL: never return a guessed score).
    """
    try:
        import pysam
    except ImportError:
        return None

    tbx = None
    try:
        tbx = pysam.TabixFile(path)
        contigs = set(tbx.contigs)
        query_contigs = []
        for c in (chrom, f"chr{chrom}", chrom[3:] if chrom.startswith("chr") else None):
            if c and c in contigs and c not in query_contigs:
                query_contigs.append(c)
        if not query_contigs:
            if not contigs:
                detail = (
                    f"SpliceAI slice at {path} reports NO CONTIGS — it is "
                    "empty, truncated or unindexed"
                )
            else:
                detail = (
                    f"SpliceAI slice at {path} does not carry contig {chrom!r} "
                    f"(it has {len(contigs)} others) for an IN-PANEL position"
                )
            log.error("%s — reporting lookup_failed (offline by design)", detail)
            return {"_unusable": detail}

        per_transcript: list[dict] = []
        overall_max = 0.0
        matched_record = False
        for c in query_contigs:
            for row in tbx.fetch(c, pos - 1, pos):
                cols = row.split("\t")
                if len(cols) < 8:
                    continue
                r_pos_s, r_ref, r_alt, r_info = cols[1], cols[3], cols[4], cols[7]
                try:
                    r_pos = int(r_pos_s)
                except ValueError:
                    continue
                if r_pos != pos:
                    continue
                if r_ref.upper() != ref or r_alt.upper() != alt:
                    continue
                for annot in _extract_spliceai_info(r_info):
                    fields = annot.split("|")
                    if len(fields) < len(_SPLICEAI_FIELDS):
                        return None
                    rec = dict(zip(_SPLICEAI_FIELDS, fields, strict=False))
                    if rec["ALLELE"].upper() != alt:
                        continue
                    matched_record = True
                    entry = {
                        "gene": rec["SYMBOL"] or None,
                        "transcript_id": None,
                        "DS_AG": _to_float(rec["DS_AG"]),
                        "DS_AL": _to_float(rec["DS_AL"]),
                        "DS_DG": _to_float(rec["DS_DG"]),
                        "DS_DL": _to_float(rec["DS_DL"]),
                        "DP_AG": _dp_to_int(rec["DP_AG"]),
                        "DP_AL": _dp_to_int(rec["DP_AL"]),
                        "DP_DG": _dp_to_int(rec["DP_DG"]),
                        "DP_DL": _dp_to_int(rec["DP_DL"]),
                    }
                    mx = _max_delta(entry)
                    if mx is not None and mx > overall_max:
                        overall_max = mx
                    entry["max_delta"] = mx
                    per_transcript.append(entry)
            if matched_record:
                break

        if not matched_record:
            return {"_no_score": True}

        return {
            "ok": True,
            "lookup_failed": False,
            "variant_id": f"{chrom}-{pos}-{ref}-{alt}",
            "assembly": "GRCh38",
            "max_delta": overall_max,
            "model_message": _interp_for(overall_max),
            "scores_per_transcript": per_transcript,
        }
    except (ValueError, OSError) as e:
        log.warning("SpliceAI local slice read error for %s-%s-%s-%s: %r",
                    chrom, pos, ref, alt, e)
        return None
    finally:
        if tbx is not None:
            try:
                tbx.close()
            except Exception:
                pass


def spliceai_live_allowed() -> bool:
    """False by default — SpliceAI is served from the local slice ONLY.

    DECISION: no live SpliceAI lookups.
    We need to get a proper offline spliceAI lookup that works." The Broad Cloud
    Run fallback is therefore OFF unless SPLICEAI_ALLOW_LIVE is explicitly set.

    The live code is kept, not deleted, so the path stays reachable and tested
    and so the decision can be reversed with an env var rather than a patch.
    Note this is INDEPENDENT of HEARTVAR_OFFLINE_STRICT, which is off for the
    twelve clients that still want a live fallback — SpliceAI opts out on its
    own rather than forcing that global flag."""
    return os.environ.get("SPLICEAI_ALLOW_LIVE", "").strip().lower() in {
        "1", "true", "yes", "on",
    }


async def _live_or_offline(
    variant_id: str | None, assembly: str | None, *,
    reason: str, failed: bool = True,
) -> dict:
    """The single choke point for every non-local outcome in ``fetch_spliceai``.

    ⚠ THIS USED TO DEGRADE TO A FALSE "no score". When offline-strict was on it
    returned ok=False with lookup_failed=False and the words "No SpliceAI score
    at this position" — the same shape as a genuine in-panel miss — for a
    MISSING SLICE, an unreadable BED, a parse failure or a GRCh37 input. Those
    are not answers about the variant; they are faults, and painting them as a
    real negative is how a total SpliceAI outage stayed invisible in production
    until 2026-09-02.

    Each caller now passes its OWN reason, and ``failed`` says whether this is a
    fault (lookup_failed=True → the loading row goes red, and BP7/PP3/BP4/PVS1
    treat it as insufficient_data) or a legitimate coverage limit
    (lookup_failed=False → off-panel, which is a true statement about scope).

    Live is consulted only when SPLICEAI_ALLOW_LIVE is set AND offline-strict is
    not forcing local-only."""
    if spliceai_live_allowed() and not offline_strict():
        return await _fetch_spliceai_live(variant_id, assembly)
    if failed:
        log.error("SpliceAI local-only: %s (variant %s)", reason, variant_id)
    return {
        "ok": False,
        "lookup_failed": failed,
        "variant_id": variant_id,
        "reason": reason,
        **({"error": reason} if failed else {}),
    }


async def fetch_spliceai(variant_id: str | None, assembly: str | None) -> dict:
    """Look up SpliceAI delta scores — LOCAL panel slice first, Broad API
    fallback.

    Signature + return-dict shape are identical to the original (live-only)
    client; downstream PP3/BP4/BP7 + the PVS1 splice arm consume the same keys
    regardless of which path answered.

    Routing:
      1. Indel (REF/ALT is "-" / empty / non-nucleotide) → skipped /
         not_applicable (handled by the live function, which keeps the
         ``_is_indel_variant_id`` check first).
      2. GRCh38 + in-panel + slice present:
           * record matches  → local success dict (ok=True),
           * no record       → {ok:False, lookup_failed:False,
                                 reason:"No SpliceAI score at this position"}.
      3. GRCh37, off-panel, slice absent, panel-membership unknown, or ANY
         local ambiguity → live Broad API (``_fetch_spliceai_live``).
    """
    if not variant_id:
        return await _fetch_spliceai_live(variant_id, assembly)

    is_indel, _ = _is_indel_variant_id(variant_id)
    if is_indel:
        return await _fetch_spliceai_live(variant_id, assembly)

    raw_build = (assembly or "GRCh38").strip()
    if raw_build.lower() != "grch38":
        return await _live_or_offline(
            variant_id, assembly,
            reason=(f"SpliceAI slice is GRCh38-only and this variant is "
                    f"{raw_build} — no local score available"))

    parsed = _parse_variant_id(variant_id)
    if parsed is None:
        return await _live_or_offline(
            variant_id, assembly,
            reason=f"could not parse variant id {variant_id!r} for a SpliceAI lookup")
    chrom, pos, ref, alt = parsed

    db = _db_path()
    if not db.is_file():
        return await _live_or_offline(
            variant_id, assembly,
            reason=f"SpliceAI slice not found at {db} — the mount is missing it")

    in_panel = _is_in_panel(chrom, pos)
    if in_panel is None:
        return await _live_or_offline(
            variant_id, assembly,
            reason=(f"cardiac-panel BED at {_bed_path()} is absent or unreadable "
                    "— cannot tell whether this position is in panel"))
    if not in_panel:
        return await _live_or_offline(
            variant_id, assembly, failed=False,
            reason=("outside the cardiac panel — the SpliceAI slice is "
                    "panel-scoped, so no score is held for this position"))

    local = await run_local(_local_lookup_sync, str(db), chrom, pos, ref, alt)
    if local is None:
        return await _live_or_offline(
            variant_id, assembly,
            reason="SpliceAI local slice read was ambiguous or failed")
    if local.get("_unusable"):
        return {
            "ok": False,
            "lookup_failed": True,
            "variant_id": variant_id,
            "error": local["_unusable"],
            "reason": "SpliceAI local slice unusable — see error",
        }
    if local.get("_no_score"):
        return {
            "ok": False,
            "lookup_failed": False,
            "variant_id": variant_id,
            "reason": "No SpliceAI score at this position",
        }
    return local
