"""phyloP100way conservation lookup for BP7's "nucleotide not highly conserved".

BP7 needs a per-base conservation score, but the dbNSFP plugin HeartVar reads
via VEP only covers nonsynonymous + splice-site SNVs — so phyloP is absent for
pure *synonymous* variants (BP7's primary population) and most deep-intronic /
UTR sites. This client backfills phyloP100way (GRCh38) for exactly those
positions, by genomic coordinate, from a bigWig:

  * LOCAL-FIRST: ``HEARTVAR_PHYLOP_PATH`` → a local hg38.phyloP100way bigWig
    (full ~9.2 GB, or a panel-scoped extract). This is the deploy path — no
    network, exact per-base values, same numbers as dbNSFP.
  * LIVE FALLBACK: the UCSC hosted bigWig over HTTPS range reads. Accurate
    (verified identical to dbNSFP), zero storage, used when no local file is
    configured. (The offline-VEP deploy supersedes this with ``--custom``.)

Reads are cached (static data → long TTL) and single-flighted. pyBigWig is an
optional C-extension dependency: if it is not installed, or the lookup fails,
the function returns ``None`` and BP7 simply no-ops on conservation (its prior
behaviour) — never raises.
"""
from __future__ import annotations

import logging
import os

from ..localio import run_local
from ._cache import EXTERNAL_CACHE, TTL_PMC
from ._offline import offline_strict

log = logging.getLogger("heartvar.conservation")

# UCSC hg38 100-way vertebrate phyloP, license-free (https://genome.ucsc.edu/license/).
_UCSC_PHYLOP_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw"
)


def _ucsc_chrom(chrom: str) -> str:
    """VEP seq_region_name ('14', 'X', 'MT') → UCSC contig ('chr14', 'chrX', 'chrM')."""
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        return c
    if c.upper() in ("MT", "M"):
        return "chrM"
    return f"chr{c}"


_PHYLOP_FAIL = object()


def _phylop_lookup_sync(chrom: str, pos: int):
    """Blocking single-base bigWig read. Wrapped in run_local by the
    caller. ``pos`` is 1-based (VEP/HGVS convention); bigWig is 0-based half-open.
    Returns a float, ``None`` (read OK but no data at that base), or
    ``_PHYLOP_FAIL`` (transient — do not cache)."""
    try:
        import pyBigWig
    except Exception:  # noqa: BLE001
        log.debug("pyBigWig not installed — conservation lookup skipped")
        return _PHYLOP_FAIL
    local_path = os.environ.get("HEARTVAR_PHYLOP_PATH")
    if not local_path and offline_strict():
        log.debug("offline-strict: HEARTVAR_PHYLOP_PATH unset — skipping UCSC phyloP fetch")
        return _PHYLOP_FAIL
    path = local_path or _UCSC_PHYLOP_URL
    bw = None
    try:
        bw = pyBigWig.open(path)
        if bw is None:
            return _PHYLOP_FAIL
        uc = _ucsc_chrom(chrom)
        if uc not in (bw.chroms() or {}):
            return None
        vals = bw.values(uc, pos - 1, pos)
        if not vals:
            return None
        v = vals[0]
        return None if v is None or v != v else float(v)
    except Exception as e:  # noqa: BLE001
        log.debug("phyloP lookup failed for %s:%s — %r", chrom, pos, e)
        return _PHYLOP_FAIL
    finally:
        if bw is not None:
            try:
                bw.close()
            except Exception:  # noqa: BLE001
                pass


async def fetch_phylop100way(chrom: str | None, pos: int | None) -> float | None:
    """phyloP100way (GRCh38) at a 1-based genomic position, or ``None`` if
    unavailable. Cached (static data) + single-flighted; never raises."""
    if not chrom or not pos:
        return None
    try:
        pos = int(pos)
    except (TypeError, ValueError):
        return None
    key = ("phylop100way", _ucsc_chrom(chrom), pos)
    val = await EXTERNAL_CACHE.get_or_set(
        key,
        lambda: run_local(_phylop_lookup_sync, chrom, pos),
        ttl=TTL_PMC,
        should_cache=lambda v: v is not _PHYLOP_FAIL,
    )
    return None if val is _PHYLOP_FAIL else val
