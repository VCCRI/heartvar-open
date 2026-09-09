"""AlphaMissense local lookup — pre-computed missense pathogenicity scores.

AlphaMissense (Cheng et al. 2023, DeepMind) is a per-missense pathogenicity
predictor with calibrated thresholds:

  * ``score > 0.564`` → likely_pathogenic — supports PP3
  * ``score < 0.340`` → likely_benign     — supports BP4
  * in-between        → ambiguous         — neither PP3 nor BP4

The hg38 single-substitution release ships as a tab-separated table,
bgzipped, with a tabix index. Columns in order:

    CHROM  POS  REF  ALT  genome  uniprot_id  transcript_id
    protein_variant  am_pathogenicity  am_class

The file uses chr-prefixed chromosomes (e.g. ``chr3``). Coordinates are
1-based. Only missense substitutions are present — frameshifts, stop
gains, splice variants, etc. have no row, which the caller surfaces as
``not_applicable: True`` rather than an error.

Path is read from the ``ALPHAMISSENSE_PATH`` env var; see ``.env.example``
for download + indexing instructions.

Non-fatal: any error returns ``{available: False, reason: …}``.
"""

from __future__ import annotations

import os
from typing import Any

from ..localio import run_local
from .ensembl_vep import _complement


def _na() -> dict:
    """Variant not in the file — expected for non-missense variants."""
    return {
        "available": True,
        "score": None,
        "classification": None,
        "protein_variant": None,
        "not_applicable": True,
    }


def _unavailable(reason: str) -> dict:
    """File missing / lookup failure — distinct from a clean miss."""
    return {
        "available": False,
        "reason": reason,
        "score": None,
        "classification": None,
        "protein_variant": None,
        "not_applicable": False,
    }


def _lookup_sync(path: str, chrom: str, pos: int, ref: str, alt: str) -> dict:
    """Synchronous tabix lookup. Wrapped in run_local by the caller
    so the (sub-millisecond, but blocking) C call doesn't stall the event
    loop."""
    try:
        import pysam
    except ImportError:
        return _unavailable("pysam not installed")

    tbx = None
    try:
        tbx = pysam.TabixFile(path)
        for row in tbx.fetch(chrom, pos - 1, pos):
            parts = row.split("\t")
            if len(parts) < 10:
                continue
            r_chrom, r_pos, r_ref, r_alt = parts[0], parts[1], parts[2], parts[3]
            if r_chrom != chrom or r_pos != str(pos):
                continue
            if r_ref != ref or r_alt != alt:
                continue
            try:
                score = float(parts[8])
            except ValueError:
                continue
            return {
                "available": True,
                "score": round(score, 3),
                "classification": parts[9].strip(),
                "protein_variant": parts[7].strip(),
                "not_applicable": False,
            }
        return _na()
    except (ValueError, OSError) as e:
        return _unavailable(f"tabix lookup error: {e!r}")
    finally:
        if tbx is not None:
            try:
                tbx.close()
            except Exception:
                pass


async def fetch_alphamissense(vep_data: dict[str, Any]) -> dict:
    """Look up the variant in the local AlphaMissense TSV.

    Args:
        vep_data: The dict returned by ``fetch_vep`` — supplies
            ``seq_region_name``, ``start``, ``allele_string``, and
            ``strand``. Strand is honoured so reverse-strand genes
            (MYH7, RAF1, etc.) match the forward-strand REF/ALT the
            TSV indexes on.

    Returns:
        {"available": True, "score": float|None, "classification": str|None,
         "protein_variant": str|None, "not_applicable": bool}
        on success — ``not_applicable=True`` for variants outside the
        single-substitution release (frameshifts, stop gains, splice
        variants, etc.), which is the expected outcome rather than an
        error.

        {"available": False, "reason": "…"} when the file is missing /
        unreadable or the VEP result is unusable.

    Never raises.
    """
    path = os.environ.get("ALPHAMISSENSE_PATH")
    if not path:
        return _unavailable("ALPHAMISSENSE_PATH not configured")
    if not os.path.isfile(path):
        return _unavailable(f"file not found: {path}")
    if not vep_data or not vep_data.get("ok"):
        return _unavailable("no VEP data")

    chrom = vep_data.get("seq_region_name")
    pos = vep_data.get("start")
    allele = vep_data.get("allele_string") or ""
    if not chrom or not pos or "/" not in allele:
        return _unavailable("missing chrom/pos/allele in VEP result")

    ref, alt = allele.split("/", 1)
    if vep_data.get("strand") == -1:
        ref, alt = _complement(ref), _complement(alt)

    chrom_q = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
    try:
        pos_i = int(pos)
    except (TypeError, ValueError):
        return _unavailable(f"invalid POS: {pos!r}")

    return await run_local(_lookup_sync, path, chrom_q, pos_i, ref, alt)
