#!/usr/bin/env python3
"""build_spliceai_db.py — build the local SpliceAI cardiac-panel VCF slice.

Tabix-slices Ensembl's precomputed masked-SNV SpliceAI VCF down to the
cardiac-panel gene spans (CHDgene + VCEP, GRCh38, padded — see
``scripts/_cardiac_panel.py``) and writes a small bgzipped + tabix-indexed
VCF to ``data/spliceai_cardiac.masked.grch38.vcf.gz`` (+ ``.tbi``).

The backend SpliceAI client (``backend/clients/spliceai.py``) reads this
slice FIRST for GRCh38 in-panel variants and only falls back to the live
Broad SpliceAI Lookup API on an off-panel / GRCh37 / slice-absent lookup.
This keeps the public deployment from hammering the Broad Cloud Run
instance on every curation.

Source (OPEN, no login, ~27 GB; the .tbi sits next to it so we can fetch
random regions without pulling the whole file):

    https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/
      variation_plugins/
      spliceai_scores.masked.snv.ensembl_mane_v1.0.grch38.vcf.gz  (+ .tbi)

The VCF INFO carries one or more comma-separated SpliceAI annotations:

    SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL

Strategy: we do NOT download the 27 GB file. Instead we open the remote
bgzipped VCF directly with ``pysam.TabixFile`` (htslib speaks HTTPS +
range requests via the remote .tbi), fetch each padded panel interval,
write the matching records to a local plain VCF, then bgzip + tabix-index
it. Only the panel regions ever cross the wire (a few MB).

Re-running is safe — the output slice + its index are overwritten every
run. The build is idempotent.

Usage::

    python3 build_spliceai_db.py                 # slice from the remote VCF
    python3 build_spliceai_db.py --source PATH   # slice from a local copy
    python3 build_spliceai_db.py --pad 5000      # override panel padding

Data source: Jaganathan et al., Cell 2019 (SpliceAI); precomputed scores
distributed by Ensembl. Please attribute Ensembl + the SpliceAI authors.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _cardiac_panel import panel_intervals  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
OUT_VCF = DATA_DIR / "spliceai_cardiac.masked.grch38.vcf.gz"

DEFAULT_SOURCE = (
    "https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/"
    "variation_plugins/"
    "spliceai_scores.masked.snv.ensembl_mane_v1.0.grch38.vcf.gz"
)

_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##INFO=<ID=SpliceAI,Number=.,Type=String,Description="
    '"SpliceAIv1.3 variant annotation. These include delta scores (DS) and '
    "delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor "
    "gain (DG), and donor loss (DL). Format: "
    'ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _open_source(source: str):
    """Open the source SpliceAI VCF (local path or remote URL) as a
    pysam.TabixFile. htslib transparently does remote range requests against
    the published .tbi for an http(s):// URL, so the 27 GB file never lands
    on disk."""
    import pysam

    return pysam.TabixFile(source)


def _contig_variants_for(tbx, chrom: str, start: int, end: int) -> list[str]:
    """Fetch the raw VCF lines overlapping [start, end] (1-based inclusive)
    for ``chrom``, trying both the bare ("3") and chr-prefixed ("chr3")
    contig spellings so the slice works regardless of the source naming.

    pysam.fetch is 0-based half-open; the panel BED carries 1-based-inclusive
    spans, so we query [start-1, end)."""
    contigs = set(tbx.contigs)
    candidates = []
    if chrom in contigs:
        candidates.append(chrom)
    chr_pref = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
    if chr_pref in contigs and chr_pref not in candidates:
        candidates.append(chr_pref)
    bare = chrom[3:] if str(chrom).startswith("chr") else chrom
    if bare in contigs and bare not in candidates:
        candidates.append(bare)

    lines: list[str] = []
    for c in candidates:
        try:
            for row in tbx.fetch(c, max(0, start - 1), end):
                lines.append(row)
        except (ValueError, OSError):
            continue
        if lines:
            break
    return lines


def _merge_spans(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Union overlapping / touching 1-based-inclusive spans, ascending.

    This is what removes the need to de-duplicate rows at all: the panel's
    padded intervals overlap (neighbouring genes, and pad=5000 on both sides),
    and the old builder absorbed that by keeping a ``seen`` set of every raw VCF
    line. Merge the spans instead and no position is ever fetched twice."""
    out: list[tuple[int, int]] = []
    for start, end in sorted(spans):
        if out and start <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], end))
        else:
            out.append((start, end))
    return out


def _source_spelling(tbx, chrom: str) -> str | None:
    """The contig spelling THIS source uses for ``chrom`` ("14" vs "chr14"),
    or None when the source does not carry it."""
    contigs = set(tbx.contigs)
    cands = [chrom, chrom if str(chrom).startswith("chr") else f"chr{chrom}"]
    if str(chrom).startswith("chr"):
        cands.append(chrom[3:])
    for c in cands:
        if c in contigs:
            return c
    return None


def build_slice(source: str, out_vcf: Path, pad: int) -> int:
    """Slice the source SpliceAI VCF to the cardiac panel and write a
    bgzipped + tabix-indexed VCF at ``out_vcf``. Returns records written.

    ⚠ MEMORY-SAFE, AND THAT IS THE WHOLE POINT OF THIS SHAPE. The previous
    version collected every raw VCF line into ``by_contig`` and kept a second
    copy in a ``seen`` de-dup set, on the stated assumption that "the panel is
    small (~200 genes) so the whole slice fits comfortably in memory". It does
    not: the correct slice is 62,003,508 records, so those two structures retain
    roughly 19 GB. Against the data-build job's ~2 GB container cap the build was
    OOM-killed — observed 2026-09-02 in run 33580187928 as
    ``Killed ... exit 137, 42s`` after only 42 seconds.

    That failure is also what caused the production outage this fixes. A killed
    build leaves a PARTIAL but structurally valid slice behind — the mount held a
    26 MB, 17-contig file dated 2026-07-21 carrying real SpliceAI INFO but no
    MYH7 or MYBPC3 records at all. Cadence ``static`` then kept it forever, and
    the client reported every missed position as an authoritative "no score", so
    no variant anywhere got a SpliceAI score and nothing said so.

    So this version holds NOTHING per record:
      * panel spans are merged per contig, which removes duplicate fetches and
        with them the entire de-dup set;
      * rows stream straight into a BGZF writer, so there is no 4.4 GB plain
        intermediate written onto the mount either;
      * merged spans are disjoint and processed in ascending order, and tabix
        returns rows position-ordered within a span, so the output is already
        sorted for indexing without an in-memory sort.

    Peak memory is one row at a time. The atomic temp-then-rename behaviour is
    unchanged, so a killed build still cannot publish a partial file — and the
    functional probe in build_all.sh's artifact_valid now rejects one that an
    earlier kill already published."""
    import pysam

    print(f"[spliceai] opening source {source}", flush=True)
    t0 = time.perf_counter()
    tbx = _open_source(source)

    intervals = panel_intervals(pad=pad)
    print(f"[spliceai] {len(intervals)} cardiac-panel intervals (pad={pad})",
          flush=True)

    by_src: dict[str, list[tuple[int, int]]] = {}
    unmapped: set[str] = set()
    for (chrom, start, end, _gene) in intervals:
        spelling = _source_spelling(tbx, chrom)
        if spelling is None:
            unmapped.add(str(chrom))
            continue
        by_src.setdefault(spelling, []).append((start, end))
    if unmapped:
        print(f"[spliceai] NOTE: source carries no contig for "
              f"{sorted(unmapped)} — skipped", flush=True)

    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    final_tbi = Path(str(out_vcf) + ".tbi")
    tmp_vcf = out_vcf.with_name(out_vcf.name + ".tmp")
    tmp_tbi = Path(str(tmp_vcf) + ".tbi")
    for stale in (tmp_vcf, tmp_tbi):
        if stale.exists():
            stale.unlink()

    total = 0
    spans_done = 0
    n_spans = sum(len(_merge_spans(v)) for v in by_src.values())
    writer = pysam.BGZFile(str(tmp_vcf), "wb")
    try:
        writer.write(_VCF_HEADER.encode("utf-8"))
        for src_contig in sorted(by_src):
            for (start, end) in _merge_spans(by_src[src_contig]):
                try:
                    rows = tbx.fetch(src_contig, max(0, start - 1), end)
                except (ValueError, OSError) as e:
                    print(f"[spliceai] WARN {src_contig}:{start}-{end} fetch "
                          f"failed ({e!r}) — skipped", flush=True)
                    continue
                for row in rows:
                    writer.write(
                        (row if row.endswith("\n") else row + "\n").encode("utf-8"))
                    total += 1
                spans_done += 1
                if spans_done % 25 == 0:
                    print(f"[spliceai] {spans_done}/{n_spans} merged spans, "
                          f"{total:,} records ({time.perf_counter() - t0:.0f} s)",
                          flush=True)
    finally:
        writer.close()
        tbx.close()

    print(f"[spliceai] wrote {total:,} records "
          f"({time.perf_counter() - t0:.1f} s); tabix indexing…", flush=True)
    if total == 0:
        raise SystemExit(
            "[spliceai] ABORT: zero records written — refusing to publish an "
            f"empty slice (temp left at {tmp_vcf})")

    pysam.tabix_index(str(tmp_vcf), preset="vcf", force=True)
    os.replace(str(tmp_tbi), str(final_tbi))
    os.replace(str(tmp_vcf), str(out_vcf))

    size_mb = out_vcf.stat().st_size / (1024 * 1024)
    print(f"[spliceai] ✓ {out_vcf} ({size_mb:.2f} MB) + .tbi", flush=True)
    return total


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        default=os.environ.get("SPLICEAI_SOURCE_VCF") or DEFAULT_SOURCE,
        help="source SpliceAI VCF (local path or remote http(s)/ftp URL with a "
             "sibling .tbi). Defaults to the Ensembl masked-SNV release.",
    )
    parser.add_argument(
        "--out",
        default=str(OUT_VCF),
        help=f"output bgzipped VCF path (default {OUT_VCF}).",
    )
    parser.add_argument(
        "--pad",
        type=int,
        default=None,
        help="bp padding each side of every panel gene span (default: the "
             "shared panel default).",
    )
    parser.add_argument("--force-download", action="store_true",
                        help="(accepted for parity; the slice is always rebuilt)")
    args = parser.parse_args(argv)

    out_vcf = Path(args.out)
    pad = args.pad if args.pad is not None else _default_pad()
    count = build_slice(args.source, out_vcf, pad)
    print(f"\nDone. {count:,} SpliceAI records sliced to {out_vcf}")
    return 0


def _default_pad() -> int:
    """The shared panel padding default (kept in sync with _cardiac_panel)."""
    from _cardiac_panel import DEFAULT_PAD
    return DEFAULT_PAD


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
