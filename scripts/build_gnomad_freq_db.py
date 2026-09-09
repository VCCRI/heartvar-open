#!/usr/bin/env python3
"""build_gnomad_freq_db.py — build the local gnomAD variant-frequency cache.

Tabix-slices the gnomAD v4.1 *sites* VCFs (exomes + genomes) to the cardiac
panel intervals and writes one row per variant into ``data/gnomad_freq.db``,
where each row's ``payload`` is the EXACT ``variant`` blob the gnomAD GraphQL
``variant`` query returns (``{variantId, rsid, exome{…}, genome{…}}``). The
backend gnomAD client serves variant frequencies from this DB in place of the
live GraphQL ``variant`` query (see ``backend/clients/gnomad.py``).

BUILD-AT-DEPLOY: the gnomAD sites VCFs are multi-TB across all chromosomes;
this script never downloads them whole — it tabix-fetches only the panel
regions (~200 cardiac genes, padded), reading remotely over HTTPS or from a
local copy. Run at deploy time once.

Source VCFs (per chromosome; GRCh38). Either the gs:// bucket, the AWS/Azure
mirrors, or the https form below:

  exomes:  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/
             vcf/exomes/gnomad.exomes.v4.1.sites.chr{N}.vcf.bgz
  genomes: https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/
             vcf/genomes/gnomad.genomes.v4.1.sites.chr{N}.vcf.bgz

  gs:// equivalents:
    gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr{N}.vcf.bgz
    gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr{N}.vcf.bgz

  AWS mirror (no-auth, often faster for tabix range reads):
    https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/exomes/…
    https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/…

pysam.TabixFile can open a remote .bgz URL directly and issue HTTP range reads
for each region — the matching ``.tbi`` index must sit alongside the .bgz
(gnomAD ships it). Pass ``--vcf-base`` to point at a local mirror directory
instead (expects the same filenames).

INFO → GraphQL mapping (per dataset = exome / genome), located BY INFO KEY:
  ac            <- AC
  an            <- AN
  af            <- AF
  ac_hom        <- nhomalt
  ac_hemi       <- AC_XY on chrX/chrY NON-PAR (there the XY allele count IS the
                   hemizygote count); 0 on autosomes + PAR. CRITICAL: on
                   autosomes AC_XY is just the male allele count (often large)
                   and is NOT a hemizygote count — the gnomAD GraphQL ac_hemi
                   field is 0 there. Verified against the live API:
                   X-100000137-A-C AC_XY=3 == GraphQL ac_hemi=3, while an
                   autosomal AC_XY=9208 maps to GraphQL ac_hemi=0. Mapping
                   AC_XY straight through would corrupt BS2 on autosomes.
  filters       <- FILTER column ([] when "." or "PASS" — the GraphQL variant
                   query returns an EMPTY list for a passing variant, not
                   ["PASS"]; non-PASS filters are split into a list).
  populations   <- [{id, ac, an, ac_hom, ac_hemi}] from AC_<grp>/AN_<grp>/
                   nhomalt_<grp>/AC_XY_<grp> for each of the 10 v4 genetic-
                   ancestry groups present (afr/amr/asj/eas/fin/nfe/sas/mid/
                   ami/remaining); per-group ac_hemi follows the same
                   sex-chrom rule. The live GraphQL also returns sex-split and
                   HGDP/1kg sub-cohorts, but the downstream curate pipeline
                   reads only the dataset-level af/ac_hom/ac_hemi/faf95 — never
                   iterates populations[] — so the 10-group list is a faithful,
                   sufficient subset.
  faf95.popmax            <- fafmax_faf95_max (the across-group max FAF95)
  faf95.popmax_population  <- fafmax_faf95_max_gen_anc (the group that maxes it)

Re-running is safe — the table is dropped and rebuilt every time.

Usage::

    python3 build_gnomad_freq_db.py                       # remote https slice
    python3 build_gnomad_freq_db.py --vcf-base /mnt/gnomad/4.1/vcf
    python3 build_gnomad_freq_db.py --chroms 1 2 7 X      # subset of chroms

Data source: gnomAD v4.1 (https://gnomad.broadinstitute.org), CC0. Please
attribute "gnomAD v4.1".
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import _dbbuild

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _cardiac_panel import panel_intervals  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DB_PATH = DATA_DIR / "gnomad_freq.db"

HTTPS_BASE = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf"
)
EXOME_TMPL = "{base}/exomes/gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz"
GENOME_TMPL = "{base}/genomes/gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz"

ANCESTRY_GROUPS = (
    "afr", "amr", "asj", "eas", "fin", "nfe", "sas", "mid", "ami", "remaining",
)

DDL = """
CREATE TABLE gnomad_freq (
    variant_id  TEXT PRIMARY KEY,
    payload     TEXT
)
"""

INSERT_SQL = "INSERT OR REPLACE INTO gnomad_freq VALUES (?, ?)"


def _info_to_dict(info: str) -> dict[str, str]:
    """Split a VCF INFO string ('AC=1;AN=2;AF=0.5;...') into a flat dict.
    Flag fields (no '=') map to '' (present)."""
    out: dict[str, str] = {}
    if not info or info == ".":
        return out
    for field in info.split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        else:
            out[field] = ""
    return out


def _maybe_int(d: dict[str, str], key: str) -> int | None:
    v = d.get(key)
    if v is None or v == "" or v == ".":
        return None
    v = v.split(",", 1)[0]
    try:
        return int(v)
    except ValueError:
        try:
            return int(float(v))
        except ValueError:
            return None


def _maybe_float(d: dict[str, str], key: str) -> float | None:
    v = d.get(key)
    if v is None or v == "" or v == ".":
        return None
    v = v.split(",", 1)[0]
    try:
        return float(v)
    except ValueError:
        return None


def _is_sex_chrom_nonpar(chrom: str, pos: int) -> bool:
    """True on the non-PAR regions of chrX/chrY, where AC_XY is the hemizygote
    count (== GraphQL ac_hemi). GRCh38 X PAR boundaries are excluded so PAR
    variants — diploid in both sexes — don't report a spurious hemi count.
    All of chrY is treated as hemizygous (gnomAD only reports XY there)."""
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if c == "Y":
        return True
    if c != "X":
        return False
    in_par1 = 10001 <= pos <= 2781479
    in_par2 = 155701383 <= pos <= 156030895
    return not (in_par1 or in_par2)


def _build_dataset_block(info: dict[str, str], filters: list[str],
                         hemi_applies: bool) -> dict | None:
    """Build one GraphQL dataset block (exome or genome) from an INFO dict.

    ``hemi_applies`` is True only on chrX/chrY non-PAR, where AC_XY is the
    hemizygote count; elsewhere ac_hemi is 0 (matching the GraphQL field).

    Returns None when the variant has no allele-number for this dataset
    (AN absent) — gnomAD's GraphQL returns null for a dataset the variant
    isn't observed in, so we mirror that (the merge layer drops a None block).
    """
    an = _maybe_int(info, "AN")
    ac = _maybe_int(info, "AC")
    if an is None and ac is None:
        return None

    populations: list[dict] = []
    for grp in ANCESTRY_GROUPS:
        grp_ac = _maybe_int(info, f"AC_{grp}")
        grp_an = _maybe_int(info, f"AN_{grp}")
        if grp_ac is None and grp_an is None:
            continue
        populations.append({
            "id": grp,
            "ac": grp_ac,
            "an": grp_an,
            "ac_hom": _maybe_int(info, f"nhomalt_{grp}"),
            "ac_hemi": (_maybe_int(info, f"AC_XY_{grp}") if hemi_applies else 0),
        })

    af = _maybe_float(info, "AF")
    ac_hom = _maybe_int(info, "nhomalt")
    ac_hemi = _maybe_int(info, "AC_XY") if hemi_applies else 0

    faf95_val = _maybe_float(info, "fafmax_faf95_max")
    if faf95_val is None:
        faf95_val = _maybe_float(info, "faf95_max")
    faf95_pop = info.get("fafmax_faf95_max_gen_anc") or info.get("faf95_max_gen_anc")
    if faf95_pop in ("", "."):
        faf95_pop = None

    return {
        "ac": ac,
        "an": an,
        "af": af,
        "ac_hom": ac_hom,
        "ac_hemi": ac_hemi,
        "filters": filters,
        "populations": populations,
        "faf95": {"popmax": faf95_val, "popmax_population": faf95_pop},
    }


def _parse_filter_col(filt: str) -> list[str]:
    """VCF FILTER column → GraphQL filters list. The gnomAD GraphQL variant
    query returns an EMPTY list for a passing variant (FILTER '.'/'PASS'), so
    we mirror that; a non-PASS FILTER ('AC0;InbreedingCoeff') splits into a
    list of the individual filter flags."""
    f = (filt or "").strip()
    if f in ("", ".", "PASS"):
        return []
    return f.split(";")


def vcf_record_to_payload_fields(
    chrom: str, pos: int, vid: str, ref: str, alt: str, filt: str, info_str: str,
) -> tuple[str, dict]:
    """Transform a single split VCF record (one ALT) into
    (variant_id, dataset_block) — the per-dataset half of the payload.

    ``variant_id`` is gnomAD's "chrom-pos-ref-alt" with no chr prefix, matching
    the client's lookup key and the GraphQL variantId form. ``vid`` is the VCF
    ID column (rsID, or '.').
    """
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    variant_id = f"{c}-{pos}-{ref}-{alt}"
    info = _info_to_dict(info_str)
    filters = _parse_filter_col(filt)
    hemi_applies = _is_sex_chrom_nonpar(chrom, pos)
    block = _build_dataset_block(info, filters, hemi_applies)
    rsid = None if vid in (".", "") else vid
    return variant_id, {"rsid": rsid, "block": block}


def merge_payload(variant_id: str, exome_half: dict | None,
                  genome_half: dict | None) -> dict:
    """Assemble the final GraphQL ``variant`` payload from the exome and/or
    genome halves produced by :func:`vcf_record_to_payload_fields`."""
    rsid = None
    if exome_half and exome_half.get("rsid"):
        rsid = exome_half["rsid"]
    elif genome_half and genome_half.get("rsid"):
        rsid = genome_half["rsid"]
    return {
        "variantId": variant_id,
        "rsid": rsid,
        "exome": exome_half["block"] if exome_half else None,
        "genome": genome_half["block"] if genome_half else None,
    }


def _open_tabix(url: str):
    """Open a (possibly remote) bgzipped VCF as a TabixFile, or None if it
    can't be opened (missing chrom file / network). Never raises."""
    import pysam

    try:
        return pysam.TabixFile(url)
    except (OSError, ValueError) as e:
        print(f"  ! could not open {url}: {e}", file=sys.stderr)
        return None


def _slice_one(tbx, chrom: str, start: int, end: int):
    """Yield raw VCF lines for [start, end] on ``chrom`` from an open
    TabixFile, trying both 'chr'-prefixed and bare contig names. tabix is
    0-based half-open; panel intervals are 1-based inclusive, so query
    [start-1, end)."""
    if tbx is None:
        return
    contigs = set(tbx.contigs)
    candidates = [c for c in (chrom, f"chr{chrom}",
                              chrom[3:] if chrom.startswith("chr") else None)
                  if c and c in contigs]
    for contig in candidates:
        try:
            for line in tbx.fetch(contig, max(0, start - 1), end):
                yield line
        except (ValueError, OSError):
            continue
        break


def build_db(vcf_base: str, chroms: list[str], db_path: Path) -> int:
    """Slice the gnomAD sites VCFs over the panel and write the freq DB.

    Returns the number of variant rows written. Variants observed in both the
    exome and genome VCFs are merged into a single payload with both blocks.
    """
    intervals = panel_intervals()
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for c, s, e, _g in intervals:
        cc = c[3:] if c.lower().startswith("chr") else c
        by_chrom.setdefault(cc, []).append((s, e))
    target_chroms = [c for c in chroms if c in by_chrom] if chroms else sorted(by_chrom)
    print(f"Panel covers {len(intervals)} intervals across "
          f"{len(by_chrom)} chromosomes; slicing {len(target_chroms)}.")

    tmp_path = _dbbuild.staging_db_path(db_path)
    conn = _dbbuild.connect(tmp_path)
    conn.execute("PRAGMA journal_mode = OFF")
    conn.execute("PRAGMA synchronous = OFF")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute(DDL)
    cur = conn.cursor()

    written = 0
    t0 = time.perf_counter()
    for chrom in target_chroms:
        regions = by_chrom[chrom]
        exome_url = EXOME_TMPL.format(base=vcf_base, chrom=chrom)
        genome_url = GENOME_TMPL.format(base=vcf_base, chrom=chrom)
        print(f"  chr{chrom}: {len(regions)} region(s)")
        tbx_ex = _open_tabix(exome_url)
        tbx_ge = _open_tabix(genome_url)
        try:
            halves: dict[str, dict] = {}
            for (start, end) in regions:
                for line in _slice_one(tbx_ex, chrom, start, end):
                    _accumulate(line, halves, "exome")
                for line in _slice_one(tbx_ge, chrom, start, end):
                    _accumulate(line, halves, "genome")
            batch = []
            for variant_id, h in halves.items():
                payload = merge_payload(variant_id, h.get("exome"), h.get("genome"))
                batch.append((variant_id, json.dumps(payload, separators=(",", ":"))))
                if len(batch) >= 20_000:
                    cur.executemany(INSERT_SQL, batch)
                    written += len(batch)
                    batch.clear()
            if batch:
                cur.executemany(INSERT_SQL, batch)
                written += len(batch)
            conn.commit()
            print(f"    ✓ chr{chrom}: {len(halves):,} variants "
                  f"({time.perf_counter() - t0:.1f} s elapsed)")
        finally:
            for t in (tbx_ex, tbx_ge):
                if t is not None:
                    try:
                        t.close()
                    except Exception:
                        pass

    print(f"  ✓ {written:,} variant rows in {time.perf_counter() - t0:.1f} s")
    conn.close()
    _dbbuild.publish(tmp_path, db_path)
    return written


def _accumulate(line: str, halves: dict[str, dict], dataset: str) -> None:
    """Parse a VCF line and stash its dataset half under its variant_id.
    Skips header lines and malformed records. Same-position multi-ALT lines
    are gnomAD-split (one ALT per line); we never split commas ourselves."""
    if not line or line.startswith("#"):
        return
    cols = line.split("\t")
    if len(cols) < 8:
        return
    chrom, pos_s, vid, ref, alt, _qual, filt, info_str = cols[:8]
    try:
        pos = int(pos_s)
    except ValueError:
        return
    if "," in alt:
        return
    variant_id, half = vcf_record_to_payload_fields(
        chrom, pos, vid, ref, alt, filt, info_str
    )
    if half["block"] is None:
        return
    slot = halves.setdefault(variant_id, {})
    slot.setdefault(dataset, half)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--vcf-base", default=HTTPS_BASE,
        help="base URL or local dir holding exomes/ and genomes/ subdirs "
             f"(default: {HTTPS_BASE})",
    )
    parser.add_argument(
        "--chroms", nargs="*", default=None,
        help="subset of chromosomes to slice (e.g. 1 2 7 X); default = all "
             "chromosomes the panel covers",
    )
    args = parser.parse_args(argv)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    count = build_db(args.vcf_base, args.chroms or [], DB_PATH)
    db_size_mb = DB_PATH.stat().st_size / (1024 * 1024)
    print(f"\nDone. {count:,} variant rows written to {DB_PATH} ({db_size_mb:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
