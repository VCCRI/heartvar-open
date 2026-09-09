#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:?usage: _build_revel_data.sh <out_dir> <assembly>}"
ASSEMBLY="${2:?usage: _build_revel_data.sh <out_dir> <assembly>}"
REVEL_ZIP_URL="${REVEL_ZIP_URL:-https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip}"

_asm_lc="$(printf '%s' "$ASSEMBLY" | tr '[:upper:]' '[:lower:]')"
OUT_TSV="$OUT_DIR/new_tabbed_revel_${_asm_lc}.tsv.gz"

if [[ "$ASSEMBLY" == "GRCh38" ]]; then pos_col=3; else pos_col=2; fi

_index_is_this_assembly() {
  local lo=23000000 hi=24000000 row v          # MYH7's neighbourhood, both builds
  row="$(tabix "$OUT_TSV" "14:${lo}-${hi}" 2>/dev/null | head -1)"
  [[ -n "$row" ]] || return 1
  v="$(printf '%s\n' "$row" | cut -f"$pos_col")"
  [[ "$v" =~ ^[0-9]+$ ]] || return 1
  (( v >= lo && v <= hi ))
}

_revel_usable() {
  [[ -s "$OUT_TSV" && -s "$OUT_TSV.tbi" ]] || return 1
  bgzip -t "$OUT_TSV" 2>/dev/null || return 1
  [[ -n "$(tabix -l "$OUT_TSV" 2>/dev/null)" ]] || return 1
  _index_is_this_assembly || return 1
  return 0
}

if [[ -z "${FORCE:-}${FORCE_REVEL:-}" ]]; then
  if _revel_usable; then
    echo "  REVEL data already built and verified at $OUT_TSV — skipping."
    exit 0
  fi
  if [[ -e "$OUT_TSV" ]]; then
    echo "  REVEL data at $OUT_TSV is present but UNUSABLE (truncated stream," >&2
    echo "  unreadable index, or indexed on the wrong assembly's position" >&2
    echo "  column) — rebuilding it rather than skipping." >&2
  fi
fi

for tool in curl unzip bgzip tabix awk sort tr head tail; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "ERROR: '$tool' not on PATH — needed to build the REVEL file." >&2
    exit 1
  }
done

mkdir -p "$OUT_DIR"

SORT_PARENT="${REVEL_SORT_TMP:-$OUT_DIR}"
mkdir -p "$SORT_PARENT"
SORT_TMP="$(mktemp -d "$SORT_PARENT/.revelsort.XXXXXX")"
tmp="$(mktemp -d "$OUT_DIR/.build.XXXXXX")"
trap 'rm -rf "$tmp" "$SORT_TMP"' EXIT

echo "  downloading $REVEL_ZIP_URL"
curl -fSL --retry 3 --retry-delay 10 "$REVEL_ZIP_URL" -o "$tmp/revel.zip"

set +o pipefail
hdr="$(unzip -p "$tmp/revel.zip" revel_with_transcript_ids | head -1 | tr ',' '\t')"
set -o pipefail
[[ -n "$hdr" ]] || { echo "ERROR: could not read the REVEL header." >&2; exit 1; }

echo "  building $OUT_TSV (assembly $ASSEMBLY, position column $pos_col)"
{
  printf '#%s\n' "$hdr"
  unzip -p "$tmp/revel.zip" revel_with_transcript_ids \
    | tail -n +2 \
    | tr ',' '\t' \
    | awk -F'\t' -v c="$pos_col" '$c != "."' \
    | sort -T "$SORT_TMP" -S "${REVEL_SORT_MEM:-1G}" \
           -k1,1 -k"${pos_col},${pos_col}n"
} | bgzip -c > "$tmp/out.tsv.gz"

tabix -f -s 1 -b "$pos_col" -e "$pos_col" "$tmp/out.tsv.gz"

bgzip -t "$tmp/out.tsv.gz"
mv -f "$tmp/out.tsv.gz" "$OUT_TSV"
mv -f "$tmp/out.tsv.gz.tbi" "$OUT_TSV.tbi"

echo "  ✓ REVEL → $OUT_TSV ($(du -h "$OUT_TSV" | cut -f1))"
