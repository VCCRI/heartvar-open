#!/usr/bin/env bash
set -euo pipefail

FORCE="${FORCE:-${FORCE_VEP:-}}"

VEP_ASSEMBLY="${VEP_ASSEMBLY:-GRCh38}"
VEP_CACHE_SPECIES="${VEP_CACHE_SPECIES:-homo_sapiens_merged}"
VEP_FASTA_SPECIES="${VEP_FASTA_SPECIES:-homo_sapiens}"
VEP_CACHE_FLAVOUR_FLAG="${VEP_CACHE_FLAVOUR_FLAG:---merged}"
VEP_BINARY="${VEP_BINARY:-vep}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VEP_DATA="${VEP_DATA:-$REPO_ROOT/data/vep}"
PLUGINS_DIR="$VEP_DATA/Plugins"
VEP_STATE_DIR="$VEP_DATA/.heartvar_vep_state"
VEP_MANIFEST="$VEP_DATA/.heartvar_vep_manifest.json"

REVEL_ZIP_URL="${REVEL_ZIP_URL:-https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip}"
REVEL_DIR="$VEP_DATA/revel"
_asm_lc="$(printf '%s' "$VEP_ASSEMBLY" | tr '[:upper:]' '[:lower:]')"
REVEL_TSV="$REVEL_DIR/new_tabbed_revel_${_asm_lc}.tsv.gz"

_vep_release() {
  local v="" bin="" constants=""
  bin="$(command -v "$VEP_BINARY" 2>/dev/null || true)"
  if [[ -n "$bin" ]]; then
    constants="$(dirname "$bin")/modules/Bio/EnsEMBL/VEP/Constants.pm"
    if [[ -r "$constants" ]]; then
      v="$(sed -n 's/.*VEP_VERSION[[:space:]]*=[[:space:]]*\([0-9][0-9]*\).*/\1/p' \
            "$constants" | head -1)" || v=""
    fi
    if [[ -z "$v" ]]; then
      v="$("$VEP_BINARY" --help 2>&1 \
           | sed -n 's/.*ensembl-vep *: *\([0-9][0-9]*\).*/\1/p' | head -1)" || v=""
    fi
  fi
  printf '%s' "${v:-${VEP_RELEASE:-113}}"
}

_free_gb() {
  local d="$1"
  while [[ -n "$d" && "$d" != "/" && ! -d "$d" ]]; do d="$(dirname "$d")"; done
  df -P -k "$d" 2>/dev/null | awk 'NR==2 {print int($4/1024/1024)}'
}

_require_gb() {               # _require_gb <dir> <gb> <what>
  local avail; avail="$(_free_gb "$1")"
  if [ "${avail:-0}" -lt "$2" ]; then
    echo "!!! $3: only ${avail:-0} GB free at $1; need >= $2 GB. Skipping it." >&2
    return 1
  fi
  echo "    free: ${avail} GB at $1 (>= $2 GB required) — OK"
  return 0
}

PROBE_VCF='14\t23429278\t.\tC\tT\t.\t.\t.'
PROBE_AA="R/Q"                 # p.Arg403Gln
PROBE_HGVSC="c.1208G>A"        # on NM_000257 / ENST00000355349

_vep_smoke() {                # _vep_smoke [extra vep args...]
  local fasta_args=()
  if [[ "$*" != *--fasta* ]]; then
    if [[ -r "$FASTA_PATH" ]]; then
      fasta_args=(--fasta "$FASTA_PATH")
    else
      fasta_args=(--use_given_ref)
    fi
  fi
  printf '%b\n' "$PROBE_VCF" \
    | "$VEP_BINARY" --offline --cache --dir_cache "$VEP_DATA" \
        --assembly "$VEP_ASSEMBLY" --json --no_stats --force_overwrite \
        $VEP_CACHE_FLAVOUR_FLAG \
        --format vcf -o STDOUT -i STDIN \
        ${fasta_args[@]+"${fasta_args[@]}"} "$@" 2>/dev/null || true
}

_check_cache_info() {
  local info="$CACHE_DIR/info.txt" v
  if [[ ! -s "$info" ]]; then
    echo "!!! $info is MISSING or empty." >&2
    echo "    The cache tree can look complete without it — it is the LAST" >&2
    echo "    member of the tarball, so a stream that stopped anywhere before" >&2
    echo "    the end leaves every chromosome present and this file absent." >&2
    echo "    VEP reads it with a bare open() and no else branch, so the" >&2
    echo "    symptom is 'ERROR: SIFT not available', not a missing file." >&2
    return 1
  fi
  v="$(awk -F'\t' '$1=="assembly" {print $2; exit}' "$info")"
  if [[ -n "$v" && "$v" != "$VEP_ASSEMBLY" ]]; then
    echo "!!! $info says assembly '$v' but this install wants '$VEP_ASSEMBLY'." >&2
    return 1
  fi
  if ! _vep_smoke --sift b --polyphen b | grep -q '"sift_prediction"'; then
    echo "!!! $info exists, but --sift b --polyphen b produced no" >&2
    echo "    sift_prediction. clients/vep_offline.py::_argv passes both on" >&2
    echo "    EVERY call and VEP throws rather than degrading, so this is fatal" >&2
    echo "    to every curation." >&2
    echo "    info.txt declares:" >&2
    for v in sift polyphen source_sift source_polyphen; do
      printf '      %-16s %s\n' "$v" \
        "$(awk -F'\t' -v k="$v" '$1==k {print (NF>1 && $2!="" ? $2 : "(empty)"); exit}' "$info")" >&2
    done
    echo "    An empty or absent declaration is how a species without the data" >&2
    echo "    says so; for human both should be populated." >&2
    return 1
  fi
  return 0
}

_check_cache() {
  _check_cache_info || return 1
  local out; out="$(_vep_smoke)"
  printf '%s' "$out" | grep -q '"missense_variant"' || return 1
  printf '%s' "$out" | grep -q "$PROBE_AA" || return 1
  return 0
}

_check_fasta() { _vep_smoke --hgvs --fasta "$FASTA_PATH" | grep -q "$PROBE_HGVSC"; }

_check_revel() {
  _vep_smoke --dir_plugins "$PLUGINS_DIR" --plugin "REVEL,file=$REVEL_TSV,no_match=1" \
    | grep -qE '"(revel|REVEL|revel_score|REVEL_score)" *:'
}

PROBE_RUNTIME_VCF='14\t23429278\t.\tC\tT\t.\t.\t.'
PROBE_RUNTIME_RESIDUE="p.Arg403Gln"
PROBE_RUNTIME_ENSP="ENSP00000347507"
PROBE_RUNTIME_MERGED_TX_RE='"transcript_id" *: *"(NM_|NR_|XM_|XR_)'

_check_runtime_call() {
  local out plugin_args=()
  if [[ -z "${SKIP_REVEL:-}" && -s "$REVEL_TSV" ]]; then
    plugin_args=(--dir_plugins "$PLUGINS_DIR" --plugin "REVEL,file=$REVEL_TSV,no_match=1")
  fi
  out="$(printf '%b\n' "$PROBE_RUNTIME_VCF" \
    | "$VEP_BINARY" --offline --cache --dir_cache "$VEP_DATA" \
        --assembly "$VEP_ASSEMBLY" --json --no_stats --force_overwrite \
        -o STDOUT -i STDIN --format vcf \
        --hgvs --symbol --mane --numbers --canonical --biotype \
        --sift b --polyphen b \
        --fasta "$FASTA_PATH" $VEP_CACHE_FLAVOUR_FLAG \
        ${plugin_args[@]+"${plugin_args[@]}"} \
        2>/dev/null || true)"
  printf '%s' "$out" | grep -q '"missense_variant"' || return 1
  printf '%s' "$out" | grep -q "$PROBE_RUNTIME_RESIDUE" || return 1
  printf '%s' "$out" | grep -q "$PROBE_RUNTIME_ENSP" || return 1
  printf '%s' "$out" | grep -qE "$PROBE_RUNTIME_MERGED_TX_RE" || return 1
  return 0
}

_runtime_diagnostic() {
  local plugin_args=()
  if [[ -z "${SKIP_REVEL:-}" && -s "$REVEL_TSV" ]]; then
    plugin_args=(--dir_plugins "$PLUGINS_DIR" --plugin "REVEL,file=$REVEL_TSV,no_match=1")
  fi
  echo "--- RUNTIME CALL diagnostic (vep stderr, HEAD first) ---" >&2
  printf '%b\n' "$PROBE_RUNTIME_VCF" \
    | "$VEP_BINARY" --offline --cache --dir_cache "$VEP_DATA" \
        --assembly "$VEP_ASSEMBLY" --json --no_stats --force_overwrite \
        -o STDOUT -i STDIN --format vcf \
        --hgvs --symbol --mane --numbers --canonical --biotype \
        --sift b --polyphen b \
        --fasta "$FASTA_PATH" $VEP_CACHE_FLAVOUR_FLAG \
        ${plugin_args[@]+"${plugin_args[@]}"} \
        2>&1 >/dev/null | head -c 900 >&2 || true
  echo "" >&2

  echo "--- RUNTIME CALL diagnostic (vep STDOUT — compare against the four" >&2
  echo "    PROBE_RUNTIME_* assertions above; these were written against a stub) ---" >&2
  printf '%b\n' "$PROBE_RUNTIME_VCF" \
    | "$VEP_BINARY" --offline --cache --dir_cache "$VEP_DATA" \
        --assembly "$VEP_ASSEMBLY" --json --no_stats --force_overwrite \
        -o STDOUT -i STDIN --format vcf \
        --hgvs --symbol --mane --numbers --canonical --biotype \
        --sift b --polyphen b \
        --fasta "$FASTA_PATH" $VEP_CACHE_FLAVOUR_FLAG \
        ${plugin_args[@]+"${plugin_args[@]}"} \
        2>/dev/null | head -c 4000 >&2 || true
  echo "" >&2
  local _out
  _out="$(printf '%b\n' "$PROBE_RUNTIME_VCF" \
    | "$VEP_BINARY" --offline --cache --dir_cache "$VEP_DATA" \
        --assembly "$VEP_ASSEMBLY" --json --no_stats --force_overwrite \
        -o STDOUT -i STDIN --format vcf \
        --hgvs --symbol --mane --numbers --canonical --biotype \
        --sift b --polyphen b \
        --fasta "$FASTA_PATH" $VEP_CACHE_FLAVOUR_FLAG \
        ${plugin_args[@]+"${plugin_args[@]}"} 2>/dev/null || true)"
  for _want in '"missense_variant"' "$PROBE_RUNTIME_RESIDUE" "$PROBE_RUNTIME_ENSP"; do
    printf '    assertion  %-34s : %s\n' "$_want" \
      "$(printf '%s' "$_out" | grep -qF "$_want" && echo FOUND || echo MISSING)" >&2
  done
  printf '    assertion  %-34s : %s\n' "RefSeq transcript_id (merged)" \
    "$(printf '%s' "$_out" | grep -qE "$PROBE_RUNTIME_MERGED_TX_RE" \
        && echo FOUND || echo "MISSING  <- cache is Ensembl-only, not merged")" >&2
  echo "    cache flavour flag : $VEP_CACHE_FLAVOUR_FLAG" >&2
  echo "    cache dir expected : $CACHE_DIR" >&2
  if [[ -s "$CACHE_DIR/info.txt" ]]; then
    echo "    info.txt           : present ($(wc -l < "$CACHE_DIR/info.txt") lines)" >&2
    echo "      source_sift      : $(awk -F'\t' '$1=="source_sift" {print $2; exit}' "$CACHE_DIR/info.txt")" >&2
    echo "      source_polyphen  : $(awk -F'\t' '$1=="source_polyphen" {print $2; exit}' "$CACHE_DIR/info.txt")" >&2
    echo "      assembly         : $(awk -F'\t' '$1=="assembly" {print $2; exit}' "$CACHE_DIR/info.txt")" >&2
  else
    echo "    info.txt           : ⚠ MISSING — THIS IS ALMOST CERTAINLY THE CAUSE." >&2
    echo "      'ERROR: SIFT not available' is what VEP says when info.txt could" >&2
    echo "      not be opened: read_cache_info_file returns an empty hash with no" >&2
    echo "      else branch, so check_sift_polyphen finds no sift and throws." >&2
    echo "      info.txt is the LAST member of the tarball, so the chromosome" >&2
    echo "      directories being present says nothing about it. Re-run with" >&2
    echo "      FORCE_VEP=1 to re-stream the cache." >&2
  fi
  ls -d "$VEP_DATA"/*/ >&2 2>&1 || true
  echo "--- end RUNTIME CALL diagnostic ---" >&2
}

_fasta_diagnostic() {
  echo "--- FASTA diagnostic (vep stderr, failure path only) ---" >&2
  printf '%b\n' "$PROBE_VCF" \
    | "$VEP_BINARY" --offline --cache --dir_cache "$VEP_DATA" \
        --assembly "$VEP_ASSEMBLY" --json --no_stats --force_overwrite \
        $VEP_CACHE_FLAVOUR_FLAG \
        --format vcf -o STDOUT -i STDIN --hgvs --fasta "$FASTA_PATH" \
        2>&1 >/dev/null | head -40 >&2 || true
  ls -l "$FASTA_PATH" "$FASTA_PATH.fai" "$FASTA_PATH.gzi" >&2 2>&1 || true
  echo "--- end FASTA diagnostic ---" >&2
}

_revel_diagnostic() {
  echo "--- REVEL diagnostic (failure path only) ---" >&2
  ls -l "$REVEL_TSV" "$REVEL_TSV.tbi" >&2 2>&1 || true
  if bgzip -t "$REVEL_TSV" 2>/dev/null; then
    echo "    bgzip -t: OK (stream complete)" >&2
  else
    echo "    bgzip -t: FAILED — truncated or not BGZF. This is a bad build," >&2
    echo "    not a plugin problem; rebuild with FORCE_REVEL=1." >&2
  fi
  echo "    tabix -l: $(tabix -l "$REVEL_TSV" 2>/dev/null | head -5 | tr '\n' ' ')" >&2
  echo "    header  : $(tabix -H "$REVEL_TSV" 2>/dev/null | head -1)" >&2
  local probe_chrom probe_pos
  probe_chrom="$(printf '%b' "$PROBE_VCF" | cut -f1)"
  probe_pos="$(printf '%b' "$PROBE_VCF" | cut -f2)"
  echo "    probe row (${probe_chrom}:${probe_pos}):" >&2
  tabix "$REVEL_TSV" "${probe_chrom}:${probe_pos}-${probe_pos}" 2>&1 | head -5 >&2 || true
  echo "    first rows in ${probe_chrom}:23000000-24000000 (cols: chr,hg19,grch38,ref,alt):" >&2
  tabix "$REVEL_TSV" "${probe_chrom}:23000000-24000000" 2>&1 \
    | head -2 | cut -f1-5 >&2 || true
  ls -l "$PLUGINS_DIR/REVEL.pm" >&2 2>&1 || true

  echo "--- REVEL probe with stderr attached ---" >&2
  printf '%b\n' "$PROBE_VCF" \
    | "$VEP_BINARY" --offline --cache --dir_cache "$VEP_DATA" \
        --assembly "$VEP_ASSEMBLY" --json --no_stats --force_overwrite \
        $VEP_CACHE_FLAVOUR_FLAG \
        --format vcf -o STDOUT -i STDIN \
        --dir_plugins "$PLUGINS_DIR" --plugin "REVEL,file=$REVEL_TSV,no_match=1" \
        2>&1 | head -c 4000 >&2 || true
  echo "" >&2
  echo "--- end REVEL diagnostic ---" >&2
}

_state_ok() {                 # _state_ok <component>
  local f="$VEP_STATE_DIR/$1.json"
  [[ -s "$f" ]] || return 1
  grep -q "\"vep_release\": *\"$REL\"" "$f" || return 1
  grep -q "\"assembly\": *\"$VEP_ASSEMBLY\"" "$f" || return 1
  return 0
}

_state_write() {              # _state_write <component>
  mkdir -p "$VEP_STATE_DIR"
  printf '{\n  "component": "%s",\n  "vep_release": "%s",\n  "assembly": "%s",\n  "completed_utc": "%s"\n}\n' \
    "$1" "$REL" "$VEP_ASSEMBLY" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
    > "$VEP_STATE_DIR/$1.json"
}

_manifest_write() {
  printf '{\n  "vep_release": "%s",\n  "assembly": "%s",\n  "species": "%s",\n  "cache_dir": "%s",\n  "fasta": "%s",\n  "revel": "%s",\n  "plugins_dir": "%s",\n  "completed_utc": "%s"\n}\n' \
    "$REL" "$VEP_ASSEMBLY" "$VEP_CACHE_SPECIES" "$CACHE_DIR" "$FASTA_PATH" \
    "$REVEL_TSV" "$PLUGINS_DIR" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
    > "$VEP_MANIFEST"
  echo "    wrote $VEP_MANIFEST"
}

REL="$(_vep_release)"
CACHE_DIR="$VEP_DATA/$VEP_CACHE_SPECIES/${REL}_${VEP_ASSEMBLY}"
FASTA_DIR="$VEP_DATA/$VEP_FASTA_SPECIES/${REL}_${VEP_ASSEMBLY}"
FASTA_PATH="$FASTA_DIR/Homo_sapiens.${VEP_ASSEMBLY}.dna.toplevel.fa.gz"

CACHE_URL="${CACHE_URL:-https://ftp.ensembl.org/pub/release-${REL}/variation/indexed_vep_cache/${VEP_CACHE_SPECIES}_vep_${REL}_${VEP_ASSEMBLY}.tar.gz}"

_EXCLUDE_ALL_VARS=(--exclude='all_vars.gz*' --exclude='*_var.gz')
if [[ -n "${VEP_CACHE_KEEP_VARIATION:-}" ]]; then
  _EXCLUDE_ALL_VARS=()
fi
if [[ -n "${VEP_CACHE_KEEP_VARIATION:-}" ]]; then
  CACHE_FLOOR_GB="${CACHE_FLOOR_GB:-28}"
else
  CACHE_FLOOR_GB="${CACHE_FLOOR_GB:-8}"
fi

_install_cache() {
  echo "==> Cache: streaming $CACHE_URL"
  if [[ ${#_EXCLUDE_ALL_VARS[@]} -gt 0 ]]; then
    echo "    excluding the variation cache (~81% of the tree; never read"
    echo "    without --check_existing — see CacheDir.pm:167)"
  else
    echo "    VEP_CACHE_KEEP_VARIATION set — installing the FULL tree (~26 GB)"
  fi
  mkdir -p "$VEP_DATA"
  curl -fSL --retry 3 --retry-delay 10 "$CACHE_URL" \
    | tar -xzf - -C "$VEP_DATA" ${_EXCLUDE_ALL_VARS[@]+"${_EXCLUDE_ALL_VARS[@]}"}
}

_install_fasta() {
  echo "==> FASTA: INSTALL.pl --AUTO f"
  INSTALL.pl --AUTO f --SPECIES "$VEP_FASTA_SPECIES" --ASSEMBLY "$VEP_ASSEMBLY" \
    --CACHEDIR "$VEP_DATA" --NO_UPDATE
}

_install_revel_plugin() {
  echo "==> REVEL plugin module: INSTALL.pl --AUTO p"
  INSTALL.pl --AUTO p --PLUGINS REVEL --PLUGINSDIR "$PLUGINS_DIR" --NO_UPDATE
}

if [[ -z "${CHECK_ONLY:-}" ]] && ! command -v INSTALL.pl >/dev/null 2>&1; then
  echo "ERROR: INSTALL.pl is not on PATH. Run this inside an Ensembl VEP image —" >&2
  echo "       the data-builder image is built from one (see Dockerfile.builder)." >&2
  exit 1
fi

echo "==> Offline VEP setup"
echo "    release  : $REL (from $VEP_BINARY)"
echo "    assembly : $VEP_ASSEMBLY ($VEP_CACHE_SPECIES, flavour flag $VEP_CACHE_FLAVOUR_FLAG)"
echo "    data dir : $VEP_DATA"
if [[ -n "${SKIP_REVEL:-}" ]]; then echo "    REVEL    : skipped"; else echo "    REVEL    : $REVEL_TSV"; fi

if [[ -n "${CHECK_ONLY:-}" ]]; then
  echo "==> CHECK ONLY — probing the existing install. Nothing will be downloaded."
  MISSING=0
  for comp in cache fasta revel runtime; do
    if [[ "$comp" == revel && -n "${SKIP_REVEL:-}" ]]; then
      echo "    revel : skipped (SKIP_REVEL set)"
      continue
    fi
    state="ABSENT"
    if _state_ok "$comp"; then state="recorded"; fi
    if [[ "$comp" == runtime ]]; then state="n/a"; fi
    works="NO"
    case "$comp" in
      cache) if _check_cache; then works="yes"; fi ;;
      fasta) if _check_fasta; then works="yes"; fi ;;
      revel) if _check_revel; then works="yes"; fi ;;
      runtime) if _check_runtime_call; then works="yes"; fi ;;
    esac
    printf '    %-5s : state=%-8s functional=%s\n' "$comp" "$state" "$works"
    if [[ "$works" != yes ]]; then
      MISSING=$((MISSING + 1))
      if [[ "$comp" == fasta ]]; then _fasta_diagnostic; fi
      if [[ "$comp" == revel ]]; then _revel_diagnostic; fi
      if [[ "$comp" == runtime ]]; then _runtime_diagnostic; fi
    fi
  done

  if [[ "$MISSING" -eq 0 ]]; then
    echo "==> Every requested component is installed and working for release $REL."
    if [[ -s "$VEP_MANIFEST" ]]; then
      echo "    Manifest present: $VEP_MANIFEST"
    else
      echo "    NO MANIFEST — components verify individually but the all-or-nothing"
      echo "    manifest was never written. A normal run writes it in minutes,"
      echo "    because every component below is already recorded and skipped."
    fi
    exit 0
  fi
  echo ""
  echo "!!! CHECK ONLY: $MISSING component(s) not usable. Nothing was installed." >&2
  echo "    'state=recorded functional=NO' means the mount changed under a" >&2
  echo "    completed component; 'state=ABSENT' means it never finished." >&2
  exit 1
fi

FAILURES=0

if [[ -z "${FORCE:-}" ]] && _state_ok cache && _check_cache; then
  echo "==> Cache already installed and working for release $REL — skipping."
elif _require_gb "$VEP_DATA" "$CACHE_FLOOR_GB" "VEP cache"; then
  _install_cache || echo "!!! Cache download/extract failed (curl/tar output above)." >&2
  if _check_cache; then
    _state_write cache
    echo "==> Cache OK (annotated MYH7 R403Q from $CACHE_DIR)."
  else
    echo "!!! Cache is NOT usable — the download failed, or vep cannot annotate" >&2
    echo "    with what landed. Not recording it, so the next build retries." >&2
    FAILURES=$((FAILURES + 1))
  fi
else
  FAILURES=$((FAILURES + 1))
fi

if [[ -z "${FORCE:-}" ]] && _state_ok fasta && _check_fasta; then
  echo "==> FASTA already installed for release $REL — skipping."
elif ! _state_ok cache; then
  echo "!!! Skipping FASTA: the cache is not usable yet." >&2
  FAILURES=$((FAILURES + 1))
elif _require_gb "$VEP_DATA" 3 "VEP FASTA"; then
  _install_fasta || echo "!!! FASTA install failed (INSTALL.pl output above)." >&2
  if _check_fasta; then
    _state_write fasta
    echo "==> FASTA OK (HGVS resolved using $FASTA_PATH)."
  else
    echo "!!! FASTA is NOT usable — install failed, or no HGVS came back. Not" >&2
    echo "    recording it, so the next build retries." >&2
    _fasta_diagnostic
    FAILURES=$((FAILURES + 1))
  fi
else
  FAILURES=$((FAILURES + 1))
fi

if [[ -n "${SKIP_REVEL:-}" ]]; then
  echo "==> REVEL: skipped (SKIP_REVEL set)."
elif [[ -z "${FORCE:-}" ]] && _state_ok revel && _check_revel; then
  echo "==> REVEL already built for release $REL — skipping."
elif _require_gb "$REVEL_DIR" 10 "REVEL"; then
  mkdir -p "$REVEL_DIR" "$PLUGINS_DIR"
  _install_revel_plugin || echo "!!! REVEL plugin module install failed." >&2
  REVEL_ZIP_URL="$REVEL_ZIP_URL" bash "$REPO_ROOT/scripts/_build_revel_data.sh" \
    "$REVEL_DIR" "$VEP_ASSEMBLY" \
    || echo "!!! REVEL data build failed (is REVEL_ZIP_URL still valid?)." >&2
  if _check_revel; then
    _state_write revel
    echo "==> REVEL OK — a REVEL score came through the CLI output."
    echo "    (that is the plugin field-name validation gating HEARTVAR_VEP_OFFLINE)"
  else
    echo "!!! REVEL is NOT usable — build failed, or no REVEL key came back." >&2
    echo "    PP3/BP4 would silently lose its calibrated signal. Compare the" >&2
    echo "    plugin's key casing against _first() in clients/vep_offline.py." >&2
    _revel_diagnostic
    FAILURES=$((FAILURES + 1))
  fi
else
  FAILURES=$((FAILURES + 1))
fi

if [[ "$FAILURES" -eq 0 ]]; then
  if _check_runtime_call; then
    echo "==> RUNTIME CALL OK — the web app's own invocation annotated"
    echo "    $(printf '%b' "$PROBE_RUNTIME_VCF" | tr '\t' ' ') -> $PROBE_RUNTIME_RESIDUE"
    echo "    from $CACHE_DIR,"
    echo "    with both transcript sets present (merged cache confirmed)."
  else
    echo "!!! THE RUNTIME CALL FAILS even though every component verified." >&2
    echo "    This is the 2026-08-26 failure mode: a verified cache is not a" >&2
    echo "    verified path. Do NOT set HEARTVAR_VEP_OFFLINE on this install." >&2
    _runtime_diagnostic
    FAILURES=$((FAILURES + 1))
  fi
fi

if [[ "$FAILURES" -eq 0 ]]; then
  _manifest_write
  echo ""
  echo "==> Offline VEP is installed and VERIFIED for release $REL."
  echo "    HEARTVAR_VEP_DATA=$VEP_DATA"
  echo "    HEARTVAR_VEP_ASSEMBLY=$VEP_ASSEMBLY"
  echo "    HEARTVAR_VEP_FASTA=$FASTA_PATH"
  if [[ -z "${SKIP_REVEL:-}" ]]; then
    echo "    HEARTVAR_VEP_REVEL=$REVEL_TSV"
    echo "    HEARTVAR_VEP_PLUGINS_DIR=$PLUGINS_DIR"
  fi
  echo "    The deployment sets HEARTVAR_VEP_OFFLINE=1, so on the Web App this"
  echo "    cache is the annotation path. Elsewhere it is inert until you set it."
  exit 0
fi

echo ""
echo "!!! Offline VEP setup finished with $FAILURES failed component(s)." >&2
echo "    No manifest written, so the next build resumes the unfinished parts." >&2
exit 1
