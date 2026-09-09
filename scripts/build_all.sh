#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$REPO_ROOT/data"
BACKEND_DATA="$REPO_ROOT/backend/data"
PY="${PYTHON:-python3}"

_HV_LOG_DIR="${DATA_DIR}/logs/db_builder"
mkdir -p "$_HV_LOG_DIR" 2>/dev/null || true
_HV_LOG_FILE="$_HV_LOG_DIR/$(date -u '+%Y-%m-%dT%H-%M-%SZ').log"
exec > >(tee "$_HV_LOG_FILE") 2>&1

SKIP_LIST="${SKIP:-}"
RESUME="${RESUME:-}"
MONTHLY="${MONTHLY:-}"
WITH_VEP="${WITH_VEP:-}"
ONLY_LIST="${ONLY:-}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip) SKIP_LIST="$2"; shift 2;;
    --skip=*) SKIP_LIST="${1#*=}"; shift;;
    --only) ONLY_LIST="$2"; shift 2;;
    --only=*) ONLY_LIST="${1#*=}"; shift;;
    --resume) RESUME=1; shift;;
    --monthly) MONTHLY=1; shift;;
    --with-vep) WITH_VEP=1; shift;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
if [[ -n "$RESUME" && -n "$MONTHLY" ]]; then
  echo "!!! NOTE  both --resume and --monthly given; --resume wins (skips all present)." >&2
fi

output_for() {
  case "$1" in
    uniprot)           echo "$DATA_DIR/uniprot.db";;
    clinvar)           echo "$DATA_DIR/clinvar.db";;
    gnomad_constraint) echo "$DATA_DIR/gnomad_constraint.db";;
    gtex)              echo "$DATA_DIR/gtex.db";;
    opentargets)       echo "$DATA_DIR/opentargets.db";;
    medgen)            echo "$DATA_DIR/medgen.db";;
    mgi)               echo "$DATA_DIR/mgi.db";;
    biogrid)           echo "$DATA_DIR/biogrid.db";;
    hpo_labels)        echo "$BACKEND_DATA/hpo_labels.json";;
    panelapp)          echo "$DATA_DIR/panelapp_aus_snapshot.json";;
    gnomad_freq)       echo "$DATA_DIR/gnomad_freq.db";;
    alphafold)         echo "$DATA_DIR/alphafold/manifest.json";;
    spliceai)          echo "$DATA_DIR/spliceai_cardiac.masked.grch38.vcf.gz";;
    fetal_heart)       echo "$DATA_DIR/fetal_heart.db";;
    hgnc_alias)        echo "$BACKEND_DATA/hgnc_alias_map.json";;
    gene_id_map)       echo "$BACKEND_DATA/gene_id_map.json.gz";;
    erepo)             echo "$BACKEND_DATA/erepo_all.tsv";;
    gencc)             echo "$BACKEND_DATA/gencc_submissions.json";;
    clingen_gv)        echo "$BACKEND_DATA/clingen_gene_validity.json";;
    vep)               echo "$DATA_DIR/vep/.heartvar_vep_manifest.json";;
    cdot)              echo "$DATA_DIR/cdot/cdot_transcripts.db";;
    *)                 echo "";;
  esac
}

cadence_for() {
  case "$1" in
    uniprot|clinvar|medgen|mgi|panelapp|hpo_labels|hgnc_alias|erepo|gencc|clingen_gv)
                       echo monthly;;
    opentargets|biogrid)
                       echo monthly;;
    gnomad_constraint|gnomad_freq|gtex|vep|cdot)
                       echo versioned;;
    alphafold|spliceai|fetal_heart)
                       echo static;;
    gene_id_map)       echo image;;
    *)                 echo monthly;;
  esac
}

artifact_valid() {
  local label="$1" out="$2"
  case "$label" in
    spliceai)
      if [[ ! -s "${out}.tbi" ]]; then
        echo "!!! $label: index ${out}.tbi is missing or empty — rebuilding"
        return 1
      fi
      if ! "$PY" - "$out" <<'PYEOF'; then
import os, sys, time
try:
    import pysam
except ImportError:                     # cannot verify -> do not force a rebuild
    sys.exit(0)

path = sys.argv[1]
# Anchors that MUST resolve in any correct cardiac slice: MYH7 and MYBPC3 are
# the two most central genes on the panel. Same trade-off cdot already accepts
# ("it resolves MYH7 R403Q before it exits 0") — two hardcoded coordinates buy
# a check that a bare "some row reads" cannot make.
#
# WHY THE WEAKER CHECK WAS NOT ENOUGH (2026-09-02): a slice on the mount passed
# .tbi + contigs + one-row-reads, and yet EVERY production lookup returned "No
# SpliceAI score at this position" — confirmed that no variant anywhere gets a
# score. A slice missing MYH7 and MYBPC3 entirely is not a coverage quirk, it is
# the wrong artifact, and only a functional probe catches that.
ANCHORS = [("14", 23433087, "MYH7"), ("11", 47331882, "MYBPC3"),
           ("14", 23423601, "MYH7")]

print(f"[spliceai:probe] path={path}")
try:
    st = os.stat(path)
    print(f"[spliceai:probe] size={st.st_size} bytes  "
          f"mtime={time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(st.st_mtime))}")
    tbi = os.stat(path + ".tbi")
    print(f"[spliceai:probe] tbi_size={tbi.st_size} bytes  "
          f"tbi_mtime={time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(tbi.st_mtime))}")
except OSError as e:
    print(f"!!! spliceai: cannot stat ({e!r})")
    sys.exit(1)

try:
    tbx = pysam.TabixFile(path)
    contigs = list(tbx.contigs)
except Exception as e:
    print(f"!!! spliceai: slice unreadable ({e!r})")
    sys.exit(1)
print(f"[spliceai:probe] contigs={len(contigs)} {sorted(contigs)[:26]}")
if not contigs:
    print("!!! spliceai: slice carries ZERO contigs — empty or truncated")
    tbx.close()
    sys.exit(1)

def summarise(row):
    """DERIVED facts only — never the row itself."""
    cols = row.split("\t")
    info = cols[7] if len(cols) > 7 else ""
    has_key = "SpliceAI=" in info
    nfields = 0
    if has_key:
        for field in info.split(";"):
            if field.startswith("SpliceAI="):
                nfields = len(field.split("=", 1)[1].split(",")[0].split("|"))
                break
    return {"cols": len(cols), "pos": cols[1] if len(cols) > 1 else "?",
            "spliceai_key": has_key, "annot_fields": nfields}

# What does the file actually contain? First readable row, as a summary.
first = None
try:
    for contig in contigs[:3]:
        for row in tbx.fetch(contig):
            first = (contig, summarise(row))
            break
        if first:
            break
except Exception as e:
    print(f"!!! spliceai: slice read failed ({e!r}) — truncated or corrupt")
    tbx.close()
    sys.exit(1)
if first is None:
    print(f"!!! spliceai: index names {len(contigs)} contigs but NO rows could "
          "be read — slice is truncated or corrupt")
    tbx.close()
    sys.exit(1)
print(f"[spliceai:probe] first_row contig={first[0]} {first[1]}")

# Functional probe: the anchors must yield SpliceAI-bearing rows.
ok_anchors = 0
for chrom, pos, gene in ANCHORS:
    hit = {"rows": 0, "with_spliceai": 0}
    for spelling in (chrom, "chr" + chrom):
        if spelling not in contigs:
            continue
        try:
            for row in tbx.fetch(spelling, pos - 1, pos):
                d = summarise(row)
                hit["rows"] += 1
                hit["with_spliceai"] += 1 if d["spliceai_key"] else 0
        except Exception as e:
            print(f"[spliceai:probe] {gene} {chrom}:{pos} fetch failed ({e!r})")
        break
    print(f"[spliceai:probe] anchor {gene} {chrom}:{pos} -> "
          f"rows={hit['rows']} with_spliceai_key={hit['with_spliceai']}")
    if hit["with_spliceai"]:
        ok_anchors += 1
tbx.close()

if ok_anchors == 0:
    print(f"!!! spliceai: NONE of the {len(ANCHORS)} cardiac anchors resolved to "
          "a SpliceAI-bearing row — this is the wrong or a partial slice, "
          "rebuilding")
    sys.exit(1)
print(f"[spliceai:probe] OK — {ok_anchors}/{len(ANCHORS)} anchors resolved")
sys.exit(0)
PYEOF
        return 1
      fi
      return 0;;
    *)
      return 0;;
  esac
}

mkdir -p "$DATA_DIR"

LOCK_DIR="$DATA_DIR/.build_all.lock"
LOCK_META="$LOCK_DIR/heartbeat"
LOCK_STALE_MINUTES=${BUILD_ALL_LOCK_STALE_MINUTES:-15}
LOCK_WAIT_SECONDS=${LOCK_WAIT_SECONDS:-1200}
LOCK_POLL_SECONDS=${LOCK_POLL_SECONDS:-30}
_LOCK_HEARTBEAT_PID=""

_mtime() { stat -c %Y "$1" 2>/dev/null || stat -f %m "$1" 2>/dev/null; }

_lock_start_heartbeat() {
  date -u +%s > "$LOCK_META" 2>/dev/null || true
  ( while :; do
      sleep 60
      date -u +%s > "$LOCK_META" 2>/dev/null || exit 0
    done ) &
  _LOCK_HEARTBEAT_PID=$!
}

_lock_release() {
  [[ -n "$_LOCK_HEARTBEAT_PID" ]] && kill "$_LOCK_HEARTBEAT_PID" 2>/dev/null
  rm -rf "$LOCK_DIR" 2>/dev/null || true
}

_lock_is_stale() {
  if [[ ! -f "$LOCK_META" ]]; then
    echo "           (no heartbeat file — abandoned, or an older script version)" >&2
    return 0
  fi
  local mt now age_min
  mt="$(_mtime "$LOCK_META")"
  if [[ -z "$mt" ]]; then
    echo "           (heartbeat unreadable)" >&2
    return 0
  fi
  now="$(date -u +%s)"; age_min=$(( (now - mt) / 60 ))
  if (( age_min >= LOCK_STALE_MINUTES )); then
    echo "           (heartbeat is ${age_min}m old, threshold ${LOCK_STALE_MINUTES}m)" >&2
    return 0
  fi
  echo "           (heartbeat is ${age_min}m old — a build appears to be LIVE)" >&2
  return 1
}

_acquire_lock() {
  local waited=0 announced=""
  while :; do
    mkdir "$LOCK_DIR" 2>/dev/null && return 0

    if _lock_is_stale 2>/dev/null; then
      echo "!!! RECLAIMING a stale lock at $LOCK_DIR (its heartbeat stopped)." >&2
      echo "    If a build really is running against this volume, abort NOW —" >&2
      echo "    concurrent writers corrupt the SQLite caches." >&2
      rm -rf "$LOCK_DIR" 2>/dev/null || true
      mkdir "$LOCK_DIR" 2>/dev/null && return 0
      echo "!!! ABORT: could not reclaim $LOCK_DIR (permissions?)." >&2
      return 1
    fi

    if (( waited >= LOCK_WAIT_SECONDS )); then
      echo "!!! ABORT: $LOCK_DIR has been held for $((LOCK_WAIT_SECONDS / 60))m and its" >&2
      echo "           heartbeat keeps advancing, so a sibling replica really is" >&2
      echo "           building. Run the data-build job with" >&2
      echo "           parallelism / replica-count = 1. To clear it by hand:" >&2
      echo "             rm -rf '$LOCK_DIR'" >&2
      return 1
    fi

    if [[ -z "$announced" ]]; then
      announced=1
      echo "!!! Lock $LOCK_DIR is held and its heartbeat looks LIVE — waiting up to" >&2
      echo "    $((LOCK_WAIT_SECONDS / 60))m. A holder killed mid-build leaves a heartbeat only" >&2
      echo "    seconds old, so waiting is what tells a live build from a dead one:" >&2
      echo "    a live one keeps refreshing it, a dead one cannot." >&2
    fi
    sleep "$LOCK_POLL_SECONDS"
    waited=$(( waited + LOCK_POLL_SECONDS ))
  done
}

if [[ -z "${BUILD_ALL_NO_LOCK:-}" ]]; then
  _acquire_lock || exit 3
  _lock_start_heartbeat
  trap '_lock_release' EXIT
  trap '_lock_release; exit 143' TERM
  trap '_lock_release; exit 130' INT
fi

declare -a FAILED=()
declare -a DONE=()
declare -a SKIPPED=()
declare -a PRESENT=()
declare -a REFRESHED=()

skipped() {
  if [[ -n "$ONLY_LIST" && ",$ONLY_LIST," != *",$1,"* ]]; then return 0; fi
  [[ ",$SKIP_LIST," == *",$1,"* ]]
}

_run_with() {                 # _run_with <runner> <label> <build-script> [args...]
  local runner="$1"; local label="$2"; shift 2
  if skipped "$label"; then
    echo ">>> SKIP  $label"
    SKIPPED+=("$label")
    return 0
  fi

  local out cad; out="$(output_for "$label")"; cad="$(cadence_for "$label")"

  local force_this=""
  if [[ "$label" == vep && -n "${CHECK_ONLY:-}${FORCE_VEP:-}" ]]; then
    force_this=1
    echo ""
    echo ">>> FORCED  $label  (${CHECK_ONLY:+CHECK_ONLY}${FORCE_VEP:+FORCE_VEP} set —" \
         "bypassing the cadence KEEP so the step actually runs)"
  fi

  if [[ -z "$force_this" && -n "$RESUME$MONTHLY" && -n "$out" && -s "$out" ]] \
     && artifact_valid "$label" "$out"; then
    local keep=""
    if [[ -n "$RESUME" ]]; then
      keep="--resume: skipping everything already built"
    elif [[ "$cad" == image ]]; then
      keep="ships in the image — regenerate + commit in a checkout instead"
    elif [[ "$cad" == versioned ]]; then
      keep="versioned upstream — rebuild on a numbered release, not monthly"
    elif [[ "$cad" == static ]]; then
      keep="static upstream — only ever built when missing"
    fi
    if [[ -n "$keep" ]]; then
      echo ">>> KEEP  $label  ($keep)"
      PRESENT+=("$label")
      return 0
    fi
    echo ""
    echo ">>> REFRESH  $label  (monthly cadence — rebuilding although present)"
    REFRESHED+=("$label")
  fi
  echo ""
  echo ">>> BUILD $label  ($*)"
  local start=$SECONDS
  if "$runner" "$@"; then
    echo "<<< OK    $label  ($((SECONDS - start))s)"
    DONE+=("$label")
  else
    echo "!!! FAIL  $label  (exit $?, $((SECONDS - start))s)" >&2
    FAILED+=("$label")
  fi
}

run() {                       # python builders (the great majority)
  _run_with "$PY" "$@"
}

run_sh() {                    # shell builders (offline VEP)
  _run_with bash "$@"
}

echo "HeartVar data build → $DATA_DIR"
echo "    flags: monthly=${MONTHLY:-0} resume=${RESUME:-0} with_vep=${WITH_VEP:-0}" \
     "only=${ONLY_LIST:-<all>} skip=${SKIP_LIST:-<none>}"
echo "python: $($PY --version 2>&1)   started: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"

run uniprot          "$SCRIPT_DIR/build_uniprot_db.py" --force-download

run clinvar          "$SCRIPT_DIR/build_clinvar_db.py"
run gnomad_constraint "$SCRIPT_DIR/build_gnomad_constraint_db.py"
run gtex             "$SCRIPT_DIR/build_gtex_db.py"
run medgen           "$SCRIPT_DIR/build_medgen_db.py"
run mgi              "$SCRIPT_DIR/build_mgi_db.py"
run biogrid          "$SCRIPT_DIR/build_biogrid_db.py"
run hpo_labels       "$SCRIPT_DIR/build_hpo_labels.py"
run panelapp         "$SCRIPT_DIR/build_panelapp_snapshot.py"

run gnomad_freq      "$SCRIPT_DIR/build_gnomad_freq_db.py"
run spliceai         "$SCRIPT_DIR/build_spliceai_db.py"
run fetal_heart      "$SCRIPT_DIR/build_fetal_heart_db.py"

run alphafold        "$SCRIPT_DIR/build_alphafold_structures.py"

run hgnc_alias       "$SCRIPT_DIR/build_hgnc_alias_db.py"
run opentargets      "$SCRIPT_DIR/build_opentargets_db.py"
run gene_id_map      "$SCRIPT_DIR/build_gene_id_map.py"
run erepo            "$SCRIPT_DIR/build_erepo_dump.py"
run gencc            "$SCRIPT_DIR/build_gencc_db.py"
run clingen_gv       "$SCRIPT_DIR/build_clingen_gene_validity.py" --force-download
for f in hgnc_alias_map.json erepo_all.tsv gencc_submissions.json clingen_gene_validity.json; do
  if [[ -f "$BACKEND_DATA/$f" ]]; then
    cp -f "$BACKEND_DATA/$f" "$DATA_DIR/$f" && echo "    copied $f → $DATA_DIR"
  fi
done

if [[ ! -f "$DATA_DIR/AlphaMissense_hg38.tsv.gz" ]]; then
  echo ""
  echo "!!! NOTE  AlphaMissense_hg38.tsv.gz is MISSING from $DATA_DIR (not auto-built)."
  echo "          Exact file: 642,961,469 bytes (= 613 MiB = 643 MB),"
  echo "          MD5 9fd167735f16a1b87da6eb3e4c25fcb5. Zenodo record 8208688:"
  echo "            curl -L -o \"$DATA_DIR/AlphaMissense_hg38.tsv.gz\" \\"
  echo "              'https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz?download=1'"
  echo "          The .tbi is NOT on Zenodo — generate it in place (the file is"
  echo "          already bgzip; do NOT gunzip/re-gzip it):"
  echo "            tabix -s 1 -b 2 -e 2 \"$DATA_DIR/AlphaMissense_hg38.tsv.gz\""
  echo "          Then point ALPHAMISSENSE_PATH at the .tsv.gz."
fi

if [[ ! -f "$DATA_DIR/hg38.phyloP100way.bw" ]]; then
  echo ""
  echo "!!! NOTE  hg38.phyloP100way.bw is MISSING from $DATA_DIR (not auto-built)."
  echo "          ~9.2 GB. Must be the FULL track — a partial extract changes BP7"
  echo "          behaviour. From UCSC:"
  echo "            curl -L -o \"$DATA_DIR/hg38.phyloP100way.bw\" \\"
  echo "              'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw'"
  echo "          Then set HEARTVAR_PHYLOP_PATH. Until then BP7 conservation is"
  echo "          read live from UCSC on every synonymous variant."
fi

if ! command -v INSTALL.pl >/dev/null 2>&1; then
  echo ""
  echo "!!! NOTE  INSTALL.pl is not on PATH, so the offline-VEP step cannot run"
  echo "          here. Expected on a dev machine. In the deploy the builder"
  echo "          image is built from ensemblorg/ensembl-vep, so it is present."
  SKIP_LIST="${SKIP_LIST:+$SKIP_LIST,}vep"
fi
run_sh cdot          "$SCRIPT_DIR/build_cdot_transcripts.sh"

if [[ -z "$WITH_VEP" ]] && ! skipped vep; then
  if [[ ! -s "$(output_for vep)" ]]; then
    echo ""
    echo "!!! NOTE  offline-VEP step SKIPPED (opt-in only — pass --with-vep)."
    echo "          No manifest at $(output_for vep), so offline VEP is NOT"
    echo "          usable yet and the app stays on live Ensembl REST. Run once"
    echo "          with --with-vep to install or finish it; per-component state"
    echo "          in vep/.heartvar_vep_state/ means a completed cache is not"
    echo "          re-downloaded."
  fi
  SKIP_LIST="${SKIP_LIST:+$SKIP_LIST,}vep"
fi
run_sh vep           "$SCRIPT_DIR/setup_offline_vep.sh"

echo ""
if [[ -n "$RESUME" ]]; then   MODE="--resume (crash retry: skip everything present)"
elif [[ -n "$MONTHLY" ]]; then MODE="--monthly (refresh monthly-cadence sources)"
else                           MODE="full refresh (rebuild everything)"
fi

"$PY" - "$DATA_DIR/build_stamp.json" "$MODE" \
  "${REFRESHED[*]:-}" "${DONE[*]:-}" "${FAILED[*]:-}" "${PRESENT[*]:-}" <<'PYEOF'
import json, sys, datetime, pathlib

path = pathlib.Path(sys.argv[1])
mode, refreshed, done, failed, present = (s.split() for s in sys.argv[2:7])
mode = " ".join(mode)
now = datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()

prior = {}
try:
    prior = json.loads(path.read_text())
except Exception:
    pass

stamp = {
    "last_run_utc": now,
    # Something was built ⇒ the mirror genuinely moved. Otherwise keep whatever
    # the last real refresh was (None if we have never recorded one).
    "last_refresh_utc": now if done else prior.get("last_refresh_utc"),
    "mode": mode,
    "refreshed": refreshed,
    "built": done,
    "kept": present,
    "failed": failed,
}
try:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(stamp, indent=2) + "\n")
    tmp.replace(path)              # atomic: the app may read this concurrently
    print(f"    wrote {path.name} (last_refresh_utc={stamp['last_refresh_utc']})")
except OSError as e:
    print(f"    WARN could not write {path}: {e}", file=sys.stderr)

PYEOF

echo ""
echo "================ build_all.sh summary ================"
echo "mode: $MODE"
echo "OK        (${#DONE[@]}): ${DONE[*]:-none}"
echo "REFRESHED (${#REFRESHED[@]}): ${REFRESHED[*]:-none}   # rebuilt although already present"
echo "KEPT      (${#PRESENT[@]}): ${PRESENT[*]:-none}   # left in place (static/versioned/image, or --resume)"
echo "SKIPPED   (${#SKIPPED[@]}): ${SKIPPED[*]:-none}   # --skip"
echo "FAILED    (${#FAILED[@]}): ${FAILED[*]:-none}"
if [[ -n "$MONTHLY" && ${#REFRESHED[@]} -eq 0 && ${#DONE[@]} -eq 0 ]]; then
  echo "!!! WARN  --monthly rebuilt NOTHING. Expected the monthly-cadence sources" >&2
  echo "          (clinvar, uniprot, medgen, mgi, panelapp, hpo_labels, hgnc_alias," >&2
  echo "          erepo, gencc, clingen_gv, opentargets, biogrid) to rebuild." >&2
fi
echo "finished: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "data dir size: $(du -sh "$DATA_DIR" 2>/dev/null | cut -f1)"
echo "data dir free:"
_df_out="$(df -h "$DATA_DIR" 2>/dev/null || true)"
if [[ -n "$_df_out" ]]; then
  printf '%s\n' "$_df_out" | sed 's/^/    /'
else
  echo "    (df returned nothing for $DATA_DIR — free space UNKNOWN on this"
  echo "     mount. The 28 GB floor for the VEP cache cannot be checked, so an"
  echo "     install may fail part-way instead of refusing up front.)"
fi

ls -t "${_HV_LOG_DIR}"/*.log 2>/dev/null | tail -n +6 | xargs -r rm -f || true

[[ ${#FAILED[@]} -eq 0 ]]
