#!/usr/bin/env bash
set -uo pipefail

DATA_DIR="${DATA_DIR:-${1:-/app/data}}"
OUT_DIR="$DATA_DIR/cdot"
CDOT_DATA_VERSION="${CDOT_DATA_VERSION:-0.2.34}"
VEP_REL="${VEP_RELEASE_NUM:-}"
if [[ -z "$VEP_REL" ]] && command -v vep >/dev/null 2>&1; then
  _c="$(dirname "$(command -v vep)")/modules/Bio/EnsEMBL/VEP/Constants.pm"
  [[ -r "$_c" ]] && VEP_REL="$(sed -n 's/.*VEP_VERSION[[:space:]]*=[[:space:]]*\([0-9][0-9]*\).*/\1/p' "$_c" | head -1)"
fi
VEP_REL="${VEP_REL:-113}"

BASE="https://github.com/SACGF/cdot/releases/download/data_v${CDOT_DATA_VERSION}"
REFSEQ_FILE="cdot-${CDOT_DATA_VERSION}.refseq.GRCh38.json.gz"
ENSEMBL_FILE="cdot-${CDOT_DATA_VERSION}.Homo_sapiens_GRCh38_Ensembl_${VEP_REL}.gtf.json.gz"

echo "==> cdot transcript tables"
echo "    cdot data version : $CDOT_DATA_VERSION"
echo "    Ensembl release   : $VEP_REL (matched to the VEP cache)"
echo "    destination       : $OUT_DIR"

mkdir -p "$OUT_DIR" || { echo "!!! cannot create $OUT_DIR" >&2; exit 1; }

avail_gb="$(df -P -k "$OUT_DIR" 2>/dev/null | awk 'NR==2 {print int($4/1024/1024)}')"
if [ "${avail_gb:-0}" -lt 1 ]; then
  echo "!!! only ${avail_gb:-0} GB free at $OUT_DIR; need >= 1 GB." >&2
  exit 1
fi

FAILED=0
fetch() {                     # fetch <filename>
  local name="$1" url="$BASE/$1" dest="$OUT_DIR/$1" tmp="$OUT_DIR/.$1.part"
  if [[ -s "$dest" ]] && gzip -t "$dest" 2>/dev/null; then
    echo "    KEEP $name (present and readable)"
    return 0
  fi
  echo "    GET  $name"
  if ! curl -fSL --retry 3 --retry-delay 10 -o "$tmp" "$url"; then
    echo "!!! download failed: $url" >&2
    rm -f "$tmp"
    return 1
  fi
  if ! gzip -t "$tmp" 2>/dev/null; then
    echo "!!! $name is not a valid gzip (truncated download)" >&2
    rm -f "$tmp"
    return 1
  fi
  mv -f "$tmp" "$dest"
  echo "         $(du -h "$dest" | cut -f1)"
  return 0
}

fetch "$REFSEQ_FILE"  || FAILED=$((FAILED + 1))
fetch "$ENSEMBL_FILE" || FAILED=$((FAILED + 1))

if [[ "$FAILED" -ne 0 ]]; then
  echo "!!! $FAILED cdot file(s) unavailable — HGVS input will keep using REST." >&2
  echo "    Egress allowlist must include github.com and objects.githubusercontent.com." >&2
  exit 1
fi

PY_BIN="${PY_BIN:-$( [[ -x /opt/venv/bin/python ]] && echo /opt/venv/bin/python || echo python3 )}"
echo "==> converting to SQLite with $PY_BIN"
if ! VEP_RELEASE_NUM="$VEP_REL" CDOT_DATA_VERSION="$CDOT_DATA_VERSION" \
     "$PY_BIN" "$(dirname "${BASH_SOURCE[0]}")/build_cdot_db.py" \
       --out "$OUT_DIR" --json-dir "$OUT_DIR"; then
  echo "!!! could not build the cdot SQLite database. The JSONs are staged but" >&2
  echo "    unusable in production (4.5 GB RSS), so this is a failure." >&2
  exit 1
fi

if command -v python3 >/dev/null 2>&1 || [[ -x /opt/venv/bin/python ]]; then
  PY_BIN="${PY_BIN:-$( [[ -x /opt/venv/bin/python ]] && echo /opt/venv/bin/python || echo python3 )}"
  echo "    verifying with $PY_BIN"
  if HEARTVAR_HGVS_RESOLVER=1 \
     HEARTVAR_CDOT_DB="$OUT_DIR/cdot_transcripts.db" \
     "$PY_BIN" - <<'VERIFY'
import os, sys
sys.path.insert(0, os.environ.get("HEARTVAR_REPO_ROOT", "/app"))
try:
    from backend.clients import hgvs_resolver as hr
except Exception as exc:
    print(f"    NOTE cannot import the resolver here ({exc!r}) — files staged, "
          f"unverified", file=sys.stderr)
    raise SystemExit(0)
r = hr.resolve("c.1208G>A", "NM_000257.4", "MYH7")
if r is None:
    print("    RESOLVER RETURNED None for MYH7 R403Q", file=sys.stderr)
    raise SystemExit(1)
got = (r.chrom, r.vep_pos, r.vep_ref, r.vep_alt)
if got != ("14", 23429278, "C", "T"):
    print(f"    WRONG COORDINATES for MYH7 R403Q: {got}", file=sys.stderr)
    raise SystemExit(1)
print(f"    OK  MYH7 R403Q -> {r.chrom}:{r.vep_pos} {r.vep_ref}>{r.vep_alt} "
      f"| gnomAD key {r.gnomad_id()}")

# ⚠⚠ THE RESOURCE CHECK, AND THIS IS THE ONE THAT WAS MISSING.
# The resolver merged in #41 was correct and UNDEPLOYABLE: cdot's
# JSONDataProvider loads all 906,754 transcripts into RAM — 4480 MB against a
# 4 GiB container. Every test passed. They passed because the unit tests stub the
# mapper and the real-data test runs in a process with the whole machine to
# itself, so nothing in the suite could see a resource cost at all.
#
# This check runs HERE because here is the only place that has the real data, in
# the real container, at the real memory limit. Same reasoning as
# setup_offline_vep.sh's _check_runtime_call: a verified component is not a
# verified deployment, and what fails is the thing running where it will run.
from backend.clients.cdot_sqlite import RSS_CEILING_MB, current_rss_mb

rss = current_rss_mb()
print(f"    RSS after provider + resolve: {rss:.0f} MB (ceiling {RSS_CEILING_MB} MB)")
if rss > RSS_CEILING_MB:
    print(f"    RESIDENT MEMORY {rss:.0f} MB EXCEEDS THE {RSS_CEILING_MB} MB "
          f"CEILING. The container has 4 GiB and the app needs the rest of it. "
          f"If this jumped to ~4.5 GB the resolver fell back to the in-memory "
          f"JSON provider — check HEARTVAR_CDOT_DB is set and the .db exists.",
          file=sys.stderr)
    raise SystemExit(1)
VERIFY
  then
    :
  else
    echo "!!! the resolver could not use the staged files. Not reporting success," >&2
    echo "    so the next build retries." >&2
    exit 1
  fi
fi

echo "==> cdot transcript database ready."
echo "    HEARTVAR_CDOT_DB=$OUT_DIR/cdot_transcripts.db"
echo "    (HEARTVAR_HGVS_RESOLVER stays UNSET until the parity gate is signed off)"
exit 0
