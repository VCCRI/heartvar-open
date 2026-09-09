#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "[entrypoint_builder] Container started at $(date -u '+%Y-%m-%dT%H:%M:%SZ')"

BUILD_ALL_ARGS="${BUILD_ALL_ARGS:---monthly}"

_IS_JOB=""
if [[ -n "${EXIT_WHEN_DONE:-}${CONTAINER_APP_JOB_NAME:-}${CONTAINER_APP_JOB_EXECUTION_NAME:-}" ]]; then
  _IS_JOB=1
fi

if [[ -z "$_IS_JOB" && -n "${ONLY:-}${SKIP:-}" ]]; then
  echo "[entrypoint_builder] ⚠ IGNORING ONLY='${ONLY:-}' SKIP='${SKIP:-}' — this is a" \
       "webapp container start (the scheduled monthly refresh), not a job dispatch," \
       "so it must cover the whole monthly set. Unset ONLY/SKIP in the app settings" \
       "to silence this; dispatch the data-build job if a narrowed run is intended."
  unset ONLY SKIP
fi

STAMP="$(cd "$SCRIPT_DIR/.." && pwd)/data/build_stamp.json"
REFRESH_DUE=1
if [[ -s "$STAMP" ]]; then
  LAST_REFRESH="$(sed -n 's/.*"last_refresh_utc"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' \
                  "$STAMP" | head -1)"
  if [[ -n "$LAST_REFRESH" && "${LAST_REFRESH:0:7}" == "$(date -u '+%Y-%m')" ]]; then
    REFRESH_DUE=0
  fi
fi
if [[ -n "${FORCE_BUILD:-}" ]]; then
  echo "[entrypoint_builder] FORCE_BUILD set — running regardless of the stamp."
  REFRESH_DUE=1
fi

_EXPLICIT_REQUEST=""
if [[ -n "${CHECK_ONLY:-}${ONLY:-}" ]]; then
  _EXPLICIT_REQUEST=1
fi
if [[ -n "$_EXPLICIT_REQUEST" && "$REFRESH_DUE" -eq 0 ]]; then
  echo "[entrypoint_builder] CHECK_ONLY/ONLY set — this is an explicit request," \
       "so the once-per-calendar-month guard does not apply."
  REFRESH_DUE=1
fi

BUILD_EXIT=0
BUILD_RAN=""
if [[ "$REFRESH_DUE" -eq 0 ]]; then
  echo "[entrypoint_builder] SKIPPING the data refresh: build_stamp.json already" \
       "records one this month (last_refresh_utc=$LAST_REFRESH). Set FORCE_BUILD=1" \
       "for a deliberate off-cycle run."
else
  BUILD_RAN=1
  echo "[entrypoint_builder] Running build_all.sh $BUILD_ALL_ARGS ..."
  # shellcheck disable=SC2086 — intentional word-splitting so BUILD_ALL_ARGS can
  bash "$SCRIPT_DIR/build_all.sh" $BUILD_ALL_ARGS
  BUILD_EXIT=$?
fi

if [[ -z "$BUILD_RAN" ]]; then
  echo "[entrypoint_builder] build_all.sh was NOT RUN (see the SKIP above)." \
       "Nothing was built, refreshed or checked."
elif [[ $BUILD_EXIT -eq 0 ]]; then
  echo "[entrypoint_builder] build_all.sh completed successfully."
else
  echo "[entrypoint_builder] build_all.sh finished with exit $BUILD_EXIT (some caches failed)." >&2
fi

_VEP_DIR="$(cd "$SCRIPT_DIR/.." && pwd)/data/vep"
if [[ -z "$_IS_JOB" ]]; then
  if [[ ! -s "$_VEP_DIR/.heartvar_vep_manifest.json" ]]; then
    echo "[entrypoint_builder] offline VEP is not complete, and this is the webapp" \
         "sidecar — NOT attempting it. Run the data-build workflow to install or" \
         "finish it."
  fi
elif [[ -n "${CHECK_ONLY:-}" ]]; then
  echo "[entrypoint_builder] CHECK_ONLY set — not attempting the offline-VEP" \
       "completion pass; the probe above is the whole run."
elif [[ ! -s "$_VEP_DIR/.heartvar_vep_manifest.json" ]]; then
  if [[ -s "$_VEP_DIR/.heartvar_vep_state/cache.json" ]]; then
    echo "[entrypoint_builder] offline-VEP install is PARTIAL (cache verified, no" \
         "manifest) — finishing it with --only vep --with-vep."
    bash "$SCRIPT_DIR/build_all.sh" --only vep --with-vep
    VEP_EXIT=$?
    echo "[entrypoint_builder] vep completion pass exited $VEP_EXIT."
    if [[ "$BUILD_EXIT" -eq 0 ]]; then BUILD_EXIT=$VEP_EXIT; fi
  else
    echo "[entrypoint_builder] offline VEP is NOT installed and there is no partial" \
         "state. Leaving it alone — a first install is a ~23 GB download and is" \
         "opt-in (run build_all.sh --only vep --with-vep deliberately)."
  fi
fi

if [[ -n "${STAY_ALIVE:-}" ]]; then
  echo "[entrypoint_builder] STAY_ALIVE set — entering idle sleep."
  exec sleep infinity
fi
if [[ -n "$_IS_JOB" ]]; then
  echo "[entrypoint_builder] Running as a Container App Job — exiting $BUILD_EXIT" \
       "so the execution reports a real status instead of Running forever."
  exit "$BUILD_EXIT"
fi

echo "[entrypoint_builder] Entering idle sleep — container staying alive."
exec sleep infinity
