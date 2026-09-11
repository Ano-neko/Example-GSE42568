#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PROJECT_DIR="$(cd "$PROJECT_DIR" && pwd -P)"

if [[ -n "${BIOINFO_ENV_PREFIX:-}" ]]; then
  R_BIN="${R_BIN:-$BIOINFO_ENV_PREFIX/bin/Rscript}"
  PY_BIN="${PY_BIN:-$BIOINFO_ENV_PREFIX/bin/python}"
else
  R_BIN="${R_BIN:-Rscript}"
  PY_BIN="${PY_BIN:-python}"
fi

require_runtime() {
  local runtime="$1"
  if [[ "$runtime" == */* ]]; then
    [[ -x "$runtime" ]] || { printf 'Required runtime is not executable: %s\n' "$runtime" >&2; exit 127; }
  else
    command -v "$runtime" >/dev/null 2>&1 || { printf 'Required command is not on PATH: %s\n' "$runtime" >&2; exit 127; }
  fi
}

require_runtime "$R_BIN"
require_runtime "$PY_BIN"
export TZ=UTC
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
cd "$PROJECT_DIR"
mkdir -p results/logs

start=$(date -u +%Y-%m-%dT%H:%M:%SZ)
printf 'Pipeline start: %s\n' "$start" | tee results/logs/pipeline.log

run_logged() {
  local log_path="$1"
  shift
  "$@" 2>&1 | tee "$log_path" | tee -a results/logs/pipeline.log
}

run_logged results/logs/00_download_inputs.log bash scripts/download_inputs.sh "$PROJECT_DIR"
run_logged results/logs/01_acquire_validate.log "$R_BIN" R/01_acquire_validate.R --project "$PROJECT_DIR"
run_logged results/logs/02_qc_annotation.log "$R_BIN" R/02_qc_annotation.R --project "$PROJECT_DIR"
run_logged results/logs/02b_qc_review.log "$R_BIN" R/02b_qc_review.R --project "$PROJECT_DIR"
run_logged results/logs/03_de_sensitivity.log "$R_BIN" R/03_de_sensitivity.R --project "$PROJECT_DIR"
run_logged results/logs/04_enrichment.log "$R_BIN" R/04_enrichment.R --project "$PROJECT_DIR"
run_logged results/logs/tests.log "$R_BIN" tests/run_tests.R --project "$PROJECT_DIR"
run_logged results/logs/05_scientific_qa.log "$R_BIN" R/05_scientific_qa.R --project "$PROJECT_DIR"
run_logged results/logs/build_reports.log "$PY_BIN" scripts/build_reports.py --project "$PROJECT_DIR"
run_logged results/logs/capture_environment.log bash scripts/capture_environment.sh "$PROJECT_DIR"
run_logged results/logs/public_release_qa.log "$PY_BIN" scripts/validate_public_release.py --project "$PROJECT_DIR"

end=$(date -u +%Y-%m-%dT%H:%M:%SZ)
printf 'Pipeline end: %s\n' "$end" | tee -a results/logs/pipeline.log
"$PY_BIN" scripts/finalize_release.py --project "$PROJECT_DIR"
