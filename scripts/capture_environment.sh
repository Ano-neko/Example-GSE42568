#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PROJECT_DIR="$(cd "$PROJECT_DIR" && pwd -P)"
ENV_PREFIX="${BIOINFO_ENV_PREFIX:-${CONDA_PREFIX:-}}"

if [[ -n "$ENV_PREFIX" ]]; then
  R_BIN="${R_BIN:-$ENV_PREFIX/bin/Rscript}"
  PY_BIN="${PY_BIN:-$ENV_PREFIX/bin/python}"
else
  R_BIN="${R_BIN:-Rscript}"
  PY_BIN="${PY_BIN:-python}"
fi

if [[ -n "${CONDA_EXE:-}" && -x "$CONDA_EXE" ]]; then
  CONDA_BIN="$CONDA_EXE"
else
  CONDA_BIN="${CONDA_BIN:-conda}"
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
require_runtime "$CONDA_BIN"
cd "$PROJECT_DIR"

session_tmp="$(mktemp)"
trap 'find "$session_tmp" -maxdepth 0 -type f -delete' EXIT
"$R_BIN" -e 'sessionInfo()' > "$session_tmp"
if [[ -n "$ENV_PREFIX" ]]; then
  "$PY_BIN" -c 'from pathlib import Path; import sys; Path(sys.argv[3]).write_text(Path(sys.argv[1]).read_text().replace(sys.argv[2], "${BIOINFO_ENV_PREFIX}"))' \
    "$session_tmp" "$ENV_PREFIX" reports/sessionInfo.txt
  "$CONDA_BIN" env export --prefix "$ENV_PREFIX" --no-builds | sed '/^prefix: /d' > config/environment.yml
  "$CONDA_BIN" list --prefix "$ENV_PREFIX" --explicit > config/conda-explicit-linux-64.txt
else
  cp "$session_tmp" reports/sessionInfo.txt
  "$CONDA_BIN" env export --no-builds | sed '/^prefix: /d' > config/environment.yml
  "$CONDA_BIN" list --explicit > config/conda-explicit-linux-64.txt
fi
"$PY_BIN" -m pip freeze > config/python-freeze.txt

printf 'Environment captured: R session, pinned Conda YAML/explicit export, and pip freeze.\n'
