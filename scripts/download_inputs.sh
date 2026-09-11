#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$PROJECT_DIR"
mkdir -p data/raw data/metadata

MATRIX_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42568/matrix/GSE42568_series_matrix.txt.gz"
MATRIX_FILE="data/raw/GSE42568_series_matrix.txt.gz"
MATRIX_SIZE=22565943
MATRIX_SHA256="43b93f72f1d81838dcd2a5b6ccc3bc1875e46230dad0d7346236773a99e32b76"
GPL_URL="https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
GPL_FILE="data/raw/GPL570.annot.gz"

download_one() {
  local url="$1" output="$2" headers="$3"
  if [[ -s "$output" ]]; then
    return 0
  fi
  curl -fL --retry 12 --retry-all-errors --retry-delay 10 \
    --connect-timeout 30 -D "$headers" "$url" -o "$output"
}

download_matrix_four_ranges() {
  local part_dir="${MATRIX_FILE}.parts.$$"
  local parts=4 chunk start end i name
  chunk=$(( (MATRIX_SIZE + parts - 1) / parts ))
  mkdir -p "$part_dir"
  cleanup_parts() {
    find "$part_dir" -type f -delete 2>/dev/null || true
    rmdir "$part_dir" 2>/dev/null || true
  }
  trap cleanup_parts EXIT
  for ((i=0; i<parts; i++)); do
    start=$((i * chunk))
    end=$((start + chunk - 1))
    (( end >= MATRIX_SIZE )) && end=$((MATRIX_SIZE - 1))
    printf -v name '%03d.part' "$i"
    curl -fsSL --retry 20 --retry-all-errors --retry-delay 10 \
      --connect-timeout 30 -r "${start}-${end}" "$MATRIX_URL" \
      -o "$part_dir/$name" &
  done
  wait
  cat "$part_dir"/*.part > "$MATRIX_FILE"
  cleanup_parts
  trap - EXIT
}

if [[ ! -s "$MATRIX_FILE" ]]; then
  download_matrix_four_ranges
fi

actual_size=$(stat -c %s "$MATRIX_FILE")
[[ "$actual_size" -eq "$MATRIX_SIZE" ]] || {
  printf 'Matrix size mismatch: expected %s, observed %s\n' "$MATRIX_SIZE" "$actual_size" >&2
  exit 2
}
gzip -t "$MATRIX_FILE"
observed_sha=$(sha256sum "$MATRIX_FILE" | awk '{print $1}')
[[ "$observed_sha" == "$MATRIX_SHA256" ]] || {
  printf 'Matrix SHA-256 mismatch: expected %s, observed %s\n' "$MATRIX_SHA256" "$observed_sha" >&2
  exit 3
}

download_one "$GPL_URL" "$GPL_FILE" data/metadata/GPL570_annotation.headers.txt
gzip -t "$GPL_FILE"

curl -fsSI "$MATRIX_URL" > data/metadata/GSE42568_series_matrix.headers.txt || true
sha256sum "$MATRIX_FILE" "$GPL_FILE" > data/metadata/download_checksums.sha256
