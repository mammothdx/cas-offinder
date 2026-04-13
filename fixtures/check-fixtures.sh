#!/usr/bin/env bash
set -euo pipefail

binary="${1:-./cas-offinder}"
device="${2:-C}"
root="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${root}/.." && pwd)"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

run_case() {
  local case_dir="$1"
  local out_file="$tmpdir/$(basename "$case_dir").out"
  local actual_file="$tmpdir/$(basename "$case_dir").actual"

  (
    cd "$repo_root"
    "$binary" "$case_dir/input.txt" "$device" "$out_file" >/dev/null
  )
  awk -F'\t' 'BEGIN{OFS=" "} !/^#/ {print $1, $2, $7, $9, $6, $3, $4}' "$out_file" | sort -u >"$actual_file"

  while read -r id bulge_type direction bulge_size index seq_rna seq_dna; do
    [[ -z "${id}" || "${id}" == \#* ]] && continue
    if ! grep -Fqx "${id} ${bulge_type} ${direction} ${bulge_size} ${index} ${seq_rna} ${seq_dna}" "$actual_file"; then
      echo "Fixture failure in ${case_dir}: missing ${id} ${bulge_type} ${direction} ${bulge_size} ${index} ${seq_rna} ${seq_dna}" >&2
      echo "Actual rows:" >&2
      cat "$actual_file" >&2
      return 1
    fi
  done <"$case_dir/expected.txt"
}

run_case "$root/3prime"
run_case "$root/5prime"
run_case "$root/pam-both"

echo "All fixture expectations were observed for device ${device}."
