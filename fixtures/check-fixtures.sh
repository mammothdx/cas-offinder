#!/usr/bin/env bash
set -euo pipefail

binary="${1:-./cas-offinder}"
device="${2:-C}"
root="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${root}/.." && pwd)"
tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/cas-offinder.XXXXXX")"
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
run_case "$root/pam-both-nonly"

echo "All fixture expectations were observed for device ${device}."

# Phase 2: brute-force validation (exact set comparison)
if command -v python3 >/dev/null 2>&1; then
  bf_script="$root/bruteforce_check.py"
  if [ -f "$bf_script" ]; then
    for case_dir in "$root/3prime" "$root/5prime" "$root/pam-both" "$root/pam-both-nonly"; do
      out_file="$tmpdir/$(basename "$case_dir").out"
      python3 "$bf_script" "$case_dir" "$out_file" || exit 1
    done
    echo "Brute-force validation passed for device ${device}."
  fi
fi
