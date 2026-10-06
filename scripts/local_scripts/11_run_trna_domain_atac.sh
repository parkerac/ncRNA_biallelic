#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

BW_URL="https://s3.us-east-1.amazonaws.com/ox.paediatrics.sanders.genome-browser-public/Narjes_tracks/SHSY5Y.k1m1l75best_treat_afterfiting_all.RPM_hg38_chrMain.bw"
BED="$REPO_DIR/hg38-tRNA-domains_atac.bed"
OUT="$REPO_DIR/SHSY5Y_ATAC.tsv"
TMP_BED="$(mktemp "${TMPDIR:-/tmp}/trna_domains_atac.XXXXXX.bed")"
TMP_TRNA_BED="$(mktemp "${TMPDIR:-/tmp}/trna_whole_atac.XXXXXX.bed")"
TMP_DOMAIN_OUT="$(mktemp "${TMPDIR:-/tmp}/trna_domains_atac.XXXXXX.tsv")"
TMP_TRNA_OUT="$(mktemp "${TMPDIR:-/tmp}/trna_whole_atac.XXXXXX.tsv")"
trap 'rm -f "$TMP_BED" "$TMP_TRNA_BED" "$TMP_DOMAIN_OUT" "$TMP_TRNA_OUT"' EXIT

tail -n +2 "$BED" > "$TMP_BED"

awk 'BEGIN { FS = OFS = "\t" }
  NR == 1 { next }
  {
    sub(/\r$/, "", $4)
    label = $4
    sub(/\|.*/, "", label)
    key = $1 SUBSEP label
    if (!(key in seen)) {
      seen[key] = 1
      chrom[key] = $1
      name[key] = label
      start[key] = $2
      end[key] = $3
      order[++n] = key
    }
    if ($2 < start[key]) start[key] = $2
    if ($3 > end[key]) end[key] = $3
  }
  END {
    for (i = 1; i <= n; i++) {
      key = order[i]
      print chrom[key], start[key], end[key], name[key]
    }
  }' "$BED" > "$TMP_TRNA_BED"

bigWigAverageOverBed "$BW_URL" "$TMP_BED" "$TMP_DOMAIN_OUT"
bigWigAverageOverBed "$BW_URL" "$TMP_TRNA_BED" "$TMP_TRNA_OUT"
cat "$TMP_DOMAIN_OUT" "$TMP_TRNA_OUT" > "$OUT"
