#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

BW_URL="https://s3.us-east-1.amazonaws.com/ox.paediatrics.sanders.genome-browser-public/Narjes_tracks/SHSY5Y.k1m1l75best_treat_afterfiting_all.RPM_hg38_chrMain.bw"
BED="$REPO_DIR/hg38-tRNA-domains_atac.bed"
OUT="$REPO_DIR/SHSY5Y_ATAC.tsv"
TMP_BED="$(mktemp "${TMPDIR:-/tmp}/trna_domains_atac.XXXXXX.bed")"
trap 'rm -f "$TMP_BED"' EXIT

tail -n +2 "$BED" > "$TMP_BED"

bigWigAverageOverBed "$BW_URL" "$TMP_BED" "$OUT"
