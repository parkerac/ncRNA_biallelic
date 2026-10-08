#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

BW_URL="https://s3.us-east-1.amazonaws.com/ox.paediatrics.sanders.genome-browser-public/Narjes_tracks/SHSY5Y.k1m1l75best_treat_afterfiting_all.RPM_hg38_chrMain.bw"
BED="$REPO_DIR/hg38-tRNA-domains_atac.bed"
PRIMATE_PHYLOP_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP17way/hg38.phyloP17way.bw"
VERTEBRATE_PHYLOP_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/trna_ucsc_data.XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT

DOMAIN_BED="$TMP_DIR/trna_domains.bed"
TRNA_BED="$TMP_DIR/trna_whole.bed"
tail -n +2 "$BED" > "$DOMAIN_BED"

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
  }' "$BED" > "$TRNA_BED"

summarize_bigwig() {
  local url="$1"
  local out="$2"
  local chrom
  local chrom_bed="$TMP_DIR/chrom.bed"
  local chrom_out="$TMP_DIR/chrom.tsv"

  echo "Writing $(basename "$out")" >&2
  : > "$out"
  for chrom in $(cut -f1 "$DOMAIN_BED" | uniq); do
    awk -v chrom="$chrom" 'BEGIN { FS = OFS = "\t" } $1 == chrom' "$DOMAIN_BED" > "$chrom_bed"
    echo "  domains $chrom" >&2
    bigWigAverageOverBed "$url" "$chrom_bed" "$chrom_out"
    cat "$chrom_out" >> "$out"
  done
  for chrom in $(cut -f1 "$TRNA_BED" | uniq); do
    awk -v chrom="$chrom" 'BEGIN { FS = OFS = "\t" } $1 == chrom' "$TRNA_BED" > "$chrom_bed"
    echo "  whole-tRNA $chrom" >&2
    bigWigAverageOverBed "$url" "$chrom_bed" "$chrom_out"
    cat "$chrom_out" >> "$out"
  done
}

summarize_bigwig "$BW_URL" "$REPO_DIR/SHSY5Y_ATAC.tsv"
summarize_bigwig "$PRIMATE_PHYLOP_URL" "$REPO_DIR/hg38_phyloP17way_tRNA_domains.tsv"
summarize_bigwig "$VERTEBRATE_PHYLOP_URL" "$REPO_DIR/hg38_phyloP100way_tRNA_domains.tsv"
