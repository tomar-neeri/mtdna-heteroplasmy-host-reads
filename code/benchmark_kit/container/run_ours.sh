#!/usr/bin/env bash
# Runs INSIDE the container. Applies the in-house caller to the mixture BAMs with
# the cohort's fitted error model, using the same thresholds as the Snakefile.
#
#   run_ours.sh <bam_dir> <out_dir> <summary_dir> <pon|nopon>
#
#   pon    as published: substitution model + leave-one-out panel of normals
#   nopon  substitution model only
#
# Why both: the panel of normals sums alt reads from EVERY cohort sample,
# including samples homoplasmic for the allele. Mixture truth sites are, by
# construction, positions where cohort samples differ homoplasmically, so the
# published caller is expected to lose sensitivity there. Running both arms
# measures that trade-off instead of hiding it.
set -euo pipefail
BAMDIR=$1; OUT=$2; SUM=$3; MODE=${4:-pon}
mkdir -p "$OUT"
[ -s "$SUM/error_model.json" ] || { echo "missing $SUM/error_model.json" >&2; exit 1; }
PON=()
if [ "$MODE" = pon ]; then
  [ -s "$SUM/site_background.tsv" ] || { echo "missing $SUM/site_background.tsv" >&2; exit 1; }
  PON=(--site-background "$SUM/site_background.tsv" --min-loo-samples 3)
fi
shopt -s nullglob
n=0
for bam in "$BAMDIR"/*.bam; do
  s=$(basename "$bam" .bam)
  [ -s "$OUT/$s.calls.tsv" ] && continue
  python3 /opt/scripts/call_heteroplasmy.py \
      --bam "$bam" --ref "${RCRS:-/opt/refs/rCRS.fasta}" \
      --model "$SUM/error_model.json" "${PON[@]}" \
      --out "$OUT/$s.calls.tsv" --summary "$OUT/$s.summary.json" \
      --contig chrM --min-depth 50 --min-baseq 20 --min-mapq 30 \
      --min-alt-reads 4 --vaf-floor 0.01 --vaf-homoplasmy 0.95 --sig-level 0.05 \
      2> "$OUT/$s.log"
  n=$((n+1))
done
echo "  in-house caller ($MODE): $n new mixture(s) called -> $OUT"
