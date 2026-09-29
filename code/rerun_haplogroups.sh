#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# rerun_haplogroups.sh — apply the v1.6 haplogroup fix to an existing run
#
# Reruns only: haplogroup.py -> lineage_matrix.py -> collect_variants.py,
# then draws Figure 3 c-e. Uses the v1.5 image with the new scripts mounted
# over the baked-in ones, so no rebuild and no 6-hour rerun is needed.
#
#   bash rerun_haplogroups.sh
#
# Put haplogroup.py and plot_haplogroups.py in $SRC/scripts/ first.
# ---------------------------------------------------------------------------
set -euo pipefail

SRC="${SRC:-/mnt/d/Cancer_microbiome/Pipeline_container}"
WORK="${WORK:-$HOME/mtdna-cohort47}"
REFS="${REFS:-$HOME/ref}"
IMAGE="${IMAGE:-mtdna-pipeline:1.5.0}"
FIGDIR="${FIGDIR:-$WORK/results/figures}"

cd "$WORK"
[ -d results/chrM ] || { echo "no results/chrM in $WORK"; exit 1; }
for f in haplogroup.py plot_haplogroups.py; do
  [ -f "$SRC/scripts/$f" ] || { echo "missing $SRC/scripts/$f"; exit 1; }
done

# keep the v1.5 outputs for comparison
BK="results/summary/pre_v1.6_$(date +%Y%m%d_%H%M)"
mkdir -p "$BK"
for f in haplogroups.tsv heteroplasmy_lineage_flags.tsv lineage_matrix.tsv \
         variants_annotated.tsv variant_recurrence.tsv; do
  [ -f "results/summary/$f" ] && cp "results/summary/$f" "$BK/"
done
echo "  v1.5 outputs backed up to $BK"

run(){ docker run --rm --user "$(id -u):$(id -g)" \
         -v "$WORK":/work -w /work -v "$REFS":/refs:ro \
         -v "$SRC/scripts/haplogroup.py":/opt/scripts/haplogroup.py:ro \
         -v "$SRC/scripts/plot_haplogroups.py":/opt/scripts/plot_haplogroups.py:ro \
         "$IMAGE" "$@"; }

echo; echo "== haplogroup.py (BAM typing)"
run python3 /opt/scripts/haplogroup.py \
    --calls results/calls --bams $(ls results/chrM/*.chrM.bam) \
    --out-dir results/summary --contig chrM \
    --min-baseq 20 --min-mapq 30 --min-depth-call 50 --vaf-homoplasmy 0.95

echo; echo "== lineage_matrix.py"
run python3 /opt/scripts/lineage_matrix.py \
    --calls results/calls --haplogroups results/summary/haplogroups.tsv \
    --out results/summary/lineage_matrix.tsv --vaf-homoplasmy 0.95

echo; echo "== collect_variants.py"
run python3 /opt/scripts/collect_variants.py \
    --calls results/calls --ref /opt/refs/rCRS.fasta --outdir results/summary \
    --common-af 0.001 --lineage-flags results/summary/heteroplasmy_lineage_flags.tsv \
    --gnomad /refs/gnomad_mt.tsv

echo; echo "== Figure 3 c-e"
mkdir -p "$FIGDIR"
run python3 /opt/scripts/plot_haplogroups.py \
    --summary results/summary --out "results/figures/Figure3_cde"

echo; echo "== what changed"
if [ -f "$BK/haplogroups.tsv" ]; then
  join -t $'\t' <(tail -n +2 "$BK/haplogroups.tsv" | cut -f1,2 | sort) \
                <(tail -n +2 results/summary/haplogroups.tsv | cut -f1,2 | sort) \
    | awk -F'\t' '$2!=$3{print "  "$1": "$2" -> "$3; n++} END{print "  "n+0" sample(s) changed"}'
fi
