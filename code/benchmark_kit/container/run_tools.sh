#!/usr/bin/env bash
# Runs INSIDE the container. Calls variants with Mutect2 (mitochondria mode) and
# mutserve on every BAM in a directory. Resumable: finished samples are skipped.
#
#   run_tools.sh <bam_dir> <out_dir> <hg38|chrM> [threads]
#
#   hg38  BAMs whose header lists all hg38 contigs (the pipeline's final chrM BAMs)
#   chrM  BAMs aligned to the mitochondrial reference alone (Pass-1 recruitment)
#
# Settings follow the GATK mitochondria pipeline (Laricchia et al. 2022) and the
# mutserve defaults, except that minimum mapping quality is 30 in both, to match
# the in-house caller. Variant filtering by allele fraction is NOT applied here;
# compare.py applies the same thresholds to every caller.
set -euo pipefail
BAMDIR=$1; OUT=$2; REFMODE=$3; T=${4:-4}
R=/bench/ref
case "$REFMODE" in
  hg38) M2REF=$R/hg38.fa ;;
  chrM) M2REF=$R/chrM.fa ;;
  *) echo "ref mode must be hg38 or chrM" >&2; exit 2 ;;
esac
MS_JAR=${MUTSERVE_JAR:-/opt/tools/mutserve/mutserve.jar}
mkdir -p "$OUT"/{mutect2,mutserve,logs,tmp}

shopt -s nullglob
bams=("$BAMDIR"/*.bam)
[ ${#bams[@]} -gt 0 ] || { echo "no BAMs in $BAMDIR" >&2; exit 1; }
echo "  ${#bams[@]} BAMs in $BAMDIR  ->  $OUT  (Mutect2 reference: $REFMODE)"

i=0
for bam in "${bams[@]}"; do
  i=$((i+1))
  s=$(basename "$bam"); s=${s%.bam}; s=${s%.chrM}; s=${s%.cand}; s=${s%.mix}
  m2=$OUT/mutect2/$s.filtered.vcf.gz
  ms=$OUT/mutserve/$s.txt
  printf "  [%2d/%d] %-28s" "$i" "${#bams[@]}" "$s"

  # --- Mutect2, mitochondria mode -------------------------------------------
  if [ ! -s "$m2" ]; then
    # the pipeline's BAMs carry no read group; Mutect2 requires one
    rg=$OUT/tmp/$s.rg.bam
    samtools addreplacerg -r "@RG\tID:$s\tSM:$s\tPL:ILLUMINA" -o "$rg" "$bam" 2> "$OUT/logs/$s.rg.log"
    samtools index "$rg"
    gatk --java-options "-Xmx4g" Mutect2 \
        -R "$M2REF" -I "$rg" -L chrM \
        --mitochondria-mode \
        --max-mnp-distance 0 \
        --max-reads-per-alignment-start 75 \
        --minimum-mapping-quality 30 \
        --annotation StrandBiasBySample \
        -O "$OUT/tmp/$s.unfiltered.vcf.gz" > "$OUT/logs/$s.mutect2.log" 2>&1 \
      || true
    [ -s "$OUT/tmp/$s.unfiltered.vcf.gz" ] && gatk --java-options "-Xmx4g" FilterMutectCalls \
        -R "$M2REF" -V "$OUT/tmp/$s.unfiltered.vcf.gz" \
        --mitochondria-mode \
        -O "$m2" >> "$OUT/logs/$s.mutect2.log" 2>&1 \
      || true
    if [ ! -s "$m2" ]; then
      echo; echo "  Mutect2 failed for $s — last lines of $OUT/logs/$s.mutect2.log:" >&2
      tail -8 "$OUT/logs/$s.mutect2.log" >&2; exit 1
    fi
    rm -f "$rg" "$rg.bai"
    printf " mutect2 ok"
  else
    printf " mutect2 (done)"
  fi

  # --- mutserve ---------------------------------------------------------------
  if [ ! -s "$ms" ]; then
    java -Xmx4g -jar "$MS_JAR" call \
        --reference "$R/chrM.fa" --contig-name chrM \
        --output "$ms" \
        --level 0.01 --baseQ 20 --mapQ 30 --threads "$T" \
        --no-ansi \
        "$bam" > "$OUT/logs/$s.mutserve.log" 2>&1 \
      || { echo; echo "  mutserve failed for $s — see $OUT/logs/$s.mutserve.log" >&2; tail -5 "$OUT/logs/$s.mutserve.log" >&2; exit 1; }
    printf "  mutserve ok\n"
  else
    printf "  mutserve (done)\n"
  fi
done
echo "  finished: $OUT"
