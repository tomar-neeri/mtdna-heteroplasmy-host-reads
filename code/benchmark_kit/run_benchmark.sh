#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_benchmark.sh — benchmark the in-house heteroplasmy caller against
# Mutect2 (mitochondria mode) and mutserve. Run from WSL, outside the container.
#
#   bash run_benchmark.sh check        # verify paths, files and tools (fast)
#   bash run_benchmark.sh setup        # build reference index/dictionaries
#   bash run_benchmark.sh cohort       # both tools on the 47 final BAMs
#   bash run_benchmark.sh mixtures     # build M x R mixtures, run all callers
#   bash run_benchmark.sh unfiltered   # regenerate Pass-1 BAMs, run both tools
#   bash run_benchmark.sh report       # tables, figure, markdown report
#   bash run_benchmark.sh all          # everything, in that order
#
# Every step is resumable: finished samples are skipped on re-run.
# ---------------------------------------------------------------------------
set -euo pipefail

# --- paths (override any of these on the command line: VAR=... bash run_benchmark.sh) ---
IMAGE="${IMAGE:-mtdna-pipeline:1.5.0}"
FINAL_BAMS="${FINAL_BAMS:-/mnt/d/Cancer_microbiome/RESULTS_FINAL/chrM}"     # post-2b chrM BAMs
SUMMARY="${SUMMARY:-/mnt/d/Cancer_microbiome/RESULTS_FINAL/summary}"       # pipeline summary tables
REFS="${REFS:-$HOME/ref}"                                                  # hg38.fa + bwa index
BENCH="${BENCH:-$HOME/mtdna-cohort47/benchmark}"                           # outputs; MUST be ext4
WORK="${WORK:-$HOME/mtdna-cohort47}"                                       # Snakemake work dir
RAW="${RAW:-/mnt/d/Cancer_microbiome/Input_all47}"                         # raw FASTQs
SCRATCH="${SCRATCH:-$HOME/scratch/mtdna_cohort47}"
PAIRS="${PAIRS:-425:418,414:402,389:393}"                                  # minor(M):major(R)
CORES="${CORES:-4}"
MEM="${MEM:-11g}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

hr(){ printf '%.0s─' {1..70}; echo; }
ok(){ echo "  [ok]   $1"; }
bad(){ echo "  [FAIL] $1"; exit 1; }

dk(){   # run a command in the container with the benchmark mounts
  docker run --rm --user "$(id -u):$(id -g)" \
    -v "$FINAL_BAMS":/bams:ro -v "$SUMMARY":/summary:ro -v "$REFS":/refs:ro \
    -v "$BENCH":/bench -v "$HERE/container":/bench_scripts:ro \
    --memory "$MEM" --cpus "$CORES" "$IMAGE" "$@"
}

need_summary(){
  for f in "$@"; do
    ls "$SUMMARY/$f" >/dev/null 2>&1 || bad "missing $SUMMARY/$f (set SUMMARY=... to the pipeline's results/summary)"
  done
}

step_check(){
  hr; echo "  check"; hr
  command -v docker >/dev/null || bad "docker not found"
  docker image inspect "$IMAGE" >/dev/null 2>&1 || bad "image $IMAGE not found"; ok "image $IMAGE"
  n=$(ls "$FINAL_BAMS"/*.chrM.bam 2>/dev/null | wc -l); [ "$n" -gt 0 ] || bad "no *.chrM.bam in $FINAL_BAMS"
  ok "$n final BAMs in $FINAL_BAMS"
  need_summary heteroplasmic_variants.tsv homoplasmic_variants.tsv error_model.json site_background.tsv
  ok "summary tables in $SUMMARY"
  [ -s "$REFS/hg38.fa" ] || bad "no $REFS/hg38.fa"; ok "hg38 at $REFS"
  mkdir -p "$BENCH"
  fs=$(df -T "$BENCH" | awk 'NR==2{print $2}'); [ "$fs" = ext4 ] && ok "benchmark dir on ext4" || echo "  [warn] $BENCH is on $fs"
  one=$(ls "$FINAL_BAMS"/*.chrM.bam | head -1)
  dk bash -c "samtools view -H /bams/$(basename "$one") | grep -c '^@SQ' | xargs -I{} echo '  [ok]   BAM header lists {} contigs';
              samtools view -H /bams/$(basename "$one") | grep -q '^@RG' && echo '  [info] BAM has read groups' || echo '  [info] no read groups (added on the fly for Mutect2)';
              gatk --version 2>&1 | grep -m1 -i 'gatk' | sed 's/^/  [ok]   /';
              test -s \"\$MUTSERVE_JAR\" && echo \"  [ok]   mutserve at \$MUTSERVE_JAR\" || { echo '  [FAIL] mutserve jar missing'; exit 1; };
              java -jar \"\$MUTSERVE_JAR\" call --help 2>&1 | grep -q -- '--contig-name' && echo '  [ok]   mutserve supports --contig-name' || echo '  [warn] mutserve call --help did not list --contig-name; see logs if it fails'"
}

step_setup(){ hr; echo "  setup: reference index and dictionaries"; hr; mkdir -p "$BENCH"; dk bash /bench_scripts/setup_refs.sh; }

step_cohort(){
  hr; echo "  cohort: Mutect2 + mutserve on the final BAMs"; hr
  dk bash /bench_scripts/run_tools.sh /bams /bench/cohort hg38 "$CORES"
}

step_mixtures(){
  hr; echo "  mixtures: build, then call with all four callers"; hr
  need_summary homoplasmic_variants.tsv error_model.json site_background.tsv
  dk python3 /bench_scripts/make_mixtures.py --bam-dir /bams \
      --homoplasmic /summary/homoplasmic_variants.tsv --pairs "$PAIRS" --out /bench/mixtures
  dk bash /bench_scripts/run_tools.sh /bench/mixtures/bams /bench/mixtures/tools hg38 "$CORES"
  dk bash /bench_scripts/run_ours.sh /bench/mixtures/bams /bench/mixtures/ours_pon   /summary pon
  dk bash /bench_scripts/run_ours.sh /bench/mixtures/bams /bench/mixtures/ours_nopon /summary nopon
}

step_unfiltered(){
  # Regenerate Pass-1 (rCRS-only, pre-NUMT-removal) BAMs with the pipeline's own
  # trim and recruit rules, one sample at a time so scratch never fills.
  hr; echo "  unfiltered: regenerate Pass-1 BAMs with the pipeline's rules"; hr
  [ -f "$WORK/Snakefile" ] && [ -f "$WORK/config/config.yaml" ] || bad "no Snakefile/config in $WORK"
  [ -d "$RAW" ] || bad "raw FASTQ dir not found: $RAW"
  mkdir -p "$BENCH/unfiltered_bams" "$SCRATCH"
  for bam in "$FINAL_BAMS"/*.chrM.bam; do
    s=$(basename "$bam" .chrM.bam)
    dest="$BENCH/unfiltered_bams/$s.cand.bam"
    [ -s "$dest.bai" ] && { echo "  $s (done)"; continue; }
    echo "  $s"
    docker run --rm --user "$(id -u):$(id -g)" \
      -v "$RAW":/data/raw:ro -v "$REFS":/refs:ro -v "$SCRATCH":/scratch -v "$HOME/db":/db:ro \
      -v "$WORK":/work -w /work --memory "$MEM" --cpus "$CORES" "$IMAGE" \
      snakemake -c"$CORES" --configfile config/config.yaml --notemp --rerun-triggers mtime \
                --quiet rules "/scratch/$s.cand.bam" > "$BENCH/unfiltered_bams/$s.snakemake.log" 2>&1 \
      || { tail -20 "$BENCH/unfiltered_bams/$s.snakemake.log"; bad "Pass-1 regeneration failed for $s"; }
    mv "$SCRATCH/$s.cand.bam" "$dest"; mv "$SCRATCH/$s.cand.bam.bai" "$dest.bai"
    rm -f "$SCRATCH/$s.recruit.bam" "$SCRATCH/$s.recruit.bam.bai" "$SCRATCH/$s.trim_R1.fq.gz" "$SCRATCH/$s.trim_R2.fq.gz"
  done
  docker run --rm --user "$(id -u):$(id -g)" -v "$REFS":/refs:ro -v "$BENCH":/bench -v "$HERE/container":/bench_scripts:ro \
    --memory "$MEM" --cpus "$CORES" "$IMAGE" \
    bash /bench_scripts/run_tools.sh /bench/unfiltered_bams /bench/unfiltered chrM "$CORES"
}

step_report(){
  hr; echo "  report"; hr
  args=(--ours /summary/heteroplasmic_variants.tsv --cohort /bench/cohort --out /bench/report)
  [ -d "$BENCH/unfiltered/mutect2" ] && args+=(--unfiltered /bench/unfiltered)
  if [ -f "$BENCH/mixtures/truth.tsv" ] && [ -d "$BENCH/mixtures/tools/mutect2" ]; then
    args+=(--mix /bench/mixtures --mix-tools /bench/mixtures/tools)
    [ -d "$BENCH/mixtures/ours_pon" ]   && args+=(--mix-ours /bench/mixtures/ours_pon)
    [ -d "$BENCH/mixtures/ours_nopon" ] && args+=(--mix-ours-nopon /bench/mixtures/ours_nopon)
  fi
  dk python3 /bench_scripts/compare.py "${args[@]}"
  echo "  open: $BENCH/report/benchmark_report.md"
}

case "${1:-}" in
  check) step_check ;;
  setup) step_setup ;;
  cohort) step_cohort ;;
  mixtures) step_mixtures ;;
  unfiltered) step_unfiltered ;;
  report) step_report ;;
  all) step_check; step_setup; step_cohort; step_mixtures; step_unfiltered; step_report ;;
  *) sed -n '2,16p' "$0"; exit 1 ;;
esac
