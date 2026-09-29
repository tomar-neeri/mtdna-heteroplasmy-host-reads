#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_container.sh — execute the containerized mtDNA pipeline
#
#   bash run_container.sh              # run
#   bash run_container.sh -n           # dry run
#   bash run_container.sh --kraken     # enable Pass 2b
#   bash run_container.sh -j 2         # fewer cores
#
# Handles the three things that are easy to get wrong:
#   1. Snakemake's working directory MUST be on ext4, not /mnt/d.
#      Its lock files in .snakemake/ fail on the 9p layer and it loses
#      track of running jobs (reports failure while the job is still going).
#   2. Container UID must match yours or outputs are unwritable.
#   3. Memory cap must leave headroom or bwa gets OOM-killed silently.
# ---------------------------------------------------------------------------
set -euo pipefail

# --- configuration ---------------------------------------------------------
IMAGE="mtdna-pipeline:1.5.0"
SRC="/mnt/d/Cancer_microbiome/Pipeline_container"   # where the recipe lives
WORK="$HOME/mtdna-work"                             # MUST be ext4
RAW="/mnt/d/Cancer_microbiome/Input"
REFS="$HOME/ref"
SCRATCH="$HOME/scratch/mtdna_interim"
RESULTS_COPY="${RESULTS_COPY:-/mnt/d/Cancer_microbiome/RESULTS}"   # exported here
MEM="11g"
CORES=4

SNAKE_ARGS=()
KRAKEN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --kraken) KRAKEN=1; shift ;;
    -j) CORES="$2"; shift 2 ;;
    *) SNAKE_ARGS+=("$1"); shift ;;
  esac
done

hr(){ printf '%.0s─' {1..70}; echo; }
ok(){ echo "  [ok]   $1"; }
bad(){ echo "  [FAIL] $1"; exit 1; }

hr; echo "  mtDNA pipeline — containerized run"; hr

# --- preflight -------------------------------------------------------------
command -v docker >/dev/null || bad "docker not found (enable Docker Desktop WSL integration)"
docker image inspect "$IMAGE" >/dev/null 2>&1 || bad "image $IMAGE not built — run: docker build -t $IMAGE $SRC"
ok "image $IMAGE present"

[ -d "$RAW" ] || bad "raw dir not found: $RAW"
n=$(ls "$RAW"/*_R1_001.fastq.gz 2>/dev/null | wc -l)
[ "$n" -gt 0 ] || bad "no *_R1_001.fastq.gz in $RAW"
ok "$n samples in $RAW"

[ -f "$REFS/hg38.fa.bwt" ] || bad "hg38 not indexed at $REFS/hg38.fa — run bwa index first"
ok "hg38 index present"

mkdir -p "$WORK" "$SCRATCH"
fs=$(df -T "$WORK" | awk 'NR==2{print $2}')
case "$fs" in
  ext4|xfs|btrfs) ok "work dir on $fs (correct)" ;;
  9p|drvfs) bad "work dir $WORK is on $fs — Snakemake lock files will fail there. Use a path under \$HOME." ;;
  *) echo "  [warn] work dir filesystem: $fs" ;;
esac

fs2=$(df -T "$SCRATCH" | awk 'NR==2{print $2}')
[ "$fs2" = "ext4" ] && ok "scratch on ext4" || echo "  [warn] scratch on $fs2 — alignment will be slow"

avail=$(df -BG "$SCRATCH" | awk 'NR==2{gsub("G","",$4);print $4}')
echo "  scratch free: ${avail} GB"
[ "${avail:-0}" -ge 80 ] || echo "  [warn] under 80 GB — a full cohort run peaks near there"

# --- sync the recipe into the work dir -------------------------------------
cp "$SRC/Snakefile" "$WORK/"
mkdir -p "$WORK/config" "$WORK/envs"
cp "$SRC/config/config.yaml" "$WORK/config/"
cp "$SRC/envs/mtdna.yaml"    "$WORK/envs/"
ok "recipe synced to $WORK"

# hash.k2d is memory-mapped at 7.5 GB; raise the cap whenever Pass 2b will run,
# whether enabled by flag or already true in the config file
if grep -qE '^run_kraken:[[:space:]]*true' "$SRC/config/config.yaml" 2>/dev/null; then
  KRAKEN=1
  ok "Pass 2b enabled in config"
fi

CFG_OVERRIDE=()
if [ "$KRAKEN" -eq 1 ]; then
  [ -f "$HOME/db/k2_standard_08gb/hash.k2d" ] || bad "kraken db not found at ~/db/k2_standard_08gb"
  CFG_OVERRIDE=(--config run_kraken=true kraken_db=/db/k2_standard_08gb)
  MEM="12g"    # hash.k2d is memory-mapped at 7.5 GB
  ok "Pass 2b enabled (memory raised to $MEM)"
fi

# --- run -------------------------------------------------------------------
hr; echo "  launching (${CORES} cores, ${MEM} memory cap)"; hr
cd "$WORK"

docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$RAW":/data/raw:ro \
  -v "$REFS":/refs:ro \
  -v "$SCRATCH":/scratch \
  -v "$HOME/db":/db:ro \
  -v "$WORK":/work -w /work \
  --memory "$MEM" --cpus "$CORES" \
  "$IMAGE" \
  snakemake -c"$CORES" --configfile config/config.yaml \
    "${CFG_OVERRIDE[@]}" "${SNAKE_ARGS[@]}"

# --- export results --------------------------------------------------------
if [[ ! " ${SNAKE_ARGS[*]} " =~ " -n " ]] && [ -d "$WORK/results" ]; then
  STAMP=$(date +%Y%m%d_%H%M)
  DEST="$RESULTS_COPY"
  mkdir -p "$DEST"

  hr; echo "  exporting results"; hr

  # everything except the bulky BAMs, which stay on ext4
  for sub in summary calls qc logs; do
    if [ -d "$WORK/results/$sub" ]; then
      mkdir -p "$DEST/$sub"
      cp -r "$WORK/results/$sub/." "$DEST/$sub/"
      n=$(find "$WORK/results/$sub" -type f | wc -l)
      printf "  %-10s %4s files\n" "$sub" "$n"
    fi
  done

  # BAMs are large; copy only if there is room
  if [ -d "$WORK/results/chrM" ]; then
    need=$(du -sm "$WORK/results/chrM" | cut -f1)
    free=$(df -Pm "$DEST" | awk 'NR==2{print $4}')
    if [ "$free" -gt $((need * 2)) ]; then
      mkdir -p "$DEST/chrM"; cp -r "$WORK/results/chrM/." "$DEST/chrM/"
      printf "  %-10s %4s files (%s MB)\n" "chrM" "$(ls "$WORK/results/chrM" | wc -l)" "$need"
    else
      printf "  %-10s skipped (%s MB needed, %s MB free)\n" "chrM" "$need" "$free"
    fi
  fi

  # single archive of the text results, for sharing or deposit
  ARCHIVE="$DEST/results_${STAMP}.tar.gz"
  tar -czf "$ARCHIVE" -C "$WORK/results" summary calls qc 2>/dev/null || true
  [ -f "$ARCHIVE" ] && printf "  %-10s %s\n" "archive" "$(basename "$ARCHIVE") ($(du -h "$ARCHIVE" | cut -f1))"

  # NOTE: intermediates are deliberately NOT cleaned here. Snakemake resumes by
  # checking which outputs exist; deleting them turns a restart into a full
  # reprocess. Clean manually once a run has completed:  bash cohort.sh clean
  hr; echo "  results: $DEST"; hr

  # --- headline tables -----------------------------------------------------
  show(){ [ -f "$1" ] && { echo; echo "  $(basename "$1"):"; column -t < "$1" | sed 's/^/    /' | head -${2:-15}; }; }
  show "$WORK/results/summary/cohort_overview.tsv" 20
  show "$WORK/results/summary/region_distribution.tsv" 10

  REC="$WORK/results/summary/variant_recurrence.tsv"
  if [ -f "$REC" ]; then
    n=$(awk -F'\t' 'NR>1 && $1>1' "$REC" | wc -l)
    echo
    if [ "$n" -gt 0 ]; then
      echo "  RECURRENT POSITIONS ($n) — check haplogroup linkage or artifact:"
      awk -F'\t' 'NR>1 && $1>1{printf "    %s samples  m.%s%s  median VAF %s%%  [%s]\n",$1,$2,$3,$5,$4}' "$REC"
    else
      echo "  No recurrent positions — all calls private to individuals (expected)."
    fi
  fi
  hr
fi
