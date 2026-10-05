#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# numt_ablation.sh — cohort-wide NUMT ablation (no Pass 2a)
#
#   bash numt_ablation.sh start     launch detached
#   bash numt_ablation.sh status    progress
#   bash numt_ablation.sh stop      graceful halt; resume later
#   bash numt_ablation.sh resume    continue from cache
#   bash numt_ablation.sh report    print the comparison table
#   bash numt_ablation.sh clean     remove intermediates, AFTER completion
#   bash numt_ablation.sh export    copy results to EXPORT_DIR
#
# SAFETY
# Separate work dir, scratch dir and container name from the main run. The
# published results in ~/mtdna-cohort47 are mounted READ-ONLY and are never
# written to. Nothing here can overwrite them.
#
# Requires a completed baseline run: reuses its error_model.json,
# site_background.tsv, heteroplasmy_lineage_flags.tsv and
# variants_annotated.tsv.
#
# No network access is needed at any point.
# ---------------------------------------------------------------------------
set -uo pipefail

SRC="${SRC:-/mnt/d/Cancer_microbiome/Pipeline_container}"
BASELINE="${BASELINE:-$HOME/mtdna-cohort47}"
WORK="${WORK:-$HOME/mtdna-numt-ablation}"
RAW="${RAW:-/mnt/d/Cancer_microbiome/Input_all47}"
REFS="${REFS:-$HOME/ref}"
SCRATCH="${SCRATCH:-$HOME/scratch/mtdna_numt}"
EXPORT_DIR="${EXPORT_DIR:-/mnt/d/Cancer_microbiome/RESULTS_NUMT_ABLATION}"
IMAGE="${IMAGE:-mtdna-pipeline:1.5.0}"
MEM="${MEM:-12g}"
CORES="${CORES:-4}"
NAME="mtdna_numt_ablation"

hr(){ printf '%.0s─' {1..70}; echo; }
ok(){ echo "  [ok]   $1"; }
warn(){ echo "  [warn] $1"; }
bad(){ echo "  [FAIL] $1"; }

running(){ [ -n "$(docker ps -q -f name=^${NAME}$ 2>/dev/null)" ]; }
exists(){ [ -n "$(docker ps -aq -f name=^${NAME}$ 2>/dev/null)" ]; }
n_samples(){ ls "$RAW"/*_R1_001.fastq.gz 2>/dev/null | wc -l; }

# Mount and resource flags, shared by every invocation. Kept as an array so
# docker's own flags, the image, and the container command stay separate --
# collapsing them into one "$@" makes docker read the command as the image.
mounts(){
  MOUNTS=(
    --user "$(id -u):$(id -g)"
    -v "$RAW":/data/raw:ro
    -v "$REFS":/refs:ro
    -v "$BASELINE":/baseline:ro
    -v "$SCRATCH":/scratch
    -v "$WORK":/work
    -w /work
    --memory "$MEM"
    --cpus "$CORES"
  )
}

SNAKE_ARGS=(-s Snakefile.numt_ablation --configfile config/config.yaml)

setup(){
  mkdir -p "$WORK"/{config,scripts} "$SCRATCH" "$EXPORT_DIR"
  cp -f "$SRC/Snakefile.numt_ablation" "$WORK/"
  cp -f "$SRC/scripts/compare_numt.py" "$WORK/scripts/"
  cp -f "$SRC/config/config.yaml"      "$WORK/config/"
}

unlock(){
  mounts
  docker run --rm "${MOUNTS[@]}" "$IMAGE" \
    snakemake "${SNAKE_ARGS[@]}" --unlock >/dev/null 2>&1
}

cmd_start(){
  hr; echo "  NUMT ablation — start"; hr
  running && { warn "already running; use 'status'"; return 1; }

  local n; n=$(n_samples)
  [ "$n" -gt 0 ] || { bad "no *_R1_001.fastq.gz in $RAW"; return 1; }

  local missing=0
  for f in summary/error_model.json summary/site_background.tsv \
           summary/heteroplasmy_lineage_flags.tsv summary/variants_annotated.tsv; do
    [ -f "$BASELINE/results/$f" ] || { bad "baseline missing: results/$f"; missing=1; }
  done
  [ "$missing" -eq 1 ] && { echo; echo "  The main run must finish first."; return 1; }

  [ -f "$SRC/Snakefile.numt_ablation" ] || { bad "Snakefile.numt_ablation not in $SRC"; return 1; }
  [ -f "$SRC/scripts/compare_numt.py" ] || { bad "compare_numt.py not in $SRC/scripts"; return 1; }

  # the gnomAD table drives the NUMT false-positive overlap, the most
  # informative output. Without it those columns silently read zero.
  [ -f "$REFS/gnomad_mt.tsv" ] || warn "no $REFS/gnomad_mt.tsv — gnomAD overlap columns will be empty"

  # create the dirs BEFORE measuring free space, or df reports nothing and
  # the headroom check fires spuriously
  mkdir -p "$WORK" "$SCRATCH"

  local fs; fs=$(df -T "$WORK" | awk 'NR==2{print $2}')
  case "$fs" in ext4|xfs|btrfs) ;; *) bad "work dir on $fs — Snakemake locks fail there"; return 1 ;; esac

  local free; free=$(df -BG "$SCRATCH" 2>/dev/null | awk 'NR==2{gsub("G","",$4);print $4}')
  ok "$n samples | scratch ${free:-?} GB free | work on $fs"
  ok "baseline (read-only): $BASELINE"
  if [ -n "${free:-}" ] && [ "$free" -lt 60 ]; then
    warn "under 60 GB — trimmed FASTQs peak near 40 GB with $CORES jobs"
  fi
  ok "hg38 not required for this run"

  setup
  exists && docker rm -f "$NAME" >/dev/null 2>&1
  unlock

  mounts
  echo
  ok "launching (${CORES} cores, ${MEM} cap) — detached as '$NAME'"
  docker run -d --name "$NAME" "${MOUNTS[@]}" "$IMAGE" \
    snakemake "${SNAKE_ARGS[@]}" -c"$CORES" \
      --config scratch_dir=/scratch out_dir=results \
               baseline_summary=/baseline/results/summary \
      --rerun-incomplete --keep-going >/dev/null

  sleep 4
  if running; then
    ok "started"
  else
    bad "failed to start — last output:"
    docker logs "$NAME" 2>&1 | tail -15 | sed 's/^/    /'
    return 1
  fi
  hr
  echo "  watch:   bash numt_ablation.sh status"
  echo "  results: bash numt_ablation.sh report"
  hr
}

cmd_status(){
  hr; echo "  NUMT ablation — status"; hr
  if running;   then echo; ok "RUNNING"
  elif exists;  then echo; warn "STOPPED (use 'resume')"
  else          echo; warn "not started"; fi

  local n; n=$(n_samples)
  echo
  printf "  %-28s %s\n" "samples"            "$n"
  printf "  %-28s %s\n" "no-NUMT chrM BAMs"  "$(ls "$WORK"/results/chrM_nonumt/*.bam 2>/dev/null | wc -l)/$n"
  printf "  %-28s %s\n" "calls (refit)"      "$(ls "$WORK"/results/calls_refit/*.calls.tsv 2>/dev/null | wc -l)/$n"
  printf "  %-28s %s\n" "calls (fixed)"      "$(ls "$WORK"/results/calls_fixed/*.calls.tsv 2>/dev/null | wc -l)/$n"
  for f in refit/error_model.json refit/site_background.tsv \
           refit/variants_annotated.tsv fixed/variants_annotated.tsv \
           numt_ablation_comparison.tsv; do
    [ -f "$WORK/results/summary/$f" ] && printf "  %-28s %s\n" "$f" "done" \
                                      || printf "  %-28s %s\n" "$f" "pending"
  done

  echo
  if exists; then
    local done total
    done=$(docker logs "$NAME" 2>&1 | grep -c "Finished job")
    total=$(docker logs "$NAME" 2>&1 | grep -oP '^total\s+\K[0-9]+' | tail -1)
    if [ -n "${total:-}" ] && [ "$total" -gt 0 ]; then
      printf "  progress  %s/%s jobs (%d%%)\n" "$done" "$total" $((done*100/total))
    fi
    echo
    echo "  last activity:"
    docker logs --tail 5 "$NAME" 2>&1 | sed 's/^/    /'
  fi
  echo
  printf "  %-28s %s\n" "scratch" "$(du -sh "$SCRATCH" 2>/dev/null | cut -f1)"
  printf "  %-28s %s\n" "free"    "$(df -h "$SCRATCH" 2>/dev/null | awk 'NR==2{print $4}')"
  hr
}

cmd_stop(){
  hr; echo "  stop"; hr
  running || { warn "not running"; return 1; }
  ok "sending stop (current jobs finish, up to 60s)"
  docker stop -t 60 "$NAME" >/dev/null
  unlock
  ok "stopped and unlocked — 'resume' continues from cache"
  hr
}

cmd_resume(){
  hr; echo "  resume"; hr
  running && { warn "already running"; return 1; }
  exists || { warn "no previous run; use 'start'"; return 1; }
  docker rm -f "$NAME" >/dev/null 2>&1
  cmd_start
}

cmd_report(){
  local f="$WORK/results/summary/numt_ablation_comparison.tsv"
  [ -f "$f" ] || { bad "not finished — no $f yet"; return 1; }
  hr; echo "  NUMT ablation — comparison"; hr
  column -t -s $'\t' "$f" 2>/dev/null || cat "$f"
  hr
  echo "  per-sample: $WORK/results/summary/numt_ablation_per_sample.tsv"
  hr
}

cmd_clean(){
  hr; echo "  clean intermediates"; hr
  running && { bad "run is active — stop it first"; return 1; }
  [ -f "$WORK/results/summary/numt_ablation_comparison.tsv" ] || {
    bad "not completed; cleaning now would force a full reprocess"; return 1; }
  rm -rf "${SCRATCH:?}"/*
  ok "scratch cleared (chrM_nonumt BAMs and all results kept)"
  hr
}

cmd_export(){
  hr; echo "  export -> $EXPORT_DIR"; hr
  [ -d "$WORK/results" ] || { bad "no results"; return 1; }
  mkdir -p "$EXPORT_DIR"
  for sub in summary calls_refit calls_fixed qc logs; do
    [ -d "$WORK/results/$sub" ] && { mkdir -p "$EXPORT_DIR/$sub"
      cp -r "$WORK/results/$sub/." "$EXPORT_DIR/$sub/"
      printf "  %-14s %4s files\n" "$sub" "$(find "$WORK/results/$sub" -type f | wc -l)"; }
  done
  cp "$WORK/config/config.yaml" "$EXPORT_DIR/" 2>/dev/null
  docker inspect --format='{{.Id}}' "$IMAGE" > "$EXPORT_DIR/IMAGE_ID.txt" 2>/dev/null
  ok "exported"
  hr
}

case "${1:-status}" in
  start)  cmd_start ;;
  status) cmd_status ;;
  stop)   cmd_stop ;;
  resume) cmd_resume ;;
  report) cmd_report ;;
  clean)  cmd_clean ;;
  export) cmd_export ;;
  *) echo "usage: bash numt_ablation.sh [start|status|stop|resume|report|clean|export]" ;;
esac
