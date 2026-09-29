#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# cohort.sh — run the full cohort with pause/resume
#
#   bash cohort.sh start     launch (detached, survives terminal close)
#   bash cohort.sh status    progress, current stage, ETA
#   bash cohort.sh pause     freeze instantly (machine must stay on)
#   bash cohort.sh resume    unfreeze, or restart from where it stopped
#   bash cohort.sh stop      graceful halt; safe to reboot, resume later
#   bash cohort.sh clean     remove intermediates AFTER the run has finished
#
# WHY INTERMEDIATES ARE KEPT
# Snakemake resumes by checking which outputs already exist. Deleting trimmed
# FASTQs or recruitment BAMs mid-run turns any restart into a full reprocess --
# 15 hours instead of minutes. Cleanup is therefore manual and only safe once
# the run has completed.
# ---------------------------------------------------------------------------
set -uo pipefail

SRC="${SRC:-/mnt/d/Cancer_microbiome/Pipeline_container}"
WORK="${WORK:-$HOME/mtdna-cohort47}"
RAW="${RAW:-/mnt/d/Cancer_microbiome/Input_all47}"
REFS="${REFS:-$HOME/ref}"
KDB_HOST="${KDB_HOST:-$HOME/db}"
SCRATCH="${SCRATCH:-$HOME/scratch/mtdna_cohort47}"
EXPORT_DIR="${EXPORT_DIR:-/mnt/d/Cancer_microbiome/RESULTS_3}"
IMAGE="${IMAGE:-mtdna-pipeline:1.5.1}"
MEM="${MEM:-12g}"
CORES="${CORES:-4}"
NAME="mtdna_cohort47"
LOG="$WORK/cohort.log"

hr(){ printf '%.0s─' {1..70}; echo; }
ok(){ echo "  [ok]   $1"; }
warn(){ echo "  [warn] $1"; }
bad(){ echo "  [FAIL] $1"; }

running(){ [ -n "$(docker ps -q -f name=^${NAME}$ 2>/dev/null)" ]; }
paused(){ [ "$(docker inspect -f '{{.State.Status}}' "$NAME" 2>/dev/null)" = "paused" ]; }
exists(){ [ -n "$(docker ps -aq -f name=^${NAME}$ 2>/dev/null)" ]; }

n_samples(){ ls "$RAW"/*_R1_001.fastq.gz 2>/dev/null | wc -l; }

# ---------------------------------------------------------------------------
setup(){
  mkdir -p "$WORK"/{config,envs} "$SCRATCH" "$EXPORT_DIR"
  cp -f "$SRC/Snakefile" "$WORK/"
  cp -f "$SRC/config/config.yaml" "$WORK/config/"
  cp -f "$SRC/envs/mtdna.yaml" "$WORK/envs/" 2>/dev/null
}

cmd_start(){
  hr; echo "  cohort run — start"; hr
  if running; then warn "already running; use 'status'"; return 1; fi
  if paused;  then warn "paused; use 'resume'"; return 1; fi

  local n; n=$(n_samples)
  [ "$n" -gt 0 ] || { bad "no *_R1_001.fastq.gz in $RAW"; return 1; }

  # every R1 needs its mate, or Snakemake fails partway through
  local miss=0
  for f in "$RAW"/*_R1_001.fastq.gz; do
    [ -f "${f/_R1_/_R2_}" ] || { bad "no mate for $(basename "$f")"; miss=1; }
  done
  [ "$miss" -eq 1 ] && return 1

  mkdir -p "$WORK" "$SCRATCH"
  [ -f "$REFS/hg38.fa.bwt" ] || { bad "hg38 not indexed at $REFS"; return 1; }
  local fs; fs=$(df -T "$WORK" | awk 'NR==2{print $2}')
  case "$fs" in ext4|xfs|btrfs) ;; *) bad "work dir on $fs — Snakemake locks fail there"; return 1 ;; esac

  local free; free=$(df -BG "$SCRATCH" | awk 'NR==2{gsub("G","",$4);print $4}')
  ok "$n samples | scratch ${free} GB free | work on $fs"
  [ "${free:-0}" -lt 150 ] && warn "under 150 GB — a 47-sample run peaks near there"

  setup
  exists && docker rm -f "$NAME" >/dev/null 2>&1

  # clear a stale lock from a previous kill
  docker run --rm --user "$(id -u):$(id -g)" \
    -v "$RAW":/data/raw:ro -v "$REFS":/refs:ro -v "$WORK":/work -w /work \
    "$IMAGE" snakemake --configfile config/config.yaml --unlock >/dev/null 2>&1

  echo
  ok "launching (${CORES} cores, ${MEM} cap) — detached as '$NAME'"
  docker run -d --name "$NAME" \
    --user "$(id -u):$(id -g)" \
    -v "$RAW":/data/raw:ro \
    -v "$REFS":/refs:ro \
    -v "$KDB_HOST":/db:ro \
    -v "$SCRATCH":/scratch \
    -v "$WORK":/work -w /work \
    --memory "$MEM" --cpus "$CORES" \
    "$IMAGE" \
    snakemake -c"$CORES" --configfile config/config.yaml \
              --rerun-incomplete --keep-going >/dev/null

  sleep 3
  running && ok "started" || bad "failed to start — docker logs $NAME"
  hr
  echo "  watch:  bash cohort.sh status"
  echo "  pause:  bash cohort.sh pause"
  hr
}

# ---------------------------------------------------------------------------
cmd_status(){
  hr; echo "  cohort run — status"; hr
  if paused;      then echo; warn "PAUSED (frozen; 'resume' to continue)"
  elif running;   then echo; ok "RUNNING"
  elif exists;    then echo; warn "STOPPED (use 'resume' to continue from cache)"
  else            echo; warn "not started"; fi

  local n; n=$(n_samples)
  echo
  printf "  %-24s %s\n" "samples"        "$n"
  printf "  %-24s %s\n" "chrM BAMs done" "$(ls "$WORK"/results/chrM/*.bam 2>/dev/null | wc -l)/$n"
  printf "  %-24s %s\n" "calls done"     "$(ls "$WORK"/results/calls/*.calls.tsv 2>/dev/null | wc -l)/$n"
  for f in site_background.tsv error_model.json heteroplasmy_summary.tsv; do
    [ -f "$WORK/results/summary/$f" ] && printf "  %-24s %s\n" "$f" "done" \
                                      || printf "  %-24s %s\n" "$f" "pending"
  done

  echo
  local done total
  done=$(docker logs "$NAME" 2>&1 | grep -c "Finished job"); done=${done:-0}
  total=$(docker logs "$NAME" 2>&1 | grep -oP '^total\s+\K[0-9]+' | tail -1)
  if [ -n "${total:-}" ] && [ "$total" -gt 0 ]; then
    printf "  progress  %s/%s jobs (%d%%)\n" "$done" "$total" $((done*100/total))
  fi
  echo
  echo "  last activity:"
  docker logs --tail 4 "$NAME" 2>&1 | sed 's/^/    /'
  echo
  printf "  %-24s %s\n" "scratch" "$(du -sh "$SCRATCH" 2>/dev/null | cut -f1)"
  printf "  %-24s %s\n" "free"    "$(df -h "$SCRATCH" | awk 'NR==2{print $4}')"
  hr
}

# ---------------------------------------------------------------------------
cmd_pause(){
  hr; echo "  pause"; hr
  running || { warn "not running"; return 1; }
  docker pause "$NAME" >/dev/null && ok "frozen — all processes suspended"
  echo
  echo "  CPU and I/O released; memory still held."
  echo "  Survives terminal close, NOT a reboot or Docker restart."
  echo "  For a reboot use 'stop' instead."
  hr
}

cmd_resume(){
  hr; echo "  resume"; hr
  if paused; then
    docker unpause "$NAME" >/dev/null && ok "unfrozen, continuing"
  elif running; then
    warn "already running"
  elif exists; then
    ok "restarting from cache — completed outputs are skipped"
    docker rm -f "$NAME" >/dev/null 2>&1
    cmd_start
  else
    warn "no previous run; use 'start'"
  fi
  hr
}

cmd_stop(){
  hr; echo "  stop"; hr
  running || paused || { warn "not running"; return 1; }
  paused && docker unpause "$NAME" >/dev/null 2>&1
  ok "sending stop (current jobs finish, up to 60s)"
  docker stop -t 60 "$NAME" >/dev/null
  docker run --rm --user "$(id -u):$(id -g)" \
    -v "$RAW":/data/raw:ro -v "$REFS":/refs:ro -v "$WORK":/work -w /work \
    "$IMAGE" snakemake --configfile config/config.yaml --unlock >/dev/null 2>&1
  ok "stopped and unlocked"
  echo
  echo "  Intermediates kept in $SCRATCH — 'resume' continues from here."
  echo "  Safe to reboot."
  hr
}

# ---------------------------------------------------------------------------
cmd_clean(){
  hr; echo "  clean intermediates"; hr
  running || paused && { bad "run is active — stop it first"; return 1; }
  [ -f "$WORK/results/summary/heteroplasmy_summary.tsv" ] || {
    bad "run has not completed; cleaning now would force a full reprocess"; return 1; }
  local b a
  b=$(df -Pm "$SCRATCH" | awk 'NR==2{print $4}')
  rm -rf "${SCRATCH:?}"/*
  a=$(df -Pm "$SCRATCH" | awk 'NR==2{print $4}')
  ok "reclaimed $(( (a-b)/1024 )) GB"
  hr
}

cmd_export(){
  hr; echo "  export results -> $EXPORT_DIR"; hr
  [ -d "$WORK/results" ] || { bad "no results"; return 1; }
  mkdir -p "$EXPORT_DIR"
  for sub in summary calls qc logs chrM; do
    [ -d "$WORK/results/$sub" ] && { mkdir -p "$EXPORT_DIR/$sub"
      cp -r "$WORK/results/$sub/." "$EXPORT_DIR/$sub/"
      printf "  %-10s %4s files\n" "$sub" "$(find "$WORK/results/$sub" -type f | wc -l)"; }
  done
  docker inspect --format='{{.Id}}' "$IMAGE" > "$EXPORT_DIR/IMAGE_ID.txt" 2>/dev/null
  cp "$WORK/config/config.yaml" "$WORK/Snakefile" "$EXPORT_DIR/" 2>/dev/null
  docker run --rm "$IMAGE" cat /opt/package_versions.txt > "$EXPORT_DIR/package_versions.txt" 2>/dev/null
  tar -czf "$EXPORT_DIR/results_$(date +%Y%m%d_%H%M).tar.gz" \
      -C "$WORK/results" summary calls qc 2>/dev/null
  ok "exported with image id and config"
  hr
}

case "${1:-status}" in
  start)  cmd_start ;;
  status) cmd_status ;;
  pause)  cmd_pause ;;
  resume) cmd_resume ;;
  stop)   cmd_stop ;;
  clean)  cmd_clean ;;
  export) cmd_export ;;
  *) echo "usage: bash cohort.sh [start|status|pause|resume|stop|clean|export]" ;;
esac
