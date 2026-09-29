#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# manage.sh — update the pipeline and keep disk usage under control
#
#   bash manage.sh status              what is installed, what disk is used
#   bash manage.sh update              stage new files, back up old, rebuild
#   bash manage.sh clean               reclaim scratch, old images, old backups
#   bash manage.sh archive             snapshot current results, then clean
#
# Files are picked up from STAGE (default: the project root plus Windows
# Downloads). Download into either and run `update`; it sorts them into the
# right subfolders, versions what it replaces, and rebuilds the image.
# ---------------------------------------------------------------------------
set -uo pipefail

SRC="${SRC:-/mnt/d/Cancer_microbiome/Pipeline_container}"
WORK="${WORK:-$HOME/mtdna-cohort47}"
SCRATCH="${SCRATCH:-$HOME/scratch/mtdna_cohort47}"
ARCHIVE_DIR="${ARCHIVE_DIR:-/mnt/d/Cancer_microbiome/ARCHIVE}"
STAGE="${STAGE:-$SRC}"
WIN_DL="${WIN_DL:-/mnt/c/Users/$(powershell.exe -NoProfile -Command '$env:USERNAME' 2>/dev/null | tr -d '\r\n')/Downloads}"
KEEP_BACKUPS="${KEEP_BACKUPS:-5}"
KEEP_IMAGES="${KEEP_IMAGES:-2}"

SCRIPTS=(pooled_error_model.py call_heteroplasmy.py collect_variants.py
         kraken_filter.py site_background.py
         haplogroup.py lineage_matrix.py)
ROOT_FILES=(Dockerfile Snakefile run_container.sh manage.sh)

hr(){ printf '%.0s─' {1..70}; echo; }
ok(){ echo "  [ok]   $1"; }
warn(){ echo "  [warn] $1"; }
bad(){ echo "  [FAIL] $1"; }
sz(){ du -sh "$1" 2>/dev/null | cut -f1 || echo "-"; }

version_of(){ grep -oP 'image\.version="\K[^"]+' "$SRC/Dockerfile" 2>/dev/null || echo "?"; }

# ---------------------------------------------------------------------------
cmd_status(){
  hr; echo "  pipeline status"; hr
  echo
  echo "Installed version : $(version_of)"
  echo "Project           : $SRC"
  echo "Work dir          : $WORK"
  echo
  echo "Scripts"
  for f in "${SCRIPTS[@]}"; do
    [ -f "$SRC/scripts/$f" ] && ok "$f  ($(date -r "$SRC/scripts/$f" +%Y-%m-%d\ %H:%M))" \
                             || bad "$f MISSING"
  done
  echo
  echo "Docker images"
  docker images --format '  {{.Repository}}:{{.Tag}}  {{.Size}}  {{.CreatedSince}}' \
    2>/dev/null | grep mtdna || warn "none built"
  echo
  echo "Disk"
  printf "  %-26s %8s\n" "scratch"        "$(sz "$SCRATCH")"
  printf "  %-26s %8s\n" "work results"   "$(sz "$WORK/results")"
  printf "  %-26s %8s\n" "project"        "$(sz "$SRC")"
  printf "  %-26s %8s\n" "backups"        "$(sz "$SRC/versions")"
  printf "  %-26s %8s\n" "archive"        "$(sz "$ARCHIVE_DIR")"
  echo
  df -h "$SCRATCH" | awk 'NR==1||NR==2{printf "  %s\n",$0}'
  hr
}

# ---------------------------------------------------------------------------
cmd_update(){
  hr; echo "  update pipeline"; hr
  local stamp found=0
  stamp=$(date +%Y%m%d_%H%M)
  mkdir -p "$SRC"/{scripts,config,envs,versions/"$stamp"}

  echo
  echo "Scanning for new files..."
  for dir in "$STAGE" "$WIN_DL"; do
    [ -d "$dir" ] || continue
    for f in "${SCRIPTS[@]}" "${ROOT_FILES[@]}" config.yaml mtdna.yaml; do
      local new="$dir/$f"
      [ -f "$new" ] || continue
      # destination by file type
      local dest
      case "$f" in
        *.py)        dest="$SRC/scripts/$f" ;;
        config.yaml) dest="$SRC/config/$f" ;;
        mtdna.yaml)  dest="$SRC/envs/$f" ;;
        *)           dest="$SRC/$f" ;;
      esac
      # skip if this IS the destination and nothing newer exists
      [ "$new" = "$dest" ] && continue
      if [ -f "$dest" ] && cmp -s "$new" "$dest"; then
        rm -f "$new"; continue
      fi
      [ -f "$dest" ] && cp "$dest" "$SRC/versions/$stamp/"
      mv -f "$new" "$dest"
      sed -i 's/\r$//' "$dest" 2>/dev/null
      ok "$f -> ${dest#$SRC/}"
      found=$((found+1))
    done
  done

  # files downloaded straight into the project root still need sorting
  for f in "${SCRIPTS[@]}"; do
    [ -f "$SRC/$f" ] && { [ -f "$SRC/scripts/$f" ] && cp "$SRC/scripts/$f" "$SRC/versions/$stamp/"
                          mv -f "$SRC/$f" "$SRC/scripts/$f"; sed -i 's/\r$//' "$SRC/scripts/$f"
                          ok "$f -> scripts/"; found=$((found+1)); }
  done
  [ -f "$SRC/config.yaml" ] && { cp "$SRC/config/config.yaml" "$SRC/versions/$stamp/" 2>/dev/null
                                 mv -f "$SRC/config.yaml" "$SRC/config/"; sed -i 's/\r$//' "$SRC/config/config.yaml"
                                 ok "config.yaml -> config/"; found=$((found+1)); }

  if [ "$found" -eq 0 ]; then
    rmdir "$SRC/versions/$stamp" 2>/dev/null
    warn "no new files found in $STAGE or $WIN_DL"
  else
    echo; ok "$found file(s) updated, previous versions in versions/$stamp"
  fi

  # --- verify completeness before building ---
  echo; echo "Verifying..."
  local missing=0
  for f in "${SCRIPTS[@]}"; do
    [ -f "$SRC/scripts/$f" ] || { bad "scripts/$f missing"; missing=1; }
  done
  [ -f "$SRC/refs/rCRS.fasta" ] || { bad "refs/rCRS.fasta missing"; missing=1; }
  [ "$missing" -eq 1 ] && { bad "cannot build"; return 1; }
  ok "all required files present"

  # --- build ---
  local ver; ver=$(version_of)
  echo; echo "Building mtdna-pipeline:$ver ..."
  if docker build -t "mtdna-pipeline:$ver" "$SRC" >/dev/null 2>&1; then
    ok "image built"
  else
    bad "build failed — rerun without suppression:  docker build -t mtdna-pipeline:$ver $SRC"
    return 1
  fi
  docker inspect --format='  image id: {{.Id}}' "mtdna-pipeline:$ver"

  # --- sync into the work dir ---
  mkdir -p "$WORK"/{config,envs}
  cp -f "$SRC/Snakefile" "$SRC/run_container.sh" "$WORK/"
  cp -f "$SRC/config/config.yaml" "$WORK/config/"
  cp -f "$SRC/envs/mtdna.yaml" "$WORK/envs/" 2>/dev/null
  sed -i "s|mtdna-pipeline:[0-9.]*|mtdna-pipeline:$ver|" "$WORK/run_container.sh"
  ok "synced to $WORK (image pinned to $ver)"
  hr
  echo "  next:  cd $WORK && nohup bash run_container.sh > run_\$(date +%Y%m%d_%H%M).log 2>&1 &"
  hr
}

# ---------------------------------------------------------------------------
cmd_clean(){
  hr; echo "  reclaim disk"; hr
  local before after
  before=$(df -Pm "$SCRATCH" | awk 'NR==2{print $4}')

  echo
  echo "Scratch intermediates ($(sz "$SCRATCH"))"
  # trimmed FASTQs and recruitment BAMs are regenerable; chrM BAMs are not
  find "$SCRATCH" -maxdepth 1 -type f \
       \( -name '*.trim_R[12].fq.gz' -o -name '*.recruit.bam*' \
          -o -name '*.hverify.bam*' -o -name '*.cand*' -o -name '*.ns.bam' \) \
       -delete 2>/dev/null
  ok "regenerable intermediates removed"

  echo
  echo "Old backups (keeping $KEEP_BACKUPS)"
  if [ -d "$SRC/versions" ]; then
    ls -1dt "$SRC/versions"/*/ 2>/dev/null | tail -n +$((KEEP_BACKUPS+1)) \
      | while read -r d; do rm -rf "$d"; echo "    removed $(basename "$d")"; done
  fi
  ok "done"

  echo
  echo "Old docker images (keeping $KEEP_IMAGES newest)"
  docker images --format '{{.Repository}}:{{.Tag}}\t{{.CreatedAt}}' 2>/dev/null \
    | grep '^mtdna-pipeline' | sort -k2 -r | tail -n +$((KEEP_IMAGES+1)) \
    | cut -f1 | while read -r img; do
        docker rmi "$img" >/dev/null 2>&1 && echo "    removed $img"
      done
  docker image prune -f >/dev/null 2>&1
  ok "done"

  after=$(df -Pm "$SCRATCH" | awk 'NR==2{print $4}')
  hr
  echo "  reclaimed $(( (after - before) / 1024 )) GB   (free now: $((after/1024)) GB)"
  hr
}

# ---------------------------------------------------------------------------
cmd_archive(){
  hr; echo "  archive results"; hr
  [ -d "$WORK/results" ] || { bad "no results in $WORK"; return 1; }
  local ver stamp dest
  ver=$(version_of); stamp=$(date +%Y%m%d_%H%M)
  dest="$ARCHIVE_DIR/${stamp}_v${ver}"
  mkdir -p "$dest"

  for sub in summary calls qc logs; do
    [ -d "$WORK/results/$sub" ] && cp -r "$WORK/results/$sub" "$dest/"
  done
  # record exactly what produced these results
  docker inspect --format='{{.Id}}' "mtdna-pipeline:$ver" > "$dest/IMAGE_ID.txt" 2>/dev/null
  cp "$SRC/config/config.yaml" "$dest/" 2>/dev/null
  tar -czf "$dest.tar.gz" -C "$ARCHIVE_DIR" "$(basename "$dest")" 2>/dev/null \
    && rm -rf "$dest" && ok "archived: $dest.tar.gz ($(sz "$dest.tar.gz"))"

  echo
  ls -1t "$ARCHIVE_DIR"/*.tar.gz 2>/dev/null | head -5 \
    | while read -r f; do printf "    %-44s %8s\n" "$(basename "$f")" "$(sz "$f")"; done
  hr
  cmd_clean
}

case "${1:-status}" in
  status)  cmd_status ;;
  update)  cmd_update ;;
  clean)   cmd_clean ;;
  archive) cmd_archive ;;
  *) echo "usage: bash manage.sh [status|update|clean|archive]"
     echo "  status   what is installed and what disk is used"
     echo "  update   stage new files, version the old, rebuild, sync"
     echo "  clean    remove regenerable intermediates, old images, old backups"
     echo "  archive  snapshot results with the image id, then clean" ;;
esac
