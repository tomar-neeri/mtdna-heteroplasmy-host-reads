#!/usr/bin/env bash
# ---------------------------------------------------------------------------
#
# The primary analysis assigns macrohaplogroups with haplogroup.py, which runs
# offline from an embedded marker set. This answers a different question: how
# deep does classification resolve when the complete PhyloTree reference is
# available, and do the marker-based assignments hold up against it?
#
# Input is the PILEUP, not the heteroplasmy caller. Calling-stage output holds
# only variants that passed a filter tuned for low-VAF detection — roughly half
# the fixed differences from the rCRS. HaploGrep3's quality score is essentially
# the fraction of expected defining markers found, so an incomplete variant list
# depresses it regardless of data quality.
#
# Self-contained: the PhyloTree reference is fetched once into a host cache
# directory and reused offline thereafter. HaploGrep3 ships no trees/ directory,
# only a repository URL, and fails config loading outright without that cache —
# including `haplogrep3 --help`.
#
#   bash code/subclade_classify.sh                      # defaults below
#   TREE=phylotree-fu-rcrs@1.2 bash code/subclade_classify.sh
#
# Outputs, in results/haplogroup_deep/:
#   varlist.tsv   sample -> fixed difference from the rCRS
#   cohort.hsd    classifier input
#   deep_<tree>.txt   HaploGrep3 extended report
#
# These are per-participant haplotypes.
# ---------------------------------------------------------------------------
set -euo pipefail

WORK="${WORK:-$PWD}"
IMAGE="${IMAGE:-mtdna-pipeline:1.5.1}"
TREE="${TREE:-phylotree-rcrs@17.2}"
TREECACHE="${TREECACHE:-$HOME/.cache/haplogrep3-trees}"
REF="${REF:-/opt/refs/rCRS.fasta}"
HG3="/opt/tools/haplogrep3"
BAMDIR="$WORK/results/chrM"
OUT="$WORK/results/haplogroup_deep"

MIN_AF="${MIN_AF:-0.95}"      # fixed difference from the rCRS
MIN_DP="${MIN_DP:-10}"
TAG="${TREE//[^A-Za-z0-9._-]/_}"

[ -d "$BAMDIR" ] || { echo "ERROR: $BAMDIR not found - run from the work directory"; exit 1; }
mkdir -p "$OUT" "$TREECACHE"

# --- 0. tree, fetched once --------------------------------------------------
if [ "${1:-}" = "trees" ]; then
  docker run --rm --network host --user root -v "$TREECACHE:$HG3/trees" \
    -w "$HG3" "$IMAGE" ./haplogrep3 trees
  exit 0
fi

if [ ! -s "$TREECACHE/.installed-$TAG" ]; then
  echo "[0/3] fetching $TREE into $TREECACHE (once; needs network)"
  docker run --rm --network host --user root -v "$TREECACHE:$HG3/trees" \
    -w "$HG3" "$IMAGE" ./haplogrep3 install-tree "$TREE" \
    && [ -n "$(ls -A "$TREECACHE")" ] \
    && touch "$TREECACHE/.installed-$TAG" \
    || { echo "tree install failed - list valid names with: bash code/subclade_classify.sh trees"; exit 1; }
fi

# --- 1. fixed differences from the rCRS ------------------------------------
echo "[1/3] calling fixed differences from BAMs (AF>=$MIN_AF, DP>=$MIN_DP)"
docker run --rm -v "$WORK:/work" "$IMAGE" bash -c '
set -euo pipefail
REF='"$REF"'; MIN_AF='"$MIN_AF"'; MIN_DP='"$MIN_DP"'
OUT=/work/results/haplogroup_deep
: > "$OUT/varlist.tsv"; : > "$OUT/pileup.log"
for bam in /work/results/chrM/*.chrM.bam; do
  sid=$(basename "$bam" .chrM.bam)
  echo "=== $sid" >> "$OUT/pileup.log"
  # -d 2000: the default cap of 250 truncates high-depth sites and biases the
  #          allele fraction; this cohort reaches 3029x.
  # -B     : BAQ suppresses real variants on a short high-copy target.
  bcftools mpileup -f "$REF" -d 2000 -q 20 -Q 20 -B -a AD,DP -Ou "$bam" 2>>"$OUT/pileup.log" \
    | bcftools call -m -Ou 2>>"$OUT/pileup.log" \
    | bcftools norm -f "$REF" -Ou 2>>"$OUT/pileup.log" \
    | bcftools query -f "%POS\t%REF\t%ALT\t[%AD]\n" -i "TYPE=\"snp\"" \
    | awk -v s="$sid" -v af="$MIN_AF" -v dp="$MIN_DP" -F"\t" "
        { split(\$4, a, \",\"); ref=a[1]+0; alt=a[2]+0; tot=ref+alt;
          if (\$1 == 3107) next;        # N placeholder in the rCRS
          if (tot >= dp && alt/tot >= af) print s \"\t\" \$1 \$3 }" >> "$OUT/varlist.tsv"
done
wc -l < "$OUT/varlist.tsv" | xargs echo "  variant rows:"
'

# --- 2. HSD -----------------------------------------------------------------
echo "[2/3] building HSD"
python3 - "$OUT/varlist.tsv" "$OUT/cohort.hsd" <<'PY'
import sys
from collections import defaultdict
d = defaultdict(list)
for line in open(sys.argv[1]):
    if line.strip():
        sid, v = line.rstrip("\n").split("\t")
        d[sid].append(v)
with open(sys.argv[2], "w") as fh:
    fh.write("SampleId\tRange\tHaplogroup\tPolymorphisms\n")
    for sid in sorted(d):
        fh.write(f"{sid}\t1-16569\t?\t{' '.join(d[sid])}\n")
n = sorted(len(v) for v in d.values())
# No minimum-variant filter. A sample whose haplogroup sits close to the rCRS
# (H2a2a1) legitimately carries few differences, so filtering on count would
# discriminate against exactly the lineages nearest the reference.
print(f"  {len(d)} samples; variants/sample mean {sum(n)/len(n):.1f} "
      f"min {n[0]} max {n[-1]}")
PY

# --- 3. classify ------------------------------------------------------------
echo "[3/3] HaploGrep3, tree=$TREE"
docker run --rm --user root -v "$WORK:/work" -v "$TREECACHE:$HG3/trees" \
  -w "$HG3" "$IMAGE" \
  ./haplogrep3 classify --tree "$TREE" \
    --in  "/work/results/haplogroup_deep/cohort.hsd" \
    --out "/work/results/haplogroup_deep/deep_${TAG}.txt" \
    --extend-report

echo
echo "wrote $OUT/deep_${TAG}.txt"
echo "Columns: terminal haplogroup, HaploGrep3 quality (fraction of expected"
echo "defining markers found; 0.90 is the conventional acceptance threshold),"
echo "markers found, markers absent, and variants not placed on the path."
