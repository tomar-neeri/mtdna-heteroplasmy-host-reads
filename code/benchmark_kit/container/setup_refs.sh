#!/usr/bin/env bash
# Runs INSIDE the container. Builds the reference files the two external callers
# need, without writing to the read-only /refs mount.
#
#   /bench/ref/hg38.fa(+.fai,.dict)  for Mutect2 on the final BAMs (hg38 header)
#   /bench/ref/chrM.fa(+.fai,.dict)  for mutserve, and for Mutect2 on the
#                                    unfiltered Pass-1 BAMs (rCRS-only header)
set -euo pipefail
R=/bench/ref; mkdir -p "$R"

# hg38: symlink the FASTA, build index and dictionary beside the symlink
[ -e "$R/hg38.fa" ] || ln -s /refs/hg38.fa "$R/hg38.fa"
if [ ! -s "$R/hg38.fa.fai" ]; then
  if [ -s /refs/hg38.fa.fai ]; then cp /refs/hg38.fa.fai "$R/hg38.fa.fai"
  else samtools faidx "$R/hg38.fa"; fi
fi
[ -s "$R/hg38.dict" ] || gatk --java-options "-Xmx4g" CreateSequenceDictionary -R "$R/hg38.fa" -O "$R/hg38.dict"

# chrM as its own reference. hg38 chrM is the rCRS (NC_012920.1); check it.
if [ ! -s "$R/chrM.fa" ]; then
  samtools faidx "$R/hg38.fa" chrM > "$R/chrM.fa"
fi
python3 - "$R/chrM.fa" "${RCRS:-/opt/refs/rCRS.fasta}" <<'EOF'
import sys
seq = lambda p: "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()
a, b = seq(sys.argv[1]), seq(sys.argv[2])
assert len(a) == 16569, f"hg38 chrM length {len(a)}"
diff = sum(x != y for x, y in zip(a, b))
print(f"  hg38 chrM vs rCRS: {len(a)} bp, {diff} mismatches (N placeholder at 3107 expected)")
assert diff <= 1, "hg38 chrM is not the rCRS"
EOF
[ -s "$R/chrM.fa.fai" ] || samtools faidx "$R/chrM.fa"
[ -s "$R/chrM.dict" ]  || gatk CreateSequenceDictionary -R "$R/chrM.fa" -O "$R/chrM.dict"

echo "  references ready in $R"
gatk --version 2>&1 | grep -i "gatk\|htsjdk" | head -2
java -jar "${MUTSERVE_JAR:-/opt/tools/mutserve/mutserve.jar}" 2>&1 | grep -i -m1 "mutserve" || true
