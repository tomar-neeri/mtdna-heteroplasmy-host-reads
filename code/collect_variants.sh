#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# collect_variants.sh — consolidate per-sample calls into cohort tables
#
#   bash collect_variants.sh                 # uses ./results
#   WORK=~/mtdna-batch1 bash collect_variants.sh
#
# Writes to results/summary/:
#   heteroplasmic_variants.tsv   every passing call, one row per sample/position
#   homoplasmic_variants.tsv     haplogroup backbone per sample
#   variant_recurrence.tsv       positions shared across samples
#   region_distribution.tsv      counts per mtDNA region
#   variants_annotated.tsv       gene, codon and amino-acid consequence
# ---------------------------------------------------------------------------
set -euo pipefail

WORK="${WORK:-$PWD}"
CALLS="$WORK/results/calls"
OUT="$WORK/results/summary"
mkdir -p "$OUT"

ls "$CALLS"/*.calls.tsv >/dev/null 2>&1 || { echo "No call files in $CALLS"; exit 1; }
hr(){ printf '%.0s-' {1..72}; echo; }

short(){ basename "$1" .calls.tsv | sed 's/24D214-5G_//; s/_S[0-9]*$//'; }

# --- 1. all passing heteroplasmic calls ------------------------------------
{
  printf 'sample\tpos\tref\talt\tdepth\talt_reads\tvaf\tfwd_alt\trev_alt\tp_error\tstrand_bias_p\tregion\n'
  for f in "$CALLS"/*.calls.tsv; do
    s=$(short "$f")
    awk -F'\t' -v s="$s" 'NR>1 && $12=="heteroplasmic" && $15=="True" {
      print s"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' "$f"
  done
} > "$OUT/heteroplasmic_variants.tsv"

# --- 2. homoplasmic (haplogroup backbone) ----------------------------------
{
  printf 'sample\tpos\tref\talt\tdepth\tvaf\tregion\n'
  for f in "$CALLS"/*.calls.tsv; do
    s=$(short "$f")
    awk -F'\t' -v s="$s" 'NR>1 && $12=="homoplasmic" {
      print s"\t"$1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$11}' "$f"
  done
} > "$OUT/homoplasmic_variants.tsv"

# --- 3. recurrence ---------------------------------------------------------
# Genuine heteroplasmy is mostly private. A position recurring across samples
# is either haplogroup-linked leakage or a systematic artifact.
{
  printf 'n_samples\tpos\tchange\tsamples\tmedian_vaf\n'
  awk -F'\t' 'NR>1 {k=$2"\t"$3">"$4; n[k]++; s[k]=s[k]","$1; v[k]=v[k]" "$7}
    END{for(k in n){gsub(/^,/,"",s[k]); print n[k]"\t"k"\t"s[k]"\t"v[k]}}' \
    "$OUT/heteroplasmic_variants.tsv" \
  | sort -k1,1nr -k2,2n
} > "$OUT/variant_recurrence.tsv"

# --- 4. region distribution ------------------------------------------------
{
  printf 'region\tn_variants\tn_samples\tmedian_vaf\n'
  awk -F'\t' 'NR>1{n[$12]++; seen[$12"|"$1]=1; v[$12]=v[$12]" "$7}
    END{for(r in n){c=0; for(k in seen) if(index(k,r"|")==1) c++;
      print r"\t"n[r]"\t"c"\t"v[r]}}' "$OUT/heteroplasmic_variants.tsv" \
  | sort -k2,2nr
} > "$OUT/region_distribution.tsv"

# --- 5. annotate protein-coding consequences -------------------------------
python3 - "$OUT/heteroplasmic_variants.tsv" "$OUT/variants_annotated.tsv" <<'PYEOF'
import sys, csv
VMC = {}
bases = "TCAG"
aas   = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
i = 0
for b1 in bases:
    for b2 in bases:
        for b3 in bases:
            VMC[b1+b2+b3] = aas[i]; i += 1
# vertebrate mitochondrial deviations from the standard code
VMC.update({"AGA":"*","AGG":"*","ATA":"M","TGA":"W"})

GENES = [("MT-ND1",3307,4262),("MT-ND2",4470,5511),("MT-CO1",5904,7445),
         ("MT-CO2",7586,8269),("MT-ATP8",8366,8572),("MT-ATP6",8527,9207),
         ("MT-CO3",9207,9990),("MT-ND3",10059,10404),("MT-ND4L",10470,10766),
         ("MT-ND4",10760,12137),("MT-ND5",12337,14148),("MT-ND6",14149,14673),
         ("MT-CYB",14747,15887)]

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
with open(sys.argv[2], "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sample","variant","vaf","depth","region","gene","codon",
                "aa_change","effect"])
    for r in rows:
        pos = int(r["pos"]); ref, alt = r["ref"], r["alt"]
        gene = codon = aa = eff = "-"
        for g, s, e in GENES:
            if s <= pos <= e:
                gene = g
                off = pos - s; ci = off // 3; cp = off % 3
                eff = "coding"   # codon context needs the reference sequence
                codon = str(ci + 1)
                break
        else:
            if pos >= 16024 or pos <= 576:
                gene = "control_region"; eff = "non-coding"
            else:
                gene = "tRNA/rRNA"; eff = "non-coding"
        w.writerow([r["sample"], f"m.{pos}{ref}>{alt}",
                    f"{float(r['vaf'])*100:.2f}%", r["depth"], r["region"],
                    gene, codon, aa, eff])
PYEOF

# --- report ----------------------------------------------------------------
hr; echo "  Variant tables written to $OUT"; hr
n=$(( $(wc -l < "$OUT/heteroplasmic_variants.tsv") - 1 ))
m=$(( $(wc -l < "$OUT/homoplasmic_variants.tsv") - 1 ))
echo "  heteroplasmic calls : $n"
echo "  homoplasmic calls   : $m"
echo
echo "  Recurrence (positions seen in >1 sample):"
awk -F'\t' 'NR>1 && $1>1{printf "    %s samples  m.%s%s   [%s]\n",$1,$2,$3,$4}' \
    "$OUT/variant_recurrence.tsv" | head -15
rec=$(awk -F'\t' 'NR>1 && $1>1' "$OUT/variant_recurrence.tsv" | wc -l)
[ "$rec" -eq 0 ] && echo "    none — all calls are private to individuals (expected)"
echo
echo "  Region distribution:"
awk -F'\t' 'NR>1{printf "    %-18s %3s variants in %2s samples\n",$1,$2,$3}' \
    "$OUT/region_distribution.tsv"
hr
echo "  Files:"
for f in heteroplasmic_variants homoplasmic_variants variant_recurrence \
         region_distribution variants_annotated; do
  printf "    %-32s %s\n" "$f.tsv" "$(wc -l < "$OUT/$f.tsv") lines"
done
hr
