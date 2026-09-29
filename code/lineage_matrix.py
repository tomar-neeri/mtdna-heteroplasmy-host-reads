#!/usr/bin/env python3
"""
Cross-tabulate heteroplasmic calls against the cohort's homoplasmic landscape.

THE QUESTION
A low-VAF call at a position that is homoplasmic (i.e. lineage-defining) in
other samples is ambiguous: it could be private heteroplasmy, or leakage of a
lineage marker. A low-VAF call at a position homoplasmic in nobody is
unambiguously private.

This builds the mapping directly: for every heteroplasmic call, how many cohort
samples carry that position homoplasmically, and are those samples in the same
macrohaplogroup as the caller?

INTERPRETATION
  homoplasmic in 0 samples        -> private; lineage plays no part
  homoplasmic within own group    -> possible lineage leakage; inspect
  homoplasmic in other group only -> more likely a mutational hotspot shared
                                     across lineages than leakage

Usage:
  lineage_matrix.py --calls results/calls --haplogroups results/summary/haplogroups.tsv \
                    --out results/summary/lineage_matrix.tsv
"""
import argparse, csv, glob, os, re, sys
from collections import defaultdict, Counter


def short(p):
    n = os.path.basename(p).split(".calls.tsv")[0]
    return re.sub(r"_S\d+$", "", n.replace("24D214-5G_", ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--haplogroups", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--vaf-homoplasmy", type=float, default=0.95)
    a = ap.parse_args()

    hg = {}
    for r in csv.DictReader(open(a.haplogroups), delimiter="\t"):
        hg[r["sample"]] = r["macrohaplogroup"]

    hom = defaultdict(set)      # pos -> {sample}
    het = []                    # (sample, pos, ref, alt, vaf, depth)
    for f in sorted(glob.glob(os.path.join(a.calls, "*.calls.tsv"))):
        s = short(f)
        for r in csv.DictReader(open(f), delimiter="\t"):
            if r.get("PASS") != "True":
                continue
            try:
                vaf = float(r["vaf"]); pos = int(r["pos"])
            except (KeyError, ValueError):
                continue
            if vaf >= a.vaf_homoplasmy:
                hom[pos].add(s)
            else:
                het.append((s, pos, r["ref"], r["alt"], vaf, int(r["depth"])))

    groups = Counter(hg.values())
    rows, summary = [], Counter()
    for s, pos, ref, alt, vaf, dep in sorted(het, key=lambda x: -x[4]):
        mine = hg.get(s, "?")
        carriers = hom.get(pos, set())
        same = {c for c in carriers if hg.get(c) == mine}
        other = carriers - same
        if not carriers:
            verdict = "PRIVATE"
        elif same and not other:
            verdict = "lineage_own_group"
        elif other and not same:
            verdict = "hotspot_other_group"
        else:
            verdict = "hotspot_both_groups"
        summary[verdict] += 1
        rows.append([s, mine, f"m.{pos}{ref}>{alt}", f"{vaf*100:.2f}", dep,
                     len(carriers), len(same), len(other),
                     ",".join(sorted(hg.get(c, "?") for c in carriers)) or "-",
                     verdict])

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "haplogroup", "variant", "vaf_pct", "depth",
                    "n_homoplasmic_carriers", "n_same_group", "n_other_group",
                    "carrier_haplogroups", "verdict"])
        w.writerows(rows)

    print(f"  cohort: " + ", ".join(f"{k} {v}" for k, v in groups.most_common()),
          file=sys.stderr)
    print(f"  heteroplasmic calls: {len(rows)}\n", file=sys.stderr)
    print(f"  {'verdict':<22}{'n':>5}{'share':>8}", file=sys.stderr)
    print("  " + "-" * 34, file=sys.stderr)
    for k, v in summary.most_common():
        print(f"  {k:<22}{v:>5}{v/len(rows):>8.1%}", file=sys.stderr)
    amb = sum(v for k, v in summary.items() if k != "PRIVATE")
    print(f"\n  unambiguously private: {summary['PRIVATE']}/{len(rows)} "
          f"({summary['PRIVATE']/max(len(rows),1):.1%})", file=sys.stderr)
    if amb:
        print(f"  requiring inspection : {amb}", file=sys.stderr)
    print(f"\n  written: {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
