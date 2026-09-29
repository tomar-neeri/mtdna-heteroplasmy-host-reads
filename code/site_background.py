#!/usr/bin/env python3
"""
Build a site-specific background table (panel of normals) across the cohort.

WHY THIS IS NEEDED
A global or even per-substitution error rate cannot describe a position whose
local error rate is elevated. The background used to fit those models is defined
as VAF < 0.5%, so the fitted rate is capped at 0.005 by construction. A position
whose true local error is ~1.5% is therefore excluded from the fit as if it were
a variant, and then clears the threshold as a variant.

Observed directly in batch 1: m.6438T>A appeared in 8 of 11 samples at ~1.5%,
m.7549T>A in 6 of 11. If heteroplasmy were private, seeing one position in 8 of
11 samples has probability ~2e-28. These are error hotspots.

HOW IT WORKS
For every (position, alt allele) this records the summed alt and depth counts
across all samples. The caller then estimates the background for sample i by
LEAVING THAT SAMPLE OUT:

    site_rate_i = (sum_alt - alt_i) / (sum_depth - depth_i)

A private heteroplasmy leaves the other samples near zero, so its site rate stays
low and the call survives. A recurrent artifact gives every carrier a high site
rate, so all of them fail. No blacklist and no arbitrary recurrence cutoff --
the cohort estimates its own position-specific background.

Usage:
  site_background.py --bams *.bam --ref rCRS.fasta --out site_background.tsv
"""
import argparse, csv, sys
import numpy as np
import pysam

BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
IBASES = "ACGT"
L = 16569


def load_ref(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
    if len(seq) != L:
        sys.exit(f"reference length {len(seq)} != {L}")
    return seq


def pileup(bam_path, contig, min_bq, min_mq):
    counts = np.zeros((L, 4), dtype=np.int64)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for col in bam.pileup(contig, 0, L, truncate=True,
                              min_base_quality=min_bq, max_depth=1_000_000):
            p = col.reference_pos
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip:
                    continue
                if pr.alignment.mapping_quality < min_mq:
                    continue
                b = pr.alignment.query_sequence[pr.query_position]
                if b in BASES:
                    counts[p, BASES[b]] += 1
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bams", nargs="+", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--contig", default="chrM")
    ap.add_argument("--min-depth", type=int, default=50)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--detect-vaf", type=float, default=0.01,
                    help="VAF at which a sample counts as 'detecting' the allele")
    a = ap.parse_args()

    ref = load_ref(a.ref)
    n_samples = len(a.bams)
    sum_alt = np.zeros((L, 4), dtype=np.int64)
    sum_dep = np.zeros((L, 4), dtype=np.int64)
    n_detect = np.zeros((L, 4), dtype=np.int32)
    n_cov = np.zeros((L, 4), dtype=np.int32)

    for b in a.bams:
        c = pileup(b, a.contig, a.min_baseq, a.min_mapq)
        dep = c.sum(1)
        ok = dep >= a.min_depth
        for i in range(4):
            sum_alt[ok, i] += c[ok, i]
            sum_dep[ok, i] += dep[ok]
            n_cov[ok, i] += 1
            det = ok & (dep > 0) & (c[:, i] / np.maximum(dep, 1) >= a.detect_vaf)
            n_detect[det, i] += 1
        print(f"  {b}: median depth {np.median(dep[ok]) if ok.any() else 0:.0f}x",
              file=sys.stderr)

    rows = 0
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pos", "ref", "alt", "sum_alt", "sum_depth",
                    "n_samples_covered", "n_samples_detected", "pooled_vaf"])
        for p in range(L):
            if ref[p] not in BASES:
                continue
            ri = BASES[ref[p]]
            for i in range(4):
                if i == ri or sum_dep[p, i] == 0:
                    continue
                # only record positions carrying signal; the rest add nothing
                if sum_alt[p, i] == 0:
                    continue
                w.writerow([p + 1, ref[p], IBASES[i], int(sum_alt[p, i]),
                            int(sum_dep[p, i]), int(n_cov[p, i]),
                            int(n_detect[p, i]),
                            f"{sum_alt[p, i] / sum_dep[p, i]:.6g}"])
                rows += 1

    # report the worst hotspots
    flat = [(int(n_detect[p, i]), p + 1, ref[p], IBASES[i],
             sum_alt[p, i] / sum_dep[p, i])
            for p in range(L) if ref[p] in BASES
            for i in range(4)
            if i != BASES[ref[p]] and sum_dep[p, i] > 0 and n_detect[p, i] >= 2]
    flat.sort(reverse=True)
    print(f"\n  {rows:,} site/allele combinations recorded from {n_samples} samples",
          file=sys.stderr)
    if flat:
        print(f"  {len(flat)} alleles detected in >=2 samples "
              f"(candidate error hotspots):", file=sys.stderr)
        for n, p, r, alt, v in flat[:15]:
            print(f"    m.{p}{r}>{alt:<3} {n}/{n_samples} samples, "
                  f"pooled VAF {v:.3%}", file=sys.stderr)


if __name__ == "__main__":
    main()
