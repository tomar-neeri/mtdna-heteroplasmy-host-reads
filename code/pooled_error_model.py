#!/usr/bin/env python3
"""
Fit beta-binomial sequencing-error models POOLED across samples, PER SUBSTITUTION.

A single global error rate is not sufficient. Validated on batch 1 (11 samples,
66 calls): T>A accounted for 47% of calls against 2.27% of gnomAD chrM variants
-- 21-fold enriched, p = 3e-33. A>C was 3.3-fold enriched. Meanwhile C>T, which
is 21% of gnomAD variants, produced zero calls. A global threshold calibrated on
the average admits substitution types whose true error rate is far above it, and
buries the types whose rate is below it.

Strand bias cannot compensate: T>A errors proved strand-SYMMETRIC (86% of raw
candidates passed the strand test) while A>C and T>G were strand-asymmetric and
correctly rejected (23% passed). The filter's power is uneven by substitution,
so the threshold itself must be substitution-specific.

Why pooled: the error rate is a property of the sequencing chemistry, not of
an individual sample. Fitting per-sample fails predictably below ~200x median
depth, where most background positions carry zero alt reads and the MLE drives
alpha toward zero. Validation case: a low-host library (117x median) fitted 1.19e-06
(Q59, impossible for Illumina) and called 22 false heteroplasmies; pooling with
a high-host library gave 2.69e-04 (Q36) and left zero.

Usage:
  pooled_error_model.py --bams a.bam b.bam --ref rCRS.fasta --out model.json
"""
import argparse, json, sys
import numpy as np
import pysam
from scipy.stats import betabinom
from scipy.optimize import minimize

BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
L = 16569


def load_ref(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">"))
    if len(seq) != L:
        sys.exit(f"reference length {len(seq)} != {L}")
    return seq.upper()


def pileup_counts(bam_path, contig, min_bq, min_mq):
    counts = np.zeros((L, 4), dtype=np.int32)
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


IBASES = "ACGT"


def background_by_substitution(counts, ref, min_depth, bg_vaf_max):
    """Collect background observations keyed by substitution type.

    Every non-reference allele is counted separately, not just the top one.
    That gives three observations per position instead of one, and avoids the
    bias of always selecting whichever substitution happens to have the highest
    error rate at that site.
    """
    dep = counts.sum(1)
    out = {}
    for p in range(L):
        if dep[p] < min_depth or ref[p] not in BASES:
            continue
        ri = BASES[ref[p]]
        for i in range(4):
            if i == ri:
                continue
            ac = counts[p, i]
            if ac / dep[p] < bg_vaf_max:
                out.setdefault(f"{ref[p]}>{IBASES[i]}", [[], []])
                out[f"{ref[p]}>{IBASES[i]}"][0].append(ac)
                out[f"{ref[p]}>{IBASES[i]}"][1].append(dep[p])
    return {k: (np.asarray(v[0], float), np.asarray(v[1], float))
            for k, v in out.items()}


def fit(alt, dep):
    def nll(q):
        a, b = np.exp(q)
        if not (np.isfinite(a) and np.isfinite(b)) or a <= 0 or b <= 0:
            return 1e12
        v = betabinom.logpmf(alt, dep, a, b)
        return 1e12 if not np.all(np.isfinite(v)) else -v.sum()

    # Bound the concentration. Unbounded, the optimizer drives alpha+beta toward
    # 1e17 to drive overdispersion to zero; betabinom then overflows and sf()
    # returns exactly 0 for every position, silently disabling the significance
    # filter. MAX_CONC keeps the distribution evaluable.
    MAX_CONC = 1e6

    def nll_bounded(q):
        a, b = np.exp(q)
        if a + b > MAX_CONC:
            return 1e12
        return nll(q)

    best = None
    for a0, b0 in [(0.1, 100.0), (1.0, 1000.0), (0.5, 500.0), (0.05, 50.0)]:
        r = minimize(nll_bounded, [np.log(a0), np.log(b0)], method="Nelder-Mead",
                     options={"maxiter": 5000, "fatol": 1e-9, "xatol": 1e-9})
        if best is None or r.fun < best.fun:
            best = r
    return float(np.exp(best.x[0])), float(np.exp(best.x[1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bams", nargs="+", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--contig", default="chrM")
    ap.add_argument("--min-depth", type=int, default=50)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--bg-vaf-max", type=float, default=0.005)
    ap.add_argument("--min-background", type=int, default=500,
                    help="minimum background observations to fit a substitution "
                         "its own model; below this it inherits the global fit")
    a = ap.parse_args()

    ref = load_ref(a.ref)
    per_sub, per_sample = {}, {}
    all_alt, all_dep = [], []

    for b in a.bams:
        counts = pileup_counts(b, a.contig, a.min_baseq, a.min_mapq)
        bg = background_by_substitution(counts, ref, a.min_depth, a.bg_vaf_max)
        if not bg:
            print(f"  {b}: no usable background, skipped", file=sys.stderr)
            continue
        npos = 0
        for k, (alt, dep) in bg.items():
            per_sub.setdefault(k, [[], []])
            per_sub[k][0].append(alt)
            per_sub[k][1].append(dep)
            all_alt.append(alt)
            all_dep.append(dep)
            npos += len(alt)
        dep_all = np.concatenate([v[1] for v in bg.values()])
        per_sample[b] = {"background_observations": int(npos),
                         "median_depth": float(np.median(dep_all))}
        print(f"  {b}: {npos:,} background observations, "
              f"median depth {np.median(dep_all):.0f}x", file=sys.stderr)

    if not all_alt:
        sys.exit("no background data from any sample")

    # --- global model (fallback for sparse substitutions) ---
    g_alt = np.concatenate(all_alt)
    g_dep = np.concatenate(all_dep)
    ga, gb = fit(g_alt, g_dep)
    g_err = ga / (ga + gb)
    print(f"\nGLOBAL : error {g_err:.3e} (~Q{-10 * np.log10(g_err):.1f}) "
          f"over {len(g_alt):,} observations", file=sys.stderr)

    # --- per-substitution models ---
    subs = {}
    print(f"\n{'sub':<6}{'n_bg':>10}{'error':>12}{'Q':>7}{'rel.':>8}  source",
          file=sys.stderr)
    print("-" * 52, file=sys.stderr)
    for k in sorted(per_sub, key=lambda x: -sum(len(v) for v in per_sub[x][0])):
        alt = np.concatenate(per_sub[k][0])
        dep = np.concatenate(per_sub[k][1])
        if len(alt) >= a.min_background:
            sa, sb = fit(alt, dep)
            src = "fitted"
        else:
            sa, sb, src = ga, gb, "global (sparse)"
        e = sa / (sa + sb)
        subs[k] = {"alpha": float(sa), "beta": float(sb), "mean_error": float(e),
                   "phred_equivalent": float(-10 * np.log10(e)) if e > 0 else None,
                   "n_background": int(len(alt)), "source": src}
        print(f"{k:<6}{len(alt):>10,}{e:>12.3e}{-10*np.log10(e) if e>0 else 0:>7.1f}"
              f"{e/g_err:>8.1f}x  {src}", file=sys.stderr)

    hi = sorted(((v["mean_error"] / g_err, k) for k, v in subs.items()
                 if v["source"] == "fitted"), reverse=True)[:3]
    if hi and hi[0][0] > 2:
        print(f"\n  Substitutions well above the global rate: " +
              ", ".join(f"{k} ({r:.1f}x)" for r, k in hi if r > 2), file=sys.stderr)
        print("  A global threshold would admit these as false heteroplasmy.",
              file=sys.stderr)

    model = {
        "global": {"alpha": float(ga), "beta": float(gb), "mean_error": float(g_err),
                   "phred_equivalent": float(-10 * np.log10(g_err)) if g_err > 0 else None,
                   "overdispersion": float(1.0 / (ga + gb + 1.0)),
                   "use_binomial": bool(1.0 / (ga + gb + 1.0) < 1e-6),
                   "n_background": int(len(g_alt))},
        "by_substitution": subs,
        "n_samples": len(per_sample),
        "per_sample": per_sample,
        "params": {k: getattr(a, k) for k in
                   ("min_depth", "min_baseq", "min_mapq", "bg_vaf_max",
                    "min_background")},
    }
    # keep the old top-level keys so existing consumers still work
    model.update({k: model["global"][k] for k in
                  ("alpha", "beta", "mean_error", "phred_equivalent",
                   "overdispersion", "use_binomial")})

    q = model["phred_equivalent"]
    if q and (q > 45 or q < 20):
        print(f"\n  WARNING: global Q{q:.0f} outside the plausible Illumina "
              f"range (Q20-Q45).", file=sys.stderr)

    with open(a.out, "w") as fh:
        json.dump(model, fh, indent=2)


if __name__ == "__main__":
    main()
