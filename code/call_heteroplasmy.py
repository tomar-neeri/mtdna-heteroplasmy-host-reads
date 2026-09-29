#!/usr/bin/env python3
"""
Call heteroplasmy against a POOLED error model, with a hard alt-read floor.

Two filters do the work:
  - significance against the pooled beta-binomial null
  - min_alt_reads, applied regardless of what the model says

The floor is not redundant. It is what separated real signal in a high-host library
(4 candidates surviving >=5 alt reads) from noise in a low-host library (0 surviving,
though a per-sample model had reported 22).

Also reports callable fraction, because mean depth misleads: a low-host library
averaged 175x yet only 33.6% of its genome was callable at 3% VAF.

Usage:
  call_heteroplasmy.py --bam s.bam --ref rCRS.fasta --model model.json \
                       --out s.calls.tsv --summary s.summary.json
"""
import argparse, csv, json, os, sys
import numpy as np
import pysam
from scipy.stats import betabinom, binom, fisher_exact

BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
IBASES = "ACGT"
L = 16569

REGIONS = [
    ("control_region", lambda p: p >= 16024 or p <= 576),
    ("MT-RNR1",        lambda p: 648 <= p <= 1601),
    ("MT-RNR2",        lambda p: 1671 <= p <= 3229),
    ("protein_coding", lambda p: 3307 <= p <= 15887),
]


def region_of(pos):
    for name, test in REGIONS:
        if test(pos):
            return name
    return "tRNA_other"


def load_ref(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">"))
    if len(seq) != L:
        sys.exit(f"reference length {len(seq)} != {L}")
    return seq.upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--contig", default="chrM")
    ap.add_argument("--min-depth", type=int, default=50)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--min-alt-reads", type=int, default=4)
    ap.add_argument("--vaf-floor", type=float, default=0.01)
    ap.add_argument("--vaf-homoplasmy", type=float, default=0.95)
    ap.add_argument("--sig-level", type=float, default=0.05)
    ap.add_argument("--strand-bias-p", type=float, default=0.01)
    ap.add_argument("--site-background", default=None,
                    help="panel-of-normals table from site_background.py. "
                         "Enables leave-one-out site-specific thresholds, which "
                         "is what controls recurrent error hotspots.")
    ap.add_argument("--min-loo-samples", type=int, default=3,
                    help="minimum other samples needed for a site estimate")
    ap.add_argument("--mask", default="302-315,16184-16193",
                    help="1-based inclusive ranges to exclude (poly-C tracts)")
    a = ap.parse_args()

    ref = load_ref(a.ref)
    model = json.load(open(a.model))
    subs = model.get("by_substitution", {})
    gmod = model.get("global", model)

    def params_for(sub):
        """Per-substitution parameters, falling back to the global model.

        Necessary because error rates differ by orders of magnitude between
        substitution types. A global threshold admits T>A and A>C artifacts
        while burying real C>T and T>C variation.
        """
        m = subs.get(sub, gmod)
        al, be = m["alpha"], m["beta"]
        return al, be, (al + be) > 1e6, al / (al + be)

    def tail_prob(k, n, sub):
        al, be, use_binom, rate = params_for(sub)
        if use_binom:
            return float(binom.sf(k - 1, int(n), rate))
        return float(betabinom.sf(k - 1, int(n), al, be))

    def site_rate_loo(pos, ref_b, alt_b, my_alt, my_depth):
        """Background at this exact site, estimated from the OTHER samples.

        Leaving this sample out is what makes the estimate usable: a private
        heteroplasmy contributes nothing to the others, so its site rate stays
        low; a recurrent artifact raises the rate for every carrier.
        """
        e = pon.get((pos, ref_b, alt_b))
        if not e:
            return None, 0
        s_alt, s_dep, n_cov, n_det = e
        o_alt, o_dep = s_alt - my_alt, s_dep - my_depth
        n_other = n_cov - 1
        if o_dep <= 0 or n_other < a.min_loo_samples:
            return None, n_det
        return o_alt / o_dep, n_det

    # --- panel of normals -------------------------------------------------
    pon = {}
    if a.site_background and os.path.exists(a.site_background):
        with open(a.site_background, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                pon[(int(r["pos"]), r["ref"], r["alt"])] = (
                    int(r["sum_alt"]), int(r["sum_depth"]),
                    int(r["n_samples_covered"]), int(r["n_samples_detected"]))
        print(f"  panel of normals: {len(pon):,} site/allele entries",
              file=sys.stderr)

    masked = set()
    for rng in filter(None, a.mask.split(",")):
        lo, hi = (int(x) for x in rng.split("-"))
        masked.update(range(lo, hi + 1))

    counts = np.zeros((L, 4), np.int32)
    fwd = np.zeros((L, 4), np.int32)
    rev = np.zeros((L, 4), np.int32)
    with pysam.AlignmentFile(a.bam, "rb") as bam:
        for col in bam.pileup(a.contig, 0, L, truncate=True,
                              min_base_quality=a.min_baseq, max_depth=1_000_000):
            p = col.reference_pos
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip:
                    continue
                if pr.alignment.mapping_quality < a.min_mapq:
                    continue
                b = pr.alignment.query_sequence[pr.query_position]
                if b in BASES:
                    i = BASES[b]
                    counts[p, i] += 1
                    (rev if pr.alignment.is_reverse else fwd)[p, i] += 1

    dep = counts.sum(1)
    rows = []
    for p in range(L):
        pos1 = p + 1
        if dep[p] < a.min_depth or ref[p] not in BASES or pos1 in masked:
            continue
        ri = BASES[ref[p]]
        ac, ai = max((counts[p, i], i) for i in range(4) if i != ri)
        if ac < a.min_alt_reads:
            continue
        vaf = ac / dep[p]
        if vaf < a.vaf_floor:
            continue
        sub = f"{ref[p]}>{IBASES[ai]}"
        p_sub = tail_prob(ac, dep[p], sub)

        # site-specific null from the other samples; use whichever background
        # is higher, so a hotspot cannot be rescued by a low global rate
        sr, n_det = site_rate_loo(pos1, ref[p], IBASES[ai], int(ac), int(dep[p]))
        _, _, _, sub_rate = params_for(sub)
        if sr is not None and sr > sub_rate:
            p_site = float(binom.sf(ac - 1, int(dep[p]), sr))
            p_err = max(p_sub, p_site)
            null_src = "site"
        else:
            p_err = p_sub
            null_src = "substitution"
        try:
            sb = float(fisher_exact([[fwd[p, ai], fwd[p, ri]],
                                     [rev[p, ai], rev[p, ri]]])[1])
        except Exception:
            sb = float("nan")
        rows.append({
            "pos": pos1, "ref": ref[p], "alt": IBASES[ai],
            "depth": int(dep[p]), "alt_reads": int(ac), "vaf": round(vaf, 5),
            "fwd_alt": int(fwd[p, ai]), "rev_alt": int(rev[p, ai]),
            "p_error": p_err, "strand_bias_p": sb,
            "substitution": sub,
            "sub_error_rate": params_for(sub)[3],
            "site_bg_rate": (f"{sr:.3e}" if sr is not None else ""),
            "n_other_samples_detected": n_det,
            "null_source": null_src,
            "sub_model": "fitted" if sub in subs and
                         subs[sub].get("source") == "fitted" else "global",
            "region": region_of(pos1),
            "class": "homoplasmic" if vaf >= a.vaf_homoplasmy else "heteroplasmic",
            "pass_significance": p_err < a.sig_level,
            "pass_strand_bias": not (sb == sb) or sb > a.strand_bias_p,
        })

    for r in rows:
        r["PASS"] = bool(r["pass_significance"] and r["pass_strand_bias"])

    cols = ["pos", "ref", "alt", "depth", "alt_reads", "vaf", "fwd_alt", "rev_alt",
            "p_error", "strand_bias_p", "region", "class",
            "pass_significance", "pass_strand_bias", "PASS",
            "substitution", "sub_error_rate", "sub_model",
            "site_bg_rate", "n_other_samples_detected", "null_source"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # callable fraction: what VAF could actually be detected, position by position
    callable_frac = {}
    for target in (0.02, 0.03, 0.05, 0.10):
        need = a.min_alt_reads / target
        callable_frac[f"{target:.0%}"] = float((dep >= need).mean())

    het = [r for r in rows if r["class"] == "heteroplasmic"]
    from collections import Counter
    sub_counts = Counter(r["substitution"] for r in het if r["PASS"])
    summary = {
        "sample": a.bam.split("/")[-1].split(".")[0],
        "mean_depth": float(dep.mean()),
        "median_depth": float(np.median(dep)),
        "zero_coverage_bases": int((dep == 0).sum()),
        "error_model_phred": gmod.get("phred_equivalent"),
        "error_model_form": "per-substitution beta-binomial" if subs
                            else "global beta-binomial",
        "n_substitution_models": len([k for k, v in subs.items()
                                      if v.get("source") == "fitted"]),
        "panel_of_normals": bool(pon),
        "n_rejected_by_site_background": sum(
            1 for r in rows if r.get("null_source") == "site"
            and not r["pass_significance"]),
        "min_alt_reads": a.min_alt_reads,
        "callable_fraction_at_vaf": callable_frac,
        "n_homoplasmic": sum(1 for r in rows if r["class"] == "homoplasmic"),
        "n_heteroplasmic_raw": len(het),
        "n_heteroplasmic_pass": sum(1 for r in het if r["PASS"]),
        "heteroplasmic_by_substitution": dict(sub_counts.most_common()),
        "heteroplasmic_by_region": {
            k: sum(1 for r in het if r["PASS"] and r["region"] == k)
            for k, _ in REGIONS
        },
    }
    with open(a.summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"{summary['sample']}: depth {summary['mean_depth']:.0f}x | "
          f"hom {summary['n_homoplasmic']} | "
          f"het {summary['n_heteroplasmic_pass']}/{summary['n_heteroplasmic_raw']} pass | "
          f"callable@3% {callable_frac['3%']:.1%}", file=sys.stderr)


if __name__ == "__main__":
    main()
