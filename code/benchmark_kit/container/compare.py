#!/usr/bin/env python3
"""
Benchmark the in-house heteroplasmy caller against Mutect2 (mitochondria mode)
and mutserve. Needs only the Python standard library (+ matplotlib for figures).

Every caller's output is put through the SAME harmonised definition of a
heteroplasmic call, so differences reflect the callers, not the thresholds:

    SNV, caller's own FILTER = PASS, depth >= 50, 1% <= VAF < 95%,
    outside the masked poly-C tracts (m.302-315, m.16184-16193)

Three analyses, each run only if its inputs are given:

  1. cohort      the 47 final (post-NUMT-removal) BAMs. Call counts, Ti/Tv,
                 T>A burden and ceiling, threshold sweep, overlap with the
                 in-house 134 calls (and specifically the 57 T>A calls).
  2. numt        the unfiltered Pass-1 BAMs vs the final BAMs, per caller:
                 how many calls each external caller adds when NUMT reads are
                 left in, and where those calls fall.
  3. mixtures    in-silico M x R mixtures with known truth: recall by allele
                 fraction and depth, false positives, allele-fraction bias.

  compare.py --ours /summary/heteroplasmic_variants.tsv --cohort /bench/cohort \
             [--unfiltered /bench/unfiltered] \
             [--mix /bench/mixtures --mix-tools /bench/mixtures/tools \
              --mix-ours /bench/mixtures/ours_pon --mix-ours-nopon /bench/mixtures/ours_nopon] \
             --out /bench/report
"""
import argparse, csv, glob, gzip, json, os, re, statistics, sys
from collections import Counter, defaultdict

MASK = set(range(302, 316)) | set(range(16184, 16194))
TI = {"A>G", "G>A", "C>T", "T>C"}
MIN_DEPTH, VAF_LO, VAF_HI = 50, 0.01, 0.95
SWEEP = [1, 2, 3, 5, 7, 10]
CALLERS = ["ours", "mutect2", "mutserve"]
LABEL = {"ours": "In-house (published)", "ours_nopon": "In-house, no panel of normals",
         "mutect2": "Mutect2 (mito mode)", "mutserve": "mutserve"}


def sid(x):
    """'24D214-5G_387_S45' / '387' / 'mix_425in418_d200_f05_r1' -> canonical id."""
    x = os.path.basename(x)
    for suf in (".filtered.vcf.gz", ".vcf.gz", ".txt", ".calls.tsv", ".chrM.bam", ".cand.bam", ".bam"):
        if x.endswith(suf):
            x = x[: -len(suf)]
    if x.startswith("mix_"):
        return x
    m = re.search(r"(?:^|_)(\d{3})(?:_S\d+)?$", x)
    return m.group(1) if m else x


# ------------------------------------------------------------------ parsers
# each returns {sample: [dict(pos, ref, alt, af, dp, filt)]}

def parse_mutect2(d):
    out = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(d, "*.vcf.gz"))):
        s = sid(f)
        out[s]  # register sample even when it has no records
        with gzip.open(f, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                pos, ref, alts, filt = int(c[1]), c[3], c[4].split(","), c[6]
                fmt = dict(zip(c[8].split(":"), c[9].split(":")))
                afs = fmt.get("AF", "").split(",")
                dp = int(fmt["DP"]) if fmt.get("DP", ".") not in (".", "") else None
                ad = fmt.get("AD", "").split(",")
                if dp is None and len(ad) > 1:
                    dp = sum(int(x) for x in ad if x not in (".", ""))
                for i, alt in enumerate(alts):
                    if len(ref) != 1 or len(alt) != 1 or i >= len(afs) or afs[i] in (".", ""):
                        continue
                    out[s].append(dict(pos=pos, ref=ref, alt=alt, af=float(afs[i]), dp=dp or 0, filt=filt))
    return out


def parse_mutserve(d):
    out = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(d, "*.txt"))):
        s = sid(f)
        out[s]
        with open(f) as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        if not rows:
            continue
        h = {k.strip().lower(): i for i, k in enumerate(rows[0])}
        def col(*names):
            for n in names:
                if n in h:
                    return h[n]
            raise KeyError(f"{f}: none of {names} in header {rows[0]}")
        ip, ir, iv, il = col("pos", "position"), col("ref", "reference"), col("variant"), col("variantlevel", "variant_level")
        ic = col("coverage", "cov")
        iflt = h.get("filter")
        for r in rows[1:]:
            if len(r) <= max(ip, ir, iv, il, ic):
                continue
            ref, alt = r[ir].strip().upper(), r[iv].strip().upper()
            if len(ref) != 1 or len(alt) != 1 or alt not in "ACGT" or alt == ref:
                continue
            out[s].append(dict(pos=int(r[ip]), ref=ref, alt=alt, af=float(r[il]),
                               dp=int(float(r[ic])), filt=(r[iflt].strip() if iflt is not None else "PASS")))
    return out


def parse_ours_table(path):
    """Cohort table heteroplasmic_variants.tsv/.csv: every row already passed."""
    out = defaultdict(list)
    with open(path, newline="") as fh:
        first = fh.readline(); fh.seek(0)
        for r in csv.DictReader(fh, delimiter="\t" if "\t" in first else ","):
            out[sid(r["sample"])].append(dict(pos=int(r["pos"]), ref=r["ref"], alt=r["alt"],
                                              af=float(r["vaf"]), dp=int(r["depth"]), filt="PASS"))
    return out


def parse_ours_dir(d):
    """Per-sample *.calls.tsv from call_heteroplasmy.py (all candidates + PASS flag)."""
    out = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(d, "*.calls.tsv"))):
        s = sid(f)
        out[s]
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out[s].append(dict(pos=int(r["pos"]), ref=r["ref"], alt=r["alt"], af=float(r["vaf"]),
                                   dp=int(r["depth"]), filt="PASS" if r["PASS"] == "True" else "FAIL"))
    return out


# ------------------------------------------------------------------ harmonise

def harmonise(calls):
    """Apply the shared definition of a heteroplasmic call."""
    out = {}
    for s, rs in calls.items():
        best = {}
        for r in rs:
            if r["filt"] != "PASS" or r["pos"] in MASK or r["dp"] < MIN_DEPTH:
                continue
            if not (VAF_LO <= r["af"] < VAF_HI):
                continue
            k = (r["pos"], r["alt"])
            if k not in best or r["af"] > best[k]["af"]:
                best[k] = r
        out[s] = best
    return out


def any_record(calls):
    """Every PASS record at any VAF < 95% (for 'did the caller see it at all')."""
    out = defaultdict(dict)
    for s, rs in calls.items():
        for r in rs:
            if r["filt"] == "PASS" and r["af"] < VAF_HI:
                out[s][(r["pos"], r["alt"])] = r
    return out


def sub(r):
    return f"{r['ref']}>{r['alt']}"


def spectrum(recs):
    n = len(recs)
    ti = sum(sub(r) in TI for r in recs)
    ta = [r for r in recs if sub(r) == "T>A"]
    return dict(n=n, ti=ti, tv=n - ti, titv=(ti / (n - ti)) if n - ti else float("nan"),
                n_TA=len(ta), pct_TA=100 * len(ta) / n if n else 0.0,
                n_CT=sum(sub(r) == "C>T" for r in recs),
                max_TA=100 * max((r["af"] for r in ta), default=0),
                max_af=100 * max((r["af"] for r in recs), default=0),
                ge5=sum(r["af"] >= 0.05 for r in recs), ge10=sum(r["af"] >= 0.10 for r in recs))


def region(pos):
    if pos >= 16024 or pos <= 576: return "control_region"
    if 648 <= pos <= 1601: return "MT-RNR1"
    if 1671 <= pos <= 3229: return "MT-RNR2"
    return "other"


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(header); w.writerows(rows)


def fmt(x, d=2):
    return "—" if isinstance(x, float) and x != x else (f"{x:.{d}f}" if isinstance(x, float) else str(x))


# ------------------------------------------------------------------ analyses

def cohort_analysis(ours_raw, tools_raw, out, md):
    H = {k: harmonise(v) for k, v in {"ours": ours_raw, **tools_raw}.items()}
    A = {k: any_record(v) for k, v in tools_raw.items()}
    samples = sorted(set(ours_raw) | {s for v in tools_raw.values() for s in v}, key=lambda x: (len(x), x))

    # 1. summary per caller
    rows = []
    for c in H:
        recs = [r for s in H[c] for r in H[c][s].values()]
        sp = spectrum(recs)
        rows.append([LABEL[c], sp["n"], sum(1 for s in H[c] if H[c][s]), fmt(sp["titv"]), sp["n_TA"],
                     fmt(sp["pct_TA"], 1), sp["n_CT"], fmt(sp["max_TA"]), fmt(sp["max_af"]), sp["ge5"], sp["ge10"]])
    hdr = ["caller", "calls", "samples_with_calls", "Ti/Tv", "T>A", "T>A_%", "C>T", "max_T>A_VAF_%",
           "max_VAF_%", "calls_ge5%", "calls_ge10%"]
    write_tsv(f"{out}/cohort_summary.tsv", hdr, rows)
    md.append("## 1. Cohort (47 final BAMs, NUMT reads already removed)\n")
    md.append(f"Harmonised call definition: PASS, depth ≥{MIN_DEPTH}, {VAF_LO:.0%} ≤ VAF < {VAF_HI:.0%}, poly-C tracts masked. "
              "gnomAD v3.1 Ti/Tv = 3.51.\n")
    md.append(table(hdr, rows))

    # 2. threshold sweep per caller
    rows = []
    for c in H:
        recs = [r for s in H[c] for r in H[c][s].values()]
        for t in SWEEP:
            k = [r for r in recs if r["af"] >= t / 100]
            sp = spectrum(k)
            rows.append([LABEL[c], f"≥{t}", sp["n"], sp["ti"], sp["tv"], fmt(sp["titv"]), sp["n_TA"]])
    hdr = ["caller", "min_VAF_%", "calls", "transitions", "transversions", "Ti/Tv", "T>A"]
    write_tsv(f"{out}/cohort_threshold_sweep.tsv", hdr, rows)
    md.append("\n### Threshold sweep\n"); md.append(table(hdr, rows))

    # 3. overlap with the in-house call set
    ours_keys = {(s, k): r for s in H["ours"] for k, r in H["ours"][s].items()}
    groups = {"all in-house calls": lambda r: True,
              "in-house T>A calls": lambda r: sub(r) == "T>A",
              "in-house calls ≥5% VAF": lambda r: r["af"] >= 0.05,
              "in-house calls <5% VAF, not T>A": lambda r: r["af"] < 0.05 and sub(r) != "T>A"}
    rows, detail = [], []
    for g, test in groups.items():
        keys = [k for k, r in ours_keys.items() if test(r)]
        row = [g, len(keys)]
        for c in A:
            harm = sum(1 for s, k in keys if k in H[c].get(s, {}))
            seen = sum(1 for s, k in keys if k in A[c].get(s, {}))
            row += [harm, seen]
        rows.append(row)
    hdr = ["in-house subset", "n"] + [x for c in A for x in (f"{c}_called(harmonised)", f"{c}_any_PASS_record")]
    write_tsv(f"{out}/cohort_overlap_with_inhouse.tsv", hdr, rows)
    md.append("\n### How many in-house calls each external caller also makes\n")
    md.append("`called` uses the harmonised definition; `any PASS record` counts the caller reporting the allele at any VAF below 95%.\n")
    md.append(table(hdr, rows))

    for (s, k), r in sorted(ours_keys.items(), key=lambda x: (x[0][0], x[0][1][0])):
        line = [s, f"m.{k[0]}{r['ref']}>{k[1]}", sub(r), f"{100*r['af']:.2f}", r["dp"]]
        for c in A:
            t = A[c].get(s, {}).get(k)
            line += [f"{100*t['af']:.2f}" if t else "", "yes" if k in H[c].get(s, {}) else "no"]
        detail.append(line)
    write_tsv(f"{out}/cohort_inhouse_calls_per_caller.tsv",
              ["sample", "variant", "substitution", "inhouse_VAF_%", "depth"] +
              [x for c in A for x in (f"{c}_VAF_%", f"{c}_called")], detail)

    # 4. calls the external callers make that the in-house caller does not
    rows = []
    for c in A:
        extra = [r for s in H[c] for k, r in H[c][s].items() if (s, k) not in ours_keys]
        sp = spectrum(extra)
        reg = Counter(region(r["pos"]) for r in extra)
        rows.append([LABEL[c], sp["n"], fmt(sp["titv"]), sp["n_TA"], sp["ge5"],
                     reg["control_region"], reg["MT-RNR1"] + reg["MT-RNR2"], reg["other"]])
    hdr = ["caller", "calls_not_in_inhouse", "Ti/Tv", "T>A", "≥5%_VAF", "control_region", "rRNA", "other"]
    write_tsv(f"{out}/cohort_external_only_calls.tsv", hdr, rows)
    md.append("\n### Calls made by an external caller but not by the in-house caller\n")
    md.append(table(hdr, rows))
    return H


def numt_analysis(final_raw, unf_raw, out, md):
    md.append("\n## 2. NUMT arm — unfiltered Pass-1 BAMs vs final BAMs\n")
    md.append("Same caller, same samples; the only difference is whether competitive realignment to GRCh38 removed NUMT reads.\n")
    rows, banded = [], []
    for c in unf_raw:
        F, U = harmonise(final_raw[c]), harmonise(unf_raw[c])
        nF = sum(len(v) for v in F.values()); nU = sum(len(v) for v in U.values())
        added = [r for s in U for k, r in U[s].items() if k not in F.get(s, {})]
        sp = spectrum(added)
        reg = Counter(region(r["pos"]) for r in added)
        band = sum(500 <= r["pos"] <= 3499 for r in added)
        rows.append([LABEL[c], nU, nF, fmt(nU / nF if nF else float("nan")),
                     fmt(100 * (1 - nF / nU) if nU else float("nan"), 1), sp["n"], fmt(sp["titv"]),
                     sp["n_CT"], sp["ge5"], reg["MT-RNR1"] + reg["MT-RNR2"], fmt(100 * band / sp["n"] if sp["n"] else 0.0, 1)])
        win = Counter((r["pos"] // 500) * 500 for r in added)
        for w0, n in sorted(win.items()):
            rr = [r for r in added if (r["pos"] // 500) * 500 == w0]
            banded.append([LABEL[c], f"m.{w0}-{w0+499}", n, fmt(spectrum(rr)["titv"]) if n >= 20 else ""])
    hdr = ["caller", "calls_unfiltered", "calls_final", "fold", "%_removed_by_NUMT_step", "added_calls",
           "added_Ti/Tv", "added_C>T", "added_≥5%", "added_in_rRNA", "added_%_in_m.500-3499"]
    write_tsv(f"{out}/numt_summary.tsv", hdr, rows)
    write_tsv(f"{out}/numt_added_calls_by_500bp.tsv", ["caller", "window", "added_calls", "Ti/Tv(≥20 calls)"], banded)
    md.append(table(hdr, rows))
    md.append("\nIn-house reference from the ablation run: 1,081 → 134 calls (87.6% removed), added calls Ti/Tv 3.11, "
              "56% in m.500–3,499. Per-window counts are in `numt_added_calls_by_500bp.tsv`.\n")


def mixture_analysis(mix_dir, callers_raw, parents_raw, out, md):
    design = {r["mixture"]: r for r in csv.DictReader(open(f"{mix_dir}/design.tsv"), delimiter="\t")}
    truth = defaultdict(dict)
    for r in csv.DictReader(open(f"{mix_dir}/truth.tsv"), delimiter="\t"):
        truth[r["mixture"]][(int(r["pos"]), r["alt"])] = r
    # alleles explained by a parent's own heteroplasmy (any caller, cohort run)
    parent_het = defaultdict(set)
    for raw in parents_raw.values():
        for s, recs in harmonise(raw).items():
            parent_het[s] |= set(recs)

    rec_rows, fp_rows, bias = [], [], []
    agg = defaultdict(lambda: [0, 0])          # (caller, depth, f) -> [found, total]
    fp_agg = defaultdict(lambda: [0, 0, 0])    # (caller, depth, f) -> [unexplained, explained, mixtures]
    for c, raw in callers_raw.items():
        H = harmonise(raw)
        for m, d in design.items():
            if m not in H:
                continue
            f, depth = float(d["fraction"]), int(d["target_depth"])
            calls = H[m]
            minor_truth = {k: t for k, t in truth[m].items() if t["carrier"] == "minor"}
            found = [k for k in minor_truth if k in calls]
            agg[(c, depth, f)][0] += len(found); agg[(c, depth, f)][1] += len(minor_truth)
            for k in found:
                bias.append([c, m, depth, f, k[0], f"{100*calls[k]['af']:.3f}", f"{100*f:.1f}"])
            explained = parent_het[d["minor"]] | parent_het[d["major"]]
            fps = [k for k in calls if k not in truth[m]]
            unexp = [k for k in fps if k not in explained]
            fp_agg[(c, depth, f)][0] += len(unexp); fp_agg[(c, depth, f)][1] += len(fps) - len(unexp)
            fp_agg[(c, depth, f)][2] += 1
            for k in unexp:
                r = calls[k]
                fp_rows.append([c, m, depth, f, f"m.{k[0]}{r['ref']}>{k[1]}", sub(r), f"{100*r['af']:.2f}", r["dp"]])

    rows = []
    for (c, depth, f), (fd, tot) in sorted(agg.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
        un, ex, nm = fp_agg[(c, depth, f)]
        rows.append([LABEL[c], depth, f"{100*f:g}", tot, fd, fmt(100 * fd / tot if tot else float("nan"), 1),
                     fmt(un / nm if nm else float("nan"), 1), fmt(ex / nm if nm else float("nan"), 1)])
    hdr = ["caller", "target_depth", "true_VAF_%", "true_sites", "recovered", "recall_%",
           "unexplained_FP_per_mixture", "parent_explained_calls_per_mixture"]
    write_tsv(f"{out}/mixture_recall.tsv", hdr, rows)
    write_tsv(f"{out}/mixture_false_positives.tsv",
              ["caller", "mixture", "depth", "fraction", "variant", "substitution", "VAF_%", "depth_at_site"], fp_rows)
    write_tsv(f"{out}/mixture_vaf_estimates.tsv",
              ["caller", "mixture", "depth", "fraction", "pos", "observed_VAF_%", "expected_VAF_%"], bias)
    md.append("\n## 3. In-silico mixtures (known truth)\n")
    md.append("Recall is over positions homoplasmic in the minor sample only (true VAF = mixing fraction). "
              "A false positive is a harmonised call at a non-truth position that is not one of either parent's own heteroplasmies.\n")
    md.append(table(hdr, rows))
    fpt = Counter((r[0], r[5]) for r in fp_rows)
    if fp_rows:
        md.append("\n### Unexplained false positives by substitution\n")
        subs = sorted({s for _, s in fpt})
        cs = sorted({c for c, _ in fpt})
        md.append(table(["substitution"] + [LABEL[c] for c in cs], [[s] + [fpt[(c, s)] for c in cs] for s in subs]))
    return agg


def table(hdr, rows):
    s = "| " + " | ".join(map(str, hdr)) + " |\n|" + "---|" * len(hdr) + "\n"
    return s + "".join("| " + " | ".join(map(str, r)) + " |\n" for r in rows)


def figures(out, H_cohort, agg):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    COL = {"ours": "#2e7d74", "ours_nopon": "#2e7d74", "mutect2": "#c97b26", "mutserve": "#2a78d6"}
    STY = {"ours_nopon": "--"}
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "axes.spines.top": False,
                         "axes.spines.right": False, "font.size": 9})
    n = 2 + (1 if agg else 0)
    fig, axes = plt.subplots(1, n, figsize=(3.4 * n, 3.6))
    ax = axes[0]
    for c, H in H_cohort.items():
        recs = [r for s in H for r in H[s].values()]
        ys = [spectrum([r for r in recs if r["af"] >= t / 100])["titv"] for t in SWEEP]
        ax.plot(SWEEP, ys, marker="o", color=COL[c], label=LABEL[c])
    ax.axhline(3.51, color="#888", ls=":", lw=1); ax.text(10.4, 3.58, "gnomAD 3.51", ha="right", fontsize=7, color="#666")
    ax.set_xlim(0.5, 10.5); ax.set_ylim(bottom=0)
    ax.set_xlabel("minimum allele fraction (%)"); ax.set_ylabel("Ti/Tv"); ax.set_title("a   Ti/Tv by threshold", loc="left", fontweight="bold")
    ax = axes[1]
    for c, H in H_cohort.items():
        recs = [r for s in H for r in H[s].values()]
        ys = [sum(1 for r in recs if sub(r) == "T>A" and r["af"] >= t / 100) for t in SWEEP]
        ax.plot(SWEEP, ys, marker="s", color=COL[c], label=LABEL[c])
    ax.axvline(5, color="#444", ls=":", lw=1); ax.set_xlim(0.5, 10.5)
    ax.set_xlabel("minimum allele fraction (%)"); ax.set_ylabel("T>A calls"); ax.set_title("b   T>A calls retained", loc="left", fontweight="bold")
    if agg:
        ax = axes[2]
        for c in sorted({k[0] for k in agg}):
            for depth, mk in ((200, "o"), (500, "^")):
                pts = sorted((f, v[0] / v[1]) for (cc, d, f), v in agg.items() if cc == c and d == depth and v[1])
                if pts:
                    ax.plot([100 * p[0] for p in pts], [100 * p[1] for p in pts], marker=mk, color=COL[c],
                            ls=STY.get(c, "-"), alpha=0.9)
        ax.axvline(5, color="#444", ls=":", lw=1)
        ax.set_xscale("log"); ax.set_xticks([1, 2, 3, 5, 10, 20], ["1", "2", "3", "5", "10", "20"])
        ax.set_xlabel("true allele fraction (%)"); ax.set_ylabel("recall (%)"); ax.set_ylim(-3, 103)
        ax.set_title("c   Recall in mixtures", loc="left", fontweight="bold")
        ax.text(0.98, 0.04, "circle 200×, triangle 500×", transform=ax.transAxes, ha="right", fontsize=7, color="#666")
    from matplotlib.lines import Line2D
    keys = list(H_cohort) + ([c for c in sorted({k[0] for k in agg}) if c not in H_cohort] if agg else [])
    fig.legend([Line2D([], [], color=COL[c], ls=STY.get(c, "-"), marker="o") for c in keys], [LABEL[c] for c in keys],
               loc="lower center", ncol=len(keys), frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(f"{out}/benchmark_overview.pdf"); fig.savefig(f"{out}/benchmark_overview.png", dpi=200)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True, help="heteroplasmic_variants.tsv (cohort, published)")
    ap.add_argument("--cohort", required=True, help="dir with mutect2/ and mutserve/ from run_tools.sh on final BAMs")
    ap.add_argument("--unfiltered", help="dir with mutect2/ and mutserve/ from run_tools.sh on Pass-1 BAMs")
    ap.add_argument("--mix", help="mixtures dir (truth.tsv, design.tsv)")
    ap.add_argument("--mix-tools", help="dir with mutect2/ and mutserve/ run on mixture BAMs")
    ap.add_argument("--mix-ours", help="in-house calls on mixtures, with panel of normals")
    ap.add_argument("--mix-ours-nopon", help="in-house calls on mixtures, without panel of normals")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    tools = {"mutect2": parse_mutect2(f"{a.cohort}/mutect2"), "mutserve": parse_mutserve(f"{a.cohort}/mutserve")}
    for k, v in tools.items():
        if not v:
            sys.exit(f"no {k} output under {a.cohort}/{k}")
    ours = parse_ours_table(a.ours)
    md = ["# Heteroplasmy caller benchmark\n",
          f"Samples: in-house {len(ours)}, Mutect2 {len(tools['mutect2'])}, mutserve {len(tools['mutserve'])}.\n"]
    H = cohort_analysis(ours, tools, a.out, md)

    if a.unfiltered:
        unf = {"mutect2": parse_mutect2(f"{a.unfiltered}/mutect2"), "mutserve": parse_mutserve(f"{a.unfiltered}/mutserve")}
        numt_analysis(tools, {k: v for k, v in unf.items() if v}, a.out, md)

    agg = None
    if a.mix and a.mix_tools:
        callers = {"mutect2": parse_mutect2(f"{a.mix_tools}/mutect2"), "mutserve": parse_mutserve(f"{a.mix_tools}/mutserve")}
        if a.mix_ours: callers["ours"] = parse_ours_dir(a.mix_ours)
        if a.mix_ours_nopon: callers["ours_nopon"] = parse_ours_dir(a.mix_ours_nopon)
        agg = mixture_analysis(a.mix, {k: v for k, v in callers.items() if v}, {"ours": ours, **tools}, a.out, md)

    figures(a.out, H, agg)
    with open(f"{a.out}/benchmark_report.md", "w") as fh:
        fh.write("\n".join(md))
    print(f"  report: {a.out}/benchmark_report.md", file=sys.stderr)


if __name__ == "__main__":
    main()
