#!/usr/bin/env python3
"""
Figure 3 (c-e): macrohaplogroup assignment and marker typing.

Reads the v1.6 outputs of haplogroup.py:
  haplogroups.tsv          one row per sample
  haplogroup_markers.tsv   one row per sample x marker (BAM typing mode)

  (c) samples per macrohaplogroup, unassigned shown separately
  (d) markers typed (of 11) against median chrM depth, coloured by assignment
  (e) major-allele fraction at every typed marker; mixture would sit between
      the 0.80 cut-off and 1.0

Usage:
  plot_haplogroups.py --summary results/summary --out figures/Figure3_cde
"""
import argparse, csv, os, random
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, MUTE, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3de"
COL = {"M": "#2a78d6", "R": "#eb6834", "N (non-R)": "#1baf7a",
       "N (R undetermined)": "#1baf7a", "non-N non-M": "#4a3aa7",
       "unassigned": "#b8b7b1"}
ORDER = ["M", "R", "N (non-R)", "N (R undetermined)", "non-N non-M", "unassigned"]


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out", required=True, help="output path without extension")
    ap.add_argument("--min-major-af", type=float, default=0.80)
    a = ap.parse_args()
    random.seed(0)

    hg = read(os.path.join(a.summary, "haplogroups.tsv"))
    mk_path = os.path.join(a.summary, "haplogroup_markers.tsv")
    mk = read(mk_path) if os.path.exists(mk_path) else []
    n = len(hg)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
        "axes.edgecolor": MUTE, "axes.labelcolor": INK2, "xtick.color": INK2,
        "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8})
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.6),
                           gridspec_kw={"width_ratios": [1, 1.2, 1.2]})

    # (c) counts
    c = Counter(r["macrohaplogroup"] for r in hg)
    cats = [k for k in ORDER if c.get(k)] + [k for k in c if k not in ORDER]
    y = range(len(cats))[::-1]
    for yi, k in zip(y, cats):
        ax[0].barh(yi, c[k], height=0.6, color=COL.get(k, MUTE),
                   edgecolor="white", lw=2)
        ax[0].text(c[k] + n * 0.015, yi, f"{c[k]} ({c[k]/n:.1%})",
                   va="center", fontsize=8, color=INK)
    ax[0].set_yticks(list(y)); ax[0].set_yticklabels(cats)
    ax[0].set_xlim(0, max(c.values()) * 1.45)
    ax[0].set_xlabel("Samples")
    ax[0].grid(axis="x", color=GRID, lw=0.6); ax[0].set_axisbelow(True)
    ax[0].set_title(f"c   Macrohaplogroup (n = {n})", loc="left",
                    fontsize=9.5, color=INK, fontweight="bold")

    # (d) markers typed vs depth
    for r in hg:
        typed = int(r["markers_typed"].split("/")[0])
        dep = float(r["median_depth"]) if r["median_depth"] else None
        if dep is None:
            continue
        typed_j = typed + random.uniform(-0.18, 0.18)   # separate overlapping samples
        ax[1].scatter(dep, typed_j, s=46, color=COL.get(r["macrohaplogroup"], MUTE),
                      edgecolor="white", lw=1.2, zorder=3)
        if r["macrohaplogroup"] == "unassigned" or "mixture" in r["flags"]:
            ax[1].annotate(r["sample"], (dep, typed_j), xytext=(6, -3),
                           textcoords="offset points", fontsize=7.5, color=INK)
    ax[1].set_xscale("log")
    ax[1].axvline(50, color=INK2, ls=(0, (3, 3)), lw=1)
    ax[1].text(50, 11.6, " 50× calling threshold", fontsize=7.5, color=INK2)
    ax[1].set_ylim(-0.5, 12.3); ax[1].set_yticks(range(0, 12, 2))
    ax[1].set_xlabel("Median chrM depth (×, log scale)")
    ax[1].set_ylabel("Markers typed (of 11)")
    ax[1].grid(color=GRID, lw=0.6); ax[1].set_axisbelow(True)
    ax[1].set_title("d   Marker typing against depth", loc="left",
                    fontsize=9.5, color=INK, fontweight="bold")
    handles = [plt.Line2D([], [], marker="o", ls="", color=COL.get(k, MUTE),
                          mec="white", ms=7, label=k) for k in cats]
    ax[1].legend(handles=handles, fontsize=7, frameon=False, loc="lower right")

    # (e) major-allele fraction
    afs = [float(r["major_allele_fraction"]) for r in mk
           if r["major_allele_fraction"] not in ("", None)]
    if afs:
        lo = min(0.5, min(afs) - 0.02)
        bins = [lo + i * (1.0 - lo) / 40 for i in range(41)]
        ax[2].hist(afs, bins=bins, color="#2a78d6", edgecolor="white", lw=0.8)
        ax[2].set_yscale("log")
        ax[2].axvline(a.min_major_af, color=INK2, ls=(0, (3, 3)), lw=1)
        ax[2].text(a.min_major_af, 0.03, f"cut-off {a.min_major_af:.2f} ",
                   transform=ax[2].get_xaxis_transform(), ha="right", va="bottom",
                   fontsize=7.5, color=INK2)
        n_mixed = sum(1 for v in afs if v < a.min_major_af)
        ax[2].text(0.03, 0.97, f"{len(afs)} markers with reads\n{n_mixed} below cut-off",
                   transform=ax[2].transAxes, ha="left", va="top", fontsize=7.5, color=INK)
        ax[2].set_title("e   Major-allele fraction at markers",
                        loc="left", fontsize=9.5, color=INK, fontweight="bold")
    else:
        ax[2].text(0.5, 0.5, "no per-marker allele fractions\n(run haplogroup.py with --bams)",
                   ha="center", va="center", transform=ax[2].transAxes, color=INK2)
        ax[2].set_title("e   Major-allele fraction at markers", loc="left",
                        fontsize=9.5, color=INK, fontweight="bold")
    ax[2].set_xlabel("Major-allele fraction")
    ax[2].set_ylabel("Markers (log scale)")
    ax[2].grid(axis="y", color=GRID, lw=0.6); ax[2].set_axisbelow(True)

    fig.tight_layout(w_pad=2.5)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=300 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    print(f"  written: {a.out}.pdf / .png")


if __name__ == "__main__":
    main()
