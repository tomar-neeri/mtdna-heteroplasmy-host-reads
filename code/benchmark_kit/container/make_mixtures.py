#!/usr/bin/env python3
"""
Build in-silico heteroplasmy by mixing reads from two cohort samples of different
macrohaplogroups, and write the truth set. Runs INSIDE the container.

An M and an R sample differ homoplasmically at dozens of positions. Mixing a
fraction f of reads from the minor sample into the major sample creates a true
heteroplasmy at every such position: expected allele fraction f where only the
minor sample carries the alternate allele, and 1-f where only the major does.
This gives a ground truth that the 5% threshold was not derived from.

Subsampling uses samtools view -s, which keeps or drops read pairs together.
Read names cannot collide between samples (distinct clusters on the flowcell).

  make_mixtures.py --bam-dir /bams --homoplasmic /summary/homoplasmic_variants.tsv \
                   --pairs 425:418,414:402,389:393 --out /bench/mixtures
"""
import argparse, csv, os, re, subprocess, sys

FRACTIONS = [0.01, 0.02, 0.03, 0.05, 0.10, 0.20]
DEPTHS = [200, 500]
REPS = [1, 2]
MASK = set(range(302, 316)) | set(range(16184, 16194))
L = 16569


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def find_bam(bam_dir, sid):
    hits = [f for f in os.listdir(bam_dir) if f.endswith(".bam") and re.search(rf"(^|_){sid}(_|\.)", f)]
    if len(hits) != 1:
        sys.exit(f"sample {sid}: expected one BAM in {bam_dir}, found {hits}")
    return os.path.join(bam_dir, hits[0])


def read_stats(bam, contig, mapq):
    """Reads on chrM at MAPQ >= mapq, and mean aligned length (from 20k reads)."""
    n = int(sh(f"samtools view -c -q {mapq} {bam} {contig}").strip())
    lens = sh(f"samtools view -q {mapq} {bam} {contig} | head -20000 | awk '{{print length($10)}}'").split()
    mean_len = sum(map(int, lens)) / max(len(lens), 1)
    return n, mean_len


def load_homoplasmic(path):
    hom = {}
    with open(path, newline="") as fh:
        delim = "\t" if "\t" in fh.readline() else ","
        fh.seek(0)
        for r in csv.DictReader(fh, delimiter=delim):
            s = re.sub(r".*_(\d{3})_S\d+$", r"\1", r["sample"])
            hom.setdefault(s, {})[int(r["pos"])] = (r["ref"], r["alt"])
    return hom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam-dir", required=True)
    ap.add_argument("--homoplasmic", required=True, help="homoplasmic_variants.tsv/.csv from the pipeline")
    ap.add_argument("--pairs", default="425:418,414:402,389:393",
                    help="minor:major sample pairs; minor is diluted into major")
    ap.add_argument("--out", required=True)
    ap.add_argument("--contig", default="chrM")
    ap.add_argument("--min-mapq", type=int, default=30)
    a = ap.parse_args()

    os.makedirs(os.path.join(a.out, "bams"), exist_ok=True)
    hom = load_homoplasmic(a.homoplasmic)
    truth_rows, design_rows = [], []

    for pair in a.pairs.split(","):
        minor, major = pair.split(":")
        for s in (minor, major):
            if s not in hom:
                sys.exit(f"sample {s} has no homoplasmic variants in {a.homoplasmic}")
        bmin, bmaj = find_bam(a.bam_dir, minor), find_bam(a.bam_dir, major)
        nmin, lmin = read_stats(bmin, a.contig, a.min_mapq)
        nmaj, lmaj = read_stats(bmaj, a.contig, a.min_mapq)
        dmin, dmaj = nmin * lmin / L, nmaj * lmaj / L
        print(f"  pair {minor} (minor, ~{dmin:.0f}x) into {major} (major, ~{dmaj:.0f}x)", file=sys.stderr)

        # truth positions: homoplasmic in exactly one of the two samples
        hmin, hmaj = hom[minor], hom[major]
        sites = []
        for pos in sorted(set(hmin) | set(hmaj)):
            if pos in MASK:
                continue
            if pos in hmin and pos in hmaj:
                continue                      # shared (or multi-allelic): no mixture signal
            if pos in hmin:
                sites.append((pos, *hmin[pos], "minor"))
            else:
                sites.append((pos, *hmaj[pos], "major"))
        print(f"    {len(sites)} informative positions "
              f"({sum(x[3]=='minor' for x in sites)} minor-only, {sum(x[3]=='major' for x in sites)} major-only)",
              file=sys.stderr)

        for depth in DEPTHS:
            for f in FRACTIONS:
                for rep in REPS:
                    # read fractions that give f of coverage from the minor sample at the target depth
                    pmin = f * depth / dmin
                    pmaj = (1 - f) * depth / dmaj
                    if pmin > 1 or pmaj > 1:
                        print(f"    skip depth {depth} f {f}: needs more reads than available", file=sys.stderr)
                        continue
                    mid = f"mix_{minor}in{major}_d{depth}_f{int(round(f*100)):02d}_r{rep}"
                    outbam = os.path.join(a.out, "bams", f"{mid}.bam")
                    if not os.path.exists(outbam + ".bai"):
                        seed = 11 * rep + int(f * 1000) + depth
                        tmin, tmaj = outbam + ".min.bam", outbam + ".maj.bam"
                        def sub(seed_, p_):   # -s INT.FRAC; a fraction of ~1 means keep everything
                            return "" if p_ >= 0.9995 else f"-s {seed_ + min(p_, 0.9994):.6f}"
                        sh(f"samtools view -b -q {a.min_mapq} {sub(seed, pmin)} {bmin} {a.contig} > {tmin}")
                        sh(f"samtools view -b -q {a.min_mapq} {sub(seed + 7, pmaj)} {bmaj} {a.contig} > {tmaj}")
                        sh(f"samtools merge -f -o {outbam}.unsorted.bam {tmin} {tmaj}")
                        sh(f"samtools sort -o {outbam} {outbam}.unsorted.bam && samtools index {outbam}")
                        for t in (tmin, tmaj, outbam + ".unsorted.bam"):
                            os.remove(t)
                    design_rows.append([mid, minor, major, depth, f, rep, f"{pmin:.5f}", f"{pmaj:.5f}"])
                    for pos, ref, alt, who in sites:
                        expected = f if who == "minor" else 1 - f
                        truth_rows.append([mid, pos, ref, alt, who, f"{expected:.4f}"])
            print(f"    depth {depth}: done", file=sys.stderr)

    with open(os.path.join(a.out, "truth.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mixture", "pos", "ref", "alt", "carrier", "expected_af"])
        w.writerows(truth_rows)
    with open(os.path.join(a.out, "design.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mixture", "minor", "major", "target_depth", "fraction", "rep", "p_minor", "p_major"])
        w.writerows(design_rows)
    print(f"  {len(design_rows)} mixtures, {len(truth_rows)} truth records -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
