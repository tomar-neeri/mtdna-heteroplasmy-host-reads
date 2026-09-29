#!/usr/bin/env python3
"""
Remove reads Kraken2 assigned to bacteria, archaea or viruses. Keep human and
unclassified.

THE TRAP THIS AVOIDS: the standard Kraken2 database contains the human genome,
so genuine mitochondrial reads are classified as Homo sapiens. A naive
"keep only --unclassified-out" filter therefore discards ~92% of real data.
Observed directly: 39,409 of 42,857 reads removed from a high-host library, top taxon
Homo sapiens.

Validation result: with a correct taxonomic filter, only 8/42,857 reads in
a high-host library and 6/21,058 in a low-host library are bacterial, despite the latter being 90.5%
microbial as a sample. Pass 2b is therefore optional and off by default.

Usage:
  kraken_filter.py --kraken-out k.out --report k.report --bam in.bam --out out.bam
"""
import argparse, subprocess, sys, tempfile, os

EXCLUDE_ROOTS = {"Bacteria", "Archaea", "Viruses"}


def taxids_under(report_path, roots):
    """Walk the indented Kraken report and collect taxids beneath the named roots."""
    excluded, stack, active_depth = set(), [], None
    for line in open(report_path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        name_field = f[5]
        depth = (len(name_field) - len(name_field.lstrip())) // 2
        name = name_field.strip()
        taxid = f[4]

        while stack and stack[-1][0] >= depth:
            stack.pop()
        stack.append((depth, name))

        if active_depth is not None and depth <= active_depth:
            active_depth = None
        if name in roots:
            active_depth = depth
        if active_depth is not None:
            excluded.add(taxid)
    return excluded


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kraken-out", required=True, help="per-read output (--output)")
    ap.add_argument("--report", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--exclude", default=",".join(sorted(EXCLUDE_ROOTS)))
    a = ap.parse_args()

    roots = set(a.exclude.split(","))
    bad_taxids = taxids_under(a.report, roots)
    print(f"  excluding {len(bad_taxids)} taxids under {sorted(roots)}", file=sys.stderr)

    drop = set()
    total = 0
    for line in open(a.kraken_out):
        f = line.split("\t")
        if len(f) < 3:
            continue
        total += 1
        if f[0] == "C" and f[2] in bad_taxids:
            drop.add(f[1])

    print(f"  {len(drop)}/{total} reads assigned to excluded taxa "
          f"({len(drop)/total*100 if total else 0:.3f}%)", file=sys.stderr)

    if not drop:
        subprocess.run(["cp", a.bam, a.out], check=True)
    else:
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
            fh.write("\n".join(sorted(drop)) + "\n")
            names = fh.name
        try:
            with open(a.out, "wb") as out:
                subprocess.run(["samtools", "view", "-b", "-N", "^" + names, a.bam],
                               stdout=out, check=True)
        finally:
            os.unlink(names)
    subprocess.run(["samtools", "index", a.out], check=True)


if __name__ == "__main__":
    main()
