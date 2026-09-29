#!/usr/bin/env python3
"""
Consolidate per-sample calls into cohort tables, with protein-coding annotation.

Outputs (into --outdir):
  heteroplasmic_variants.tsv   every passing call
  homoplasmic_variants.tsv     haplogroup backbone per sample
  variant_recurrence.tsv       positions shared across samples
  region_distribution.tsv      counts per mtDNA region
  variants_annotated.tsv       gene, codon, amino-acid consequence
  cohort_overview.tsv          one row per sample

Recurrence matters: genuine heteroplasmy is largely private to an individual.
A position appearing in several samples is either haplogroup-linked leakage or
a systematic artifact, and is flagged as such.

Optional annotation against static reference tables (--mitomap, --gnomad):
these are used for FILTERING, not clinical interpretation. A call at 3% VAF
that is also a known population polymorphism is more likely lineage leakage or
a mutational hotspot than private heteroplasmy.

Deliberately NOT annotated against ClinVar as a primary field. Pathogenic mtDNA
variants generally require 60-90% heteroplasmy to produce a phenotype, so a
"pathogenic" label on a 3% call is accurate and substantively misleading.
Population frequency is the useful signal at these levels.

Usage:
  collect_variants.py --calls results/calls --ref rCRS.fasta --outdir results/summary
  collect_variants.py ... --mitomap mitomap_polymorphisms.tsv --gnomad gnomad_mt.tsv
"""
import argparse, csv, glob, json, os, re, sys
from collections import defaultdict
from statistics import median

# vertebrate mitochondrial genetic code
_B, _A = "TCAG", ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG")
CODE = {a + b + c: _A[i] for i, (a, b, c) in enumerate(
    (x, y, z) for x in _B for y in _B for z in _B)}
CODE.update({"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W"})

GENES = [("MT-ND1", 3307, 4262), ("MT-ND2", 4470, 5511), ("MT-CO1", 5904, 7445),
         ("MT-CO2", 7586, 8269), ("MT-ATP8", 8366, 8572), ("MT-ATP6", 8527, 9207),
         ("MT-CO3", 9207, 9990), ("MT-ND3", 10059, 10404), ("MT-ND4L", 10470, 10766),
         ("MT-ND4", 10760, 12137), ("MT-ND5", 12337, 14148), ("MT-ND6", 14149, 14673),
         ("MT-CYB", 14747, 15887)]
# MT-ND6 is transcribed from the light strand
REVERSE = {"MT-ND6"}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def load_gnomad(path):
    """Load the gnomAD v3.1 chrM reduced-annotations TSV.

    Schema: chromosome position ref alt filters AC_hom AC_het AF_hom AF_het AN
            max_observed_heteroplasmy

    Four fields matter here:
      AF_hom   population frequency as a fixed (homoplasmic) variant
      AF_het   population frequency at heteroplasmy -- the relevant comparator
               for low-VAF calls like ours
      filters  gnomAD's own QC flags. "artifact_prone_site" and
               "common_low_heteroplasmy" mark positions gnomAD themselves
               identified as NUMT-misalignment false positives across 56,434
               genomes. A call of ours at such a position is suspect even after
               Pass 2a.
      max_observed_heteroplasmy  highest VAF gnomAD saw at this position
    """
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        need = {"position", "ref", "alt"}
        if not rdr.fieldnames or not need.issubset(set(rdr.fieldnames)):
            print(f"  [warn] gnomAD: unexpected columns {rdr.fieldnames}",
                  file=sys.stderr)
            return {}
        def f(v):
            try:
                return float(v)
            except (TypeError, ValueError):
                return None
        for r in rdr:
            try:
                p = int(r["position"])
            except (TypeError, ValueError):
                continue
            out[(p, r["ref"].strip().upper(), r["alt"].strip().upper())] = {
                "af_hom": f(r.get("AF_hom")),
                "af_het": f(r.get("AF_het")),
                "filters": (r.get("filters") or "").strip(),
                "max_het": f(r.get("max_observed_heteroplasmy")),
            }
    print(f"  gnomAD    : {len(out):,} variants loaded", file=sys.stderr)
    return out


def load_mitomap(path):
    """MITOMAP GenBank-frequency table (optional). Column names vary by export,
    so match on substring."""
    if not path or not os.path.exists(path):
        return {}
    out, delim = {}, ("\t" if path.endswith((".tsv", ".txt")) else ",")
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rdr = csv.DictReader(fh, delimiter=delim)
        if not rdr.fieldnames:
            return {}
        def find(want):
            for c in rdr.fieldnames:
                if want in c.lower().replace("_", "").replace(" ", ""):
                    return c
            return None
        cp, cr, ca = find("pos"), find("ref"), find("alt")
        cf = find("gbfreq") or find("freq")
        if not all([cp, cr, ca]):
            print(f"  [warn] MITOMAP: no pos/ref/alt in {rdr.fieldnames[:8]}",
                  file=sys.stderr)
            return {}
        for r in rdr:
            try:
                p = int(re.sub(r"[^0-9]", "", str(r[cp])))
            except (TypeError, ValueError):
                continue
            try:
                v = float(str(r.get(cf, "")).strip().rstrip("%"))
            except (TypeError, ValueError):
                v = None
            out[(p, str(r[cr]).strip().upper()[:1],
                 str(r[ca]).strip().upper()[:1])] = v
    print(f"  MITOMAP   : {len(out):,} entries loaded", file=sys.stderr)
    return out


def short_name(path):
    n = os.path.basename(path).split(".calls.tsv")[0]
    return re.sub(r"_S\d+$", "", n.replace("24D214-5G_", ""))


def load_ref(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
    if len(seq) != 16569:
        sys.exit(f"reference length {len(seq)} != 16569")
    return seq


def annotate(pos, ref_b, alt_b, refseq):
    for g, s, e in GENES:
        if s <= pos <= e:
            off = pos - s
            ci, cp = off // 3, off % 3
            cs = s + ci * 3
            codon = refseq[cs - 1:cs + 2]
            if len(codon) != 3:
                return g, "-", "-", "-", "coding"
            mut = list(codon); mut[cp] = alt_b; mut = "".join(mut)
            if g in REVERSE:
                codon = "".join(COMP.get(b, b) for b in reversed(codon))
                mut = "".join(COMP.get(b, b) for b in reversed(mut))
            a1, a2 = CODE.get(codon, "?"), CODE.get(mut, "?")
            if a1 == a2:
                eff = "synonymous"
            elif a2 == "*":
                eff = "nonsense"
            elif a1 == "*":
                eff = "stop_loss"
            else:
                eff = "missense"
            return g, str(ci + 1), f"{codon}>{mut}", f"{a1}{ci+1}{a2}", eff
    if pos >= 16024 or pos <= 576:
        return "control_region", "-", "-", "-", "non-coding"
    if 648 <= pos <= 1601:
        return "MT-RNR1", "-", "-", "-", "rRNA"
    if 1671 <= pos <= 3229:
        return "MT-RNR2", "-", "-", "-", "rRNA"
    return "tRNA/other", "-", "-", "-", "non-coding"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mitomap", default=None,
                    help="MITOMAP polymorphisms table (optional). Used for "
                         "population-frequency filtering, not interpretation.")
    ap.add_argument("--gnomad", default=None,
                    help="gnomAD v3.1 mtDNA table (optional).")
    ap.add_argument("--lineage-flags", default=None,
                    help="heteroplasmy_lineage_flags.tsv from haplogroup.py; "
                         "marks calls at haplogroup-defining positions")
    ap.add_argument("--common-af", type=float, default=0.001,
                    help="AF above which a call is flagged as a known "
                         "polymorphism rather than private heteroplasmy")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    refseq = load_ref(a.ref)
    # positions that define a haplogroup anywhere in the cohort: phylogenetic
    # hotspots, where a low-VAF call is more likely recurrent mutation than
    # private heteroplasmy
    lineage, haplo = set(), {}
    if a.lineage_flags and os.path.exists(a.lineage_flags):
        with open(a.lineage_flags, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                haplo[r["sample"]] = r.get("haplogroup", "")
                if r.get("lineage_position") == "YES":
                    lineage.add(int(r["pos"]))
        print(f"  lineage   : {len(lineage):,} haplogroup-defining positions",
              file=sys.stderr)

    mitomap = load_mitomap(a.mitomap)
    gnomad  = load_gnomad(a.gnomad)
    files = sorted(glob.glob(os.path.join(a.calls, "*.calls.tsv")))
    if not files:
        sys.exit(f"no *.calls.tsv in {a.calls}")

    het, hom = [], []
    for f in files:
        s = short_name(f)
        for r in csv.DictReader(open(f), delimiter="\t"):
            rec = dict(r); rec["sample"] = s
            if r["class"] == "heteroplasmic" and r["PASS"] == "True":
                het.append(rec)
            elif r["class"] == "homoplasmic":
                hom.append(rec)

    def write(name, header, rows):
        p = os.path.join(a.outdir, name)
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t"); w.writerow(header)
            for r in rows:
                w.writerow(r)
        return p

    write("heteroplasmic_variants.tsv",
          ["sample","pos","ref","alt","depth","alt_reads","vaf","fwd_alt","rev_alt",
           "p_error","strand_bias_p","region"],
          [[r["sample"],r["pos"],r["ref"],r["alt"],r["depth"],r["alt_reads"],r["vaf"],
            r["fwd_alt"],r["rev_alt"],r["p_error"],r["strand_bias_p"],r["region"]]
           for r in het])

    write("homoplasmic_variants.tsv",
          ["sample","pos","ref","alt","depth","vaf","region"],
          [[r["sample"],r["pos"],r["ref"],r["alt"],r["depth"],r["vaf"],r["region"]]
           for r in hom])

    # --- annotation ---
    ann = []
    for r in het:
        pos = int(r["pos"])
        key = (pos, r["ref"].upper(), r["alt"].upper())
        g, codon, cch, aa, eff = annotate(pos, r["ref"], r["alt"], refseq)
        mm = mitomap.get(key)
        gd = gnomad.get(key)

        af_hom = gd["af_hom"] if gd else None
        af_het = gd["af_het"] if gd else None
        gflt   = gd["filters"] if gd else ""
        maxhet = gd["max_het"] if gd else None

        # status, most-actionable first
        status = ""
        if gnomad or mitomap:
            flags = gflt.lower()
            if "artifact_prone" in flags:
                status = "GNOMAD_ARTIFACT_SITE"
            elif "common_low_heteroplasmy" in flags:
                status = "GNOMAD_NUMT_FP_SITE"
            elif af_hom is not None and af_hom >= a.common_af:
                status = "KNOWN_POLYMORPHISM"
            elif isinstance(mm, float) and mm >= a.common_af * 100:
                status = "KNOWN_POLYMORPHISM"
            elif af_het is not None and af_het > 0:
                status = "seen_at_heteroplasmy"
            elif gd is not None or mm is not None:
                status = "rare"
            else:
                status = "novel"

        fmt = lambda v, p=6: (f"{v:.{p}f}" if isinstance(v, float) else "")
        lin = "YES" if pos in lineage else ("no" if lineage else "")
        hg = haplo.get(r["sample"], "")
        ann.append([r["sample"], hg, lin, f"m.{pos}{r['ref']}>{r['alt']}",
                    f"{float(r['vaf'])*100:.2f}", r["depth"], r["region"],
                    g, codon, cch, aa, eff,
                    fmt(af_hom), fmt(af_het), fmt(maxhet, 3), gflt,
                    fmt(mm, 4) if isinstance(mm, float) else "",
                    status])
    write("variants_annotated.tsv",
          ["sample","haplogroup","lineage_position","variant","vaf_pct","depth",
           "region","gene","codon","codon_change","aa_change","effect",
           "gnomad_af_hom","gnomad_af_het","gnomad_max_heteroplasmy",
           "gnomad_filters","mitomap_gb_freq_pct","population_status"], ann)

    # --- recurrence ---
    grp = defaultdict(list)
    for r in het:
        grp[(int(r["pos"]), f"{r['ref']}>{r['alt']}")].append(r)
    rec = []
    for (pos, chg), rs in sorted(grp.items(), key=lambda kv: (-len(kv[1]), kv[0][0])):
        vafs = [float(x["vaf"]) for x in rs]
        flag = ("PRIVATE" if len(rs) == 1 else
                "RECURRENT — check haplogroup / artifact")
        rec.append([len(rs), pos, chg, ",".join(sorted(x["sample"] for x in rs)),
                    f"{median(vafs)*100:.2f}", flag])
    write("variant_recurrence.tsv",
          ["n_samples","pos","change","samples","median_vaf_pct","note"], rec)

    # --- region distribution ---
    reg = defaultdict(list)
    for r in het:
        reg[r["region"]].append(r)
    write("region_distribution.tsv",
          ["region","n_variants","n_samples","median_vaf_pct"],
          [[k, len(v), len({x["sample"] for x in v}),
            f"{median([float(x['vaf']) for x in v])*100:.2f}"]
           for k, v in sorted(reg.items(), key=lambda kv: -len(kv[1]))])

    # --- per-sample overview ---
    ov = []
    for f in files:
        s = short_name(f)
        sj = f.replace(".calls.tsv", ".summary.json")
        d = json.load(open(sj)) if os.path.exists(sj) else {}
        h = [r for r in het if r["sample"] == s]
        cf = d.get("callable_fraction_at_vaf", {})
        ov.append([s, f"{d.get('mean_depth',0):.0f}", f"{d.get('median_depth',0):.0f}",
                   d.get("n_homoplasmic",""), d.get("n_heteroplasmic_raw",""),
                   len(h), cf.get("3%",""), cf.get("5%",""),
                   sum(1 for r in h if r["region"]=="protein_coding"),
                   sum(1 for r in h if r["region"]=="control_region")])
    write("cohort_overview.tsv",
          ["sample","mean_depth","median_depth","n_homoplasmic","n_het_raw",
           "n_het_pass","callable_3pct","callable_5pct","het_coding","het_control"], ov)

    # --- report ---
    n_lin = sum(1 for r in ann if r[2] == "YES")
    if lineage:
        print(f"  at lineage positions : {n_lin}/{len(ann)} "
              f"({n_lin/max(len(ann),1):.1%}) -- phylogenetic hotspots, "
              f"lower confidence", file=sys.stderr)
    n_rec = sum(1 for r in rec if r[0] > 1)
    print(f"  samples            : {len(files)}")
    print(f"  heteroplasmic calls: {len(het)}")
    print(f"  homoplasmic calls  : {len(hom)}")
    print(f"  distinct positions : {len(grp)}")
    print(f"  recurrent (>1 sample): {n_rec}"
          + ("  <- inspect these" if n_rec else "  (all private, as expected)"))
    eff = defaultdict(int)
    for r in ann:
        eff[r[11]] += 1
    print("  effects            : " + ", ".join(f"{k} {v}" for k, v in
                                                sorted(eff.items(), key=lambda x: -x[1])))
    if mitomap or gnomad:
        pop = defaultdict(int)
        for r in ann:
            pop[r[17]] += 1
        print("  population status  : " + ", ".join(f"{k} {v}" for k, v in
                                                    sorted(pop.items(), key=lambda x: -x[1])))
        for k, msg in [
            ("GNOMAD_ARTIFACT_SITE",
             "gnomAD flags this position as artifact-prone"),
            ("GNOMAD_NUMT_FP_SITE",
             "gnomAD identified this as a NUMT false-positive site across "
             "56,434 genomes -- suspect even after Pass 2a"),
            ("KNOWN_POLYMORPHISM",
             f"known population variant above AF {a.common_af:.1%}; lineage "
             f"or hotspot rather than private heteroplasmy")]:
            n = pop.get(k, 0)
            if n:
                print(f"  -> {n} call(s): {msg}")
    else:
        print("  population status  : not annotated "
              "(pass --mitomap and/or --gnomad to enable)")


if __name__ == "__main__":
    main()
