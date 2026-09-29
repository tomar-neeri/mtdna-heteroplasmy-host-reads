#!/usr/bin/env python3
"""
Assign macrohaplogroups from the allele observed at 11 diagnostic positions,
and flag heteroplasmic calls at lineage positions.

v1.6 changes (from v1.5)
------------------------
1. Marker polarity corrected. v1.5 listed 8701G, 9540C, 10398G, 10873C, 15301A
   as "N markers". Those are the ancestral L3/M states: the rCRS (haplogroup H,
   inside N and R) carries the N-derived bases at these positions, and every
   M genome differs from the rCRS there (Rajkumar et al. 2005, BMC Evol Biol 5:26).
   A sample "hitting" that set is therefore OUTSIDE N, not in it.
2. Genuine non-R N lineages are now detectable, from the R-defining positions
   12705 and 16223 (non-R N retains 12705T/16223T; R carries the rCRS bases).
3. No assignment by default. v1.5 returned R whenever no derived marker was
   found, so a sample with no data was labelled R. v1.6 reads the base at each
   marker from the BAM and reports "unassigned" unless the markers are covered
   and point one way.
4. The L0-L2 marker list was removed: it was not verified against PhyloTree,
   and an Indian cohort needs no sub-L resolution. A sample outside N that
   lacks M markers is reported as "non-N, non-M (L3/other)".

DECISION TREE (all states read relative to the rCRS)
  N block  8701 9540 10398 10873 15301   rCRS base = inside N
    majority non-rCRS  -> outside N
        M block 489 10400 14783 15043: >=2 derived  -> M
        otherwise                                   -> non-N, non-M (L3/other)
    majority rCRS      -> inside N
        R block 12705 (primary), 16223 (support)
          rCRS base  -> R
          derived    -> N (non-R)
          uncovered  -> N (R status undetermined)
    tie / too few covered -> unassigned

A marker counts as covered when >= --min-marker-depth reads pass quality
filters AND the major allele is >= --min-major-af of them. Below that fraction
the marker is "mixed", reported separately, because intermediate fractions at
lineage markers are the signature of sample mixture.

WHY NOT HAPLOGREP3
HaploGrep3 downloads PhyloTree on first use and fails silently offline. The
macrohaplogroup needs only these 11 stable markers. Subclades are not attempted.

WHY LINEAGE POSITIONS MATTER
Haplogroup-defining positions have mutated repeatedly across the phylogeny, so
a low-VAF call there is more likely recurrent mutation than private
heteroplasmy. Such calls are FLAGGED, not removed.

Usage:
  haplogroup.py --calls results/calls --bams results/chrM/*.chrM.bam \\
                --out-dir results/summary
  (without --bams the script falls back to homoplasmic calls, which cannot
   distinguish "reference base" from "not covered"; see --help)
"""
import argparse, csv, glob, json, os, re, sys
from collections import Counter

# ---------------------------------------------------------------------------
# Markers: position -> base that differs from the rCRS (1-based coordinates)
# ---------------------------------------------------------------------------
N_BLOCK = {8701: "G", 9540: "C", 10398: "G", 10873: "C", 15301: "A"}  # non-rCRS = outside N
M_BLOCK = {489: "C", 10400: "T", 14783: "C", 15043: "A"}               # non-rCRS = M
R_BLOCK = {12705: "T", 16223: "T"}                                     # non-rCRS = non-R
RCRS = {8701: "A", 9540: "T", 10398: "A", 10873: "T", 15301: "G",
        489: "T", 10400: "C", 14783: "T", 15043: "G",
        12705: "C", 16223: "C"}
MARKERS = {**N_BLOCK, **M_BLOCK, **R_BLOCK}          # 11 positions

# Phylogenetically recurrent positions (hypervariable / homoplasy-prone)
RECURRENT = {
    16093, 16129, 16183, 16189, 16192, 16311, 16362, 16519,
    152, 195, 204, 207, 228, 263, 310, 709, 750,
    1438, 2706, 3010, 4769, 8860, 9540, 10398, 11719, 12705,
    14766, 15326, 16126, 16223, 16278, 16294,
}


def short_name(p):
    n = os.path.basename(p)
    for suf in (".calls.tsv", ".chrM.bam", ".summary.json"):
        n = n.split(suf)[0]
    return re.sub(r"_S\d+$", "", n.replace("24D214-5G_", ""))


# ---------------------------------------------------------------------------
# Marker states
# ---------------------------------------------------------------------------
def states_from_bam(bam_path, contig, min_bq, min_mq, min_depth, min_major):
    """Read the base at each marker. Returns {pos: (state, depth, major_af)}.

    state: 'ref' | 'alt' | 'other' | 'mixed' | 'nocov'
    """
    import pysam
    out = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for pos in MARKERS:
            c = Counter()
            for col in bam.pileup(contig, pos - 1, pos, truncate=True,
                                  min_base_quality=min_bq, max_depth=1_000_000,
                                  stepper="nofilter"):
                for pr in col.pileups:
                    if pr.is_del or pr.is_refskip:
                        continue
                    aln = pr.alignment
                    if (aln.mapping_quality < min_mq or aln.is_secondary
                            or aln.is_supplementary or aln.is_duplicate):
                        continue
                    c[aln.query_sequence[pr.query_position]] += 1
            dep = sum(c.values())
            if dep < min_depth:
                out[pos] = ("nocov", dep, None)
                continue
            base, n = c.most_common(1)[0]
            af = n / dep
            if af < min_major:
                out[pos] = ("mixed", dep, af)
            elif base == RCRS[pos]:
                out[pos] = ("ref", dep, af)
            elif base == MARKERS[pos]:
                out[pos] = ("alt", dep, af)
            else:
                out[pos] = ("other", dep, af)
    return out


def states_from_calls(hom, median_depth, min_depth_call):
    """Fallback when no BAMs are given.

    A homoplasmic call at a marker is 'alt'. Absence of a call is only taken
    as 'ref' when the sample's median depth reaches the calling threshold;
    otherwise the marker is 'nocov'. This cannot see per-position coverage,
    so BAM mode is preferred.
    """
    callable_ = median_depth is not None and median_depth >= min_depth_call
    out = {}
    for pos, alt in MARKERS.items():
        if pos in hom:
            out[pos] = ("alt" if hom[pos] == alt else "other", None, None)
        else:
            out[pos] = ("ref" if callable_ else "nocov", None, None)
    return out


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------
def classify(st, min_block_covered=3):
    """Returns (call, reason)."""
    n_alt = sum(1 for p in N_BLOCK if st[p][0] == "alt")
    n_ref = sum(1 for p in N_BLOCK if st[p][0] == "ref")
    m_alt = sum(1 for p in M_BLOCK if st[p][0] == "alt")
    m_cov = sum(1 for p in M_BLOCK if st[p][0] in ("ref", "alt"))

    if n_alt + n_ref < min_block_covered:
        return "unassigned", f"N block covered at {n_alt + n_ref}/5 (<{min_block_covered})"
    if n_alt == n_ref:
        return "unassigned", f"N block conflicting ({n_alt} non-rCRS vs {n_ref} rCRS)"

    if n_alt > n_ref:                                   # outside N
        if m_alt >= 2:
            return "M", f"outside N ({n_alt}/{n_alt+n_ref}); M markers {m_alt}/{m_cov}"
        if m_cov < 2:
            return "unassigned", f"outside N but M block covered at {m_cov}/4"
        return "non-N non-M", f"outside N ({n_alt}/{n_alt+n_ref}); M markers {m_alt}/{m_cov}"

    # inside N
    if m_alt >= 2:
        return "unassigned", f"inside N ({n_ref}/{n_alt+n_ref}) but M markers {m_alt}/{m_cov}"
    s12705, s16223 = st[12705][0], st[16223][0]
    if s12705 in ("ref", "alt"):
        r_state = s12705
        why = f"12705 {'rCRS' if s12705 == 'ref' else 'non-rCRS'}"
    elif s16223 in ("ref", "alt"):
        r_state = s16223
        why = f"12705 not covered; 16223 {'rCRS' if s16223 == 'ref' else 'non-rCRS'}"
    else:
        return "N (R undetermined)", f"inside N ({n_ref}/{n_alt+n_ref}); 12705 and 16223 not covered"
    hg = "R" if r_state == "ref" else "N (non-R)"
    return hg, f"inside N ({n_ref}/{n_alt+n_ref}); {why}"


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, help="directory of *.calls.tsv")
    ap.add_argument("--bams", nargs="*", default=None,
                    help="final chrM BAMs; strongly recommended")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--contig", default="chrM")
    ap.add_argument("--vaf-homoplasmy", type=float, default=0.95)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--min-marker-depth", type=int, default=10,
                    help="reads needed to type a marker (lower than the 50x "
                         "calling threshold: reading a fixed base needs less)")
    ap.add_argument("--min-major-af", type=float, default=0.80,
                    help="major-allele fraction below which a marker is 'mixed'")
    ap.add_argument("--min-depth-call", type=int, default=50,
                    help="fallback mode only: median depth for absence = rCRS")
    ap.add_argument("--min-samples-lineage", type=int, default=2)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(a.calls, "*.calls.tsv")))
    if not files:
        sys.exit(f"no *.calls.tsv in {a.calls}")

    hom, het, med = {}, {}, {}
    for f in files:
        s = short_name(f)
        hom[s], het[s] = {}, []
        for r in csv.DictReader(open(f), delimiter="\t"):
            if r.get("PASS") != "True":
                continue
            try:
                vaf = float(r["vaf"]); pos = int(r["pos"])
            except (KeyError, ValueError):
                continue
            if vaf >= a.vaf_homoplasmy:
                hom[s][pos] = r["alt"]
            else:
                het[s].append((pos, r["ref"], r["alt"], vaf))
        sj = f.replace(".calls.tsv", ".summary.json")
        med[s] = json.load(open(sj)).get("median_depth") if os.path.exists(sj) else None

    bam_of = {short_name(b): b for b in (a.bams or [])}
    mode = "bam" if bam_of else "calls-fallback"
    if mode == "calls-fallback":
        print("  [warn] no --bams: markers typed from calls; 'not called' is "
              "treated as rCRS only where median depth >= "
              f"{a.min_depth_call}x", file=sys.stderr)

    # --- lineage positions -------------------------------------------------
    shared = Counter(p for v in hom.values() for p in v)
    lineage = set(RECURRENT) | {p for p, n in shared.items()
                                if n >= a.min_samples_lineage}

    # --- assignments -------------------------------------------------------
    rows, marker_rows = [], []
    for s in sorted(hom):
        if mode == "bam" and s in bam_of:
            st = states_from_bam(bam_of[s], a.contig, a.min_baseq, a.min_mapq,
                                 a.min_marker_depth, a.min_major_af)
        else:
            st = states_from_calls(hom[s], med.get(s), a.min_depth_call)
        hg, reason = classify(st)
        typed = sum(1 for v in st.values() if v[0] in ("ref", "alt"))
        mixed = sum(1 for v in st.values() if v[0] == "mixed")
        flags = []
        if mixed >= 2:
            flags.append("possible_mixture")
        if len(hom[s]) == 0:
            flags.append("no_homoplasmic_calls")
        if med.get(s) is not None and med[s] < a.min_depth_call:
            flags.append("below_calling_depth")
        support = [f"m.{p}{MARKERS[p]}" for p in sorted(MARKERS) if st[p][0] == "alt"]
        rows.append([s, hg, reason, f"{typed}/{len(MARKERS)}", mixed,
                     len(hom[s]), "" if med.get(s) is None else f"{med[s]:.0f}",
                     ";".join(flags) or "-", " ".join(support) or "-", mode])
        for p in sorted(MARKERS):
            state, dep, af = st[p]
            block = "N" if p in N_BLOCK else "M" if p in M_BLOCK else "R"
            marker_rows.append([s, p, block, RCRS[p], MARKERS[p], state,
                                "" if dep is None else dep,
                                "" if af is None else f"{af:.4f}"])

    out = os.path.join(a.out_dir, "haplogroups.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "macrohaplogroup", "reason", "markers_typed",
                    "markers_mixed", "n_homoplasmic", "median_depth", "flags",
                    "non_rcrs_markers", "typing_mode"])
        w.writerows(rows)

    mout = os.path.join(a.out_dir, "haplogroup_markers.tsv")
    with open(mout, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "pos", "block", "rcrs_base", "non_rcrs_base",
                    "state", "depth", "major_allele_fraction"])
        w.writerows(marker_rows)

    # --- flag heteroplasmic calls at lineage positions ---------------------
    assign = {r[0]: r[1] for r in rows}
    flagged = os.path.join(a.out_dir, "heteroplasmy_lineage_flags.tsv")
    nflag = ntot = 0
    with open(flagged, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "haplogroup", "pos", "variant", "vaf",
                    "lineage_position", "n_samples_homoplasmic"])
        for s in sorted(het):
            for pos, ref, alt, vaf in het[s]:
                flag = "YES" if pos in lineage else "no"
                nflag += flag == "YES"; ntot += 1
                w.writerow([s, assign.get(s, ""), pos, f"m.{pos}{ref}>{alt}",
                            f"{vaf*100:.2f}", flag, shared.get(pos, 0)])

    # --- report ------------------------------------------------------------
    c = Counter(r[1] for r in rows)
    n = len(rows)
    print(f"  typing mode          : {mode}", file=sys.stderr)
    print(f"  samples              : {n}", file=sys.stderr)
    print("  macrohaplogroups     : "
          + ", ".join(f"{k} {v} ({v/n:.1%})" for k, v in c.most_common()),
          file=sys.stderr)
    assigned = n - c.get("unassigned", 0)
    if c.get("unassigned"):
        print(f"  -> unassigned: " + ", ".join(r[0] for r in rows
              if r[1] == "unassigned"), file=sys.stderr)
    mix = [r[0] for r in rows if "possible_mixture" in r[7]]
    print(f"  possible mixtures    : {len(mix)}"
          + (f" ({', '.join(mix)})" if mix else ""), file=sys.stderr)
    print(f"  lineage positions    : {len(lineage):,}", file=sys.stderr)
    print(f"  heteroplasmic calls at lineage positions: {nflag}/{ntot} "
          f"({nflag/max(ntot,1):.1%})", file=sys.stderr)
    print(f"\n  written: {out}\n           {mout}\n           {flagged}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
