# Heteroplasmy caller benchmark

Compares the in-house caller with **Mutect2 (mitochondria mode, GATK 4.5.0.0)** and
**mutserve (2.0.0-rc13)**. Both tools are already in `mtdna-pipeline:1.5.0`, so
nothing is installed or downloaded.

## Put it in place

```bash
cd ~/mtdna-cohort47
unzip /mnt/c/Users/<you>/Downloads/benchmark.zip    # creates ./benchmark_kit
cd benchmark_kit
bash run_benchmark.sh check
```

`check` takes under a minute. It confirms the image, the 47 BAMs, the summary
tables, hg38, the BAM header and that mutserve accepts `--contig-name`. Fix
anything it flags before going further.

Default paths (override any of them on the command line, e.g.
`SUMMARY=/mnt/d/.../summary bash run_benchmark.sh check`):

| Variable | Default |
|---|---|
| `FINAL_BAMS` | `/mnt/d/Cancer_microbiome/RESULTS_FINAL/chrM` — post-2b BAMs, the ones the 134 calls came from |
| `SUMMARY` | `/mnt/d/Cancer_microbiome/RESULTS_FINAL/summary` — needs `heteroplasmic_variants.tsv`, `homoplasmic_variants.tsv`, `error_model.json`, `site_background.tsv` |
| `REFS` | `~/ref` (hg38.fa) |
| `BENCH` | `~/mtdna-cohort47/benchmark` — all outputs; must be ext4 |
| `WORK`, `RAW`, `SCRATCH` | as in `cohort.sh`; used only by the `unfiltered` step |
| `PAIRS` | `425:418,414:402,389:393` — M sample diluted into R sample |

## Run

```bash
bash run_benchmark.sh setup        # ~5 min   hg38 .fai/.dict, chrM reference
bash run_benchmark.sh cohort       # ~2-4 h   Mutect2 + mutserve on 47 final BAMs
bash run_benchmark.sh report       # seconds  first results, cohort only
bash run_benchmark.sh mixtures     # ~2-3 h   up to 72 mixtures x 4 callers
bash run_benchmark.sh unfiltered   # ~6-10 h  regenerate Pass-1 BAMs, call both tools
bash run_benchmark.sh report       # full report
```

Every step resumes where it stopped. Run `report` after any step; it uses
whatever outputs exist. Times are estimates for 4 cores.

## What each analysis answers

**1. Cohort (final BAMs).** All three callers see the same NUMT-cleaned reads.
- Do Mutect2 and mutserve also produce a T>A-dominated, low-Ti/Tv call set at low VAF?
- How many of the in-house 57 T>A calls do they make? (`cohort_overlap_with_inhouse.tsv`)
- Does the 5% threshold clean their call sets too? (`cohort_threshold_sweep.tsv`)

**2. NUMT arm (unfiltered Pass-1 BAMs).** Each tool on reads before competitive
realignment vs after. Shows how many NUMT-derived calls each tool makes when
nothing removes NUMT reads upstream. Compare against the in-house ablation
(1,081 → 134). Pass-1 BAMs are rebuilt with the pipeline's own `trim` and
`recruit` rules (`snakemake --notemp`), one sample at a time; intermediates are
deleted as it goes.

**3. Mixtures (ground truth).** Reads from an M sample are mixed into an R sample
at 1, 2, 3, 5, 10 and 20% of coverage, at 200× and 500×, two replicates each.
Every position homoplasmic in only the M sample becomes a true heteroplasmy at a
known fraction. Gives recall, false positives and VAF accuracy per caller,
including the in-house caller with and without the panel of normals.

> **Expect the published in-house caller to lose recall in the mixtures.** Its
> panel of normals sums alt reads from every cohort sample, including those
> homoplasmic for the allele. Mixture truth sites are exactly such lineage
> positions, so the site background there is high. The `nopon` arm isolates
> this. It is a real property of the method (and why only 3/134 calls sit on
> lineage positions); report it rather than tune it away.

## Harmonised call definition

Applied identically to all callers in `compare.py`:
PASS (caller's own filters) · depth ≥ 50 · 1% ≤ VAF < 95% · SNVs only ·
masked m.302–315 and m.16184–16193.

Mutect2 settings follow the GATK mitochondria pipeline (`--mitochondria-mode`,
`--max-mnp-distance 0`, `--max-reads-per-alignment-start 75`), with minimum
MAPQ 30. The shifted-reference step for the control-region breakpoint is not
used, so Mutect2 calls within ~100 bp of m.1/m.16569 may be under-called.
mutserve uses `--level 0.01 --baseQ 20 --mapQ 30`. Read groups, which the
pipeline's BAMs lack and Mutect2 requires, are added to a temporary copy.

## Outputs (`$BENCH/report/`)

| File | Content |
|---|---|
| `benchmark_report.md` | all tables, readable |
| `benchmark_overview.pdf/.png` | Ti/Tv sweep, T>A sweep, mixture recall |
| `cohort_summary.tsv` | calls, Ti/Tv, T>A, C>T, VAF ceilings per caller |
| `cohort_threshold_sweep.tsv` | 1–10% sweep per caller (cf. Table S6) |
| `cohort_overlap_with_inhouse.tsv` | in-house calls also made by each tool |
| `cohort_inhouse_calls_per_caller.tsv` | each of the 134 calls, each tool's VAF |
| `cohort_external_only_calls.tsv` | calls only the external tools make |
| `numt_summary.tsv`, `numt_added_calls_by_500bp.tsv` | NUMT arm |
| `mixture_recall.tsv`, `mixture_false_positives.tsv`, `mixture_vaf_estimates.tsv` | mixtures |

## If something fails

Per-sample logs are in `$BENCH/<arm>/logs/`. The scripts stop on the first
failure and print the last lines of the log. Two things most likely to need a
tweak on first run:

- **mutserve flags**: verified against the mutserve 2 documentation. If
  `check` warns about `--contig-name`, run
  `docker run --rm mtdna-pipeline:1.5.0 java -jar /opt/tools/mutserve/mutserve.jar call --help`
  and adjust `container/run_tools.sh`.
- **Sanity check after `cohort`**: sample 400 m.16093T>C is 90.81% in the
  in-house set, so every caller should report it at roughly 90% in
  `cohort_inhouse_calls_per_caller.tsv`. If mutserve shows ~9% there, its level
  column is the reference allele's and the parser needs flipping.
