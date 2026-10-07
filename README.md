# Mitochondrial heteroplasmy from shotgun metagenomic host reads

Pipeline and results for *Repurposing Shotgun Metagenomic Human Host Reads for
Mitochondrial Heteroplasmy Detection*.

- **Version:** v1.7 (NUMT ablation added; see *Version history*)
- **Contact:** Siddharth Singh Tomar
- **Licence:** MIT for `code/`, CC BY 4.0 for `results/` and `figures/`
- **Related manuscript:** [add DOI or "in preparation" on release]

---

## What this is

Shotgun metagenomic sequencing of human samples produces a large fraction of
host reads that microbiome pipelines discard. Because mitochondrial DNA is
present at hundreds to thousands of copies per cell, those discarded reads
carry enough mitochondrial sequence to reconstruct the genome and, within
limits, to call heteroplasmy.

This deposit contains the pipeline that does so and the analytical limits it
establishes. The study is methodological. It makes no clinical claim: the
cohort contains no cancer cases, there are no matched normal samples, and
nasopharyngeal swabs are not the relevant tissue for the cancers that motivated
the wider project.

**Cohort:** 47 nasopharyngeal swabs from unselected adults, Nagpur, India.
Paired-end Illumina sequencing, processed through CZ ID.

---

## Principal results

| | |
|---|---|
| Mitochondrial genome recovery | complete in all 47 samples; no uncovered positions |
| Depth | mean 433×, median 226×, range 9–3,029× |
| NUMT contamination | 3.22% of recruited reads relocate to nuclear loci; omitting competitive realignment yields 1,081 calls against 134 (87.6% of calls removed) |
| Microbial cross-mapping | <0.05% of candidate read pairs; negligible after competitive human alignment |
| Error model | pooled, substitution-specific beta-binomial; Q37.0 over 2,083,234 background observations |
| Heteroplasmy calls | 134 across 47 samples; 16 clear the 5% allele-fraction limit |
| Detection limit | statistical 1.9%; **empirical 5%** |

Two distinct artifacts limit heteroplasmy detection in these data, and they
have different signatures, so neither control substitutes for the other.

**Nuclear mitochondrial insertions (NUMTs)** are removed by competitive
realignment against the complete human genome. They are *not* identifiable
from the substitution spectrum: without realignment the call set has a
Ti/Tv of 3.11, close to gnomAD's 3.51, so a spectrum-plausibility check passes
a call set that is 87.6% artifact. They are identifiable from where they
fall — 88% of the additional calls lie within 8 kb of the 16.6 kb genome.

**A dispersed, T>A-dominated error process** survives realignment and is what
the substitution spectrum does catch. It is unresolved mechanistically: not
sequence context, not homopolymer-driven, not positional. It is controlled by
an allele-fraction threshold rather than by the error model.

Heteroplasmy calls from these data should not be interpreted below
approximately 5% allele fraction.

---

## NUMT ablation

The pipeline was rerun on all 47 samples with competitive realignment to
GRCh38 removed, refitting the error model and the panel of normals on the
unfiltered alignments. Cohort-level outputs are in
`results/summary/numt_ablation_comparison.tsv`.

| | no realignment | published | |
|---|---|---|---|
| heteroplasmic calls | 1,081 | 134 | 8.07× |
| calls ≥5% allele fraction | 79 | 16 | 4.9× |
| Ti/Tv | 3.11 | 0.43 | gnomAD 3.51 |
| C>T calls | 332 | 1 | |
| T>A calls | 99 (9.2%) | 57 (42.5%) | |
| maximum T>A allele fraction | 5.59% | 4.72% | |
| calls on gnomAD NUMT false-positive sites | 0 | 0 | |

By region, the additional calls are overwhelmingly ribosomal: 596 calls in
MT-RNR1 and MT-RNR2 against 8 in the published set, a 74-fold difference,
while the control region changes only from 23 to 28 calls.

Three results from the ablation are worth stating separately.

**NUMT false positives reach high allele fraction.** m.3106C>A appears at
39.62% and 31.37% in two samples without realignment and in neither with it.
A 39.6% call would have been reported as confident and biologically real. This
runs against the common framing of NUMT contamination as a low-heteroplasmy
problem.

**The additional calls are spatially banded, and the bands separate by
spectrum.** The window m.1,000–1,499 alone holds 295 calls, 27.3% of the
total. High-Ti/Tv bands (m.500–3,499 at 5.59; m.5,500–5,999 at 9.09;
m.9,500–9,999 at 6.38) are consistent with transition-diverged nuclear
insertions. Low-Ti/Tv windows (m.12,000–12,499 at 0.21; m.13,000–13,499 at
0.85) are the T>A process. The two artifacts separate in one table.

**The 5% threshold is robust to the filter.** Above 5% in the unfiltered arm,
only 1 of 79 calls is T>A and Ti/Tv rises to 4.64. The threshold cleans up even
unfiltered data.

Two caveats. Refitting on unfiltered alignments raises the fitted background
from Q37.0 to approximately Q36.3, because NUMT mismatches present as
low-allele-fraction noise; the refitted arm therefore *understates* the NUMT
contribution, since a higher background tightens the significance threshold.
And the T>A maximum of 4.72% is a property of the filtered data at this depth
and sample size, not an absolute ceiling of the error process — it reaches
5.59% unfiltered.

---

## Contents

```
code/                   pipeline and analysis scripts
results/summary/        cohort-level result tables
figures/                manuscript figures (PDF and PNG)
supplementary/          supplementary tables
MANIFEST.tsv            path, size and SHA-256 for every file
```

### code/

| File | Purpose |
|---|---|
| `Snakefile` | workflow definition, 16 rules |
| `config.yaml` | thresholds and reference paths |
| `Dockerfile`, `mtdna.yaml` | pinned container environment |
| `pooled_error_model.py` | fits the substitution-specific beta-binomial error model |
| `site_background.py` | builds the leave-one-out panel of normals |
| `call_heteroplasmy.py` | variant calling against the applicable null |
| `haplogroup.py` | macrohaplogroup assignment from 11 diagnostic markers |
| `lineage_matrix.py` | cross-tabulates calls against the cohort homoplasmic landscape |
| `collect_variants.py` | gene, codon and amino-acid annotation; gnomAD status |
| `kraken_filter.py` | taxonomic read filtering |
| `plot_haplogroups.py` | draws Figure 3c–e |
| `numt_ablation.py` | reruns the cohort without competitive realignment and compares arms |
| `cohort.sh`, `manage.sh`, `run_container.sh` | run and maintenance wrappers |
| `rerun_haplogroups.sh` | applies the v1.6 haplogroup fix to existing outputs |
| `CONTAINERIZATION.md` | container build and usage notes |

### results/summary/

| File | Contents |
|---|---|
| `error_model.json` | fitted error rates for all 12 substitution classes |
| `site_background.tsv` | per-position background from the panel of normals |
| `region_distribution.tsv` | calls by mitochondrial region |
| `numt_ablation_comparison.tsv` | cohort, region and 500 bp window comparison of the two arms |
| `refit/error_model.json` | error model refitted on unfiltered alignments |
| `versions.txt` | tool versions and reference checksums for this run |

Per-participant tables are **not** included in this deposit. That covers
per-sample haplogroups, marker bases, annotated variant calls, coverage and
call-count summaries, the per-sample ablation arms, and `variant_recurrence.tsv`,
which lists the samples carrying each position. A mitochondrial haplotype is
identifying: it is heritable, shared with maternal relatives, and matchable
against other datasets. Sequence data (FASTQ, BAM) are likewise not deposited
and require controlled access. Requests for individual-level data should be
directed to the corresponding author and are subject to the original consent
and institutional approval.

---

## Reproducing the analysis

The pipeline runs offline in a pinned Docker container.

```bash
# build
docker build -t mtdna-pipeline:1.7.0 -f code/Dockerfile .

# run (see config.yaml for the paths it expects)
snakemake -c4 --configfile config/config.yaml
```

Approximately 6 hours for 47 samples on 4 cores with a 12 GB cap.

**References required** (not redistributed here):

- rCRS / chrM (NC_012920.1)
- GRCh38 with BWA index
- gnomAD v3.1 mitochondrial variant table
- Kraken2 standard database (optional; the microbial step is off by default)

### Three failure modes worth knowing

These cost time to diagnose and are documented in the code:

1. **Kraken2 "unclassified" filtering.** The standard database contains the
   human genome, so genuine mitochondrial reads classify as *Homo sapiens*.
   Retaining only unclassified reads discards 92% of real data. Filter by
   taxonomy instead.
2. **Beta-binomial overflow.** An unbounded concentration parameter drives
   `betabinom.sf` to return exactly 0 everywhere, silently disabling the
   significance test. The concentration is bounded at 10⁶.
3. **R1/R2 desynchronisation.** `samtools fastq` cannot pair reads from a
   coordinate-sorted BAM. Name-sort first, or `bwa mem` emits no alignments
   while exiting successfully.

---

## Version history

**v1.7** — NUMT ablation.

The cohort was rerun without competitive realignment to quantify the NUMT
contribution directly rather than by inference. Three statements from earlier
versions changed as a result. The figure for calls removed by realignment is
**87.6% cohort-wide** (1,081 → 134), superseding the 96% previously quoted from
a smaller validation set. The claim that no call falls on a gnomAD NUMT
false-positive site is true but carries no evidential weight, since no call
falls there in either arm; the spatial banding is the real independent
evidence. And the T>A allele-fraction ceiling is filter- and cohort-specific
rather than an absolute property of the error process.

No pipeline code changed in this version; the ablation is an additional arm.

**v1.6** — haplogroup assignment corrected.

v1.5 treated 8701G, 9540C, 10398G, 10873C and 15301A as haplogroup N markers.
These are the ancestral L3/M states: the rCRS lies within N and R and carries
the N-derived bases, so every M genome differs from the reference at these
positions (Rajkumar et al. 2005, *BMC Evol Biol* 5:26). The test therefore
detected samples *outside* N. v1.5 also assigned R by default when no marker
was found, so a sample with insufficient coverage was labelled R rather than
reported as untyped, and it could not detect non-R N lineages.

v1.6 reads the base at each of 11 markers directly from the alignment (≥10
reads, major allele ≥0.80), adds m.12705 and m.16223 to separate R from other
N lineages, and reports `unassigned` when markers are uncovered or conflicting.
Per-marker states, depths and allele fractions are written to
`haplogroup_markers.tsv`.

Effect on the results: two samples changed group (one M → unassigned at 16×
mean depth, one N → N outside R). Assignments are M 22, R 23, N outside R 1,
unassigned 1. No other pipeline output changed; the lineage cross-tabulation
is identical.

**v1.5** — initial complete run of all 47 samples.

---

## Citation

If you use this pipeline, please cite the manuscript and this deposit.

```
[SS Tomar, Krishna Khairnar]. Repurposing shotgun metagenomic human host reads for mitochondrial
heteroplasmy detection. [2026].
[SS Tomar, Krishna Khairnar]. Pipeline and results (v1.7). Zenodo. https://doi.org/[DOI]
```

---

## Known limitations

- **No external validation in this study.** Replication requires data in which
  host reads survive deposition. Host-scrubbing of public metagenomes is
  neither universal nor mitochondria-aware, so such data do exist in open
  archives; identifying a suitable set is outstanding work rather than an
  impossibility.
- **No orthogonal confirmation** of individual calls, and no technical
  replicate libraries.
- **No benchmark** against mutserve or Mutect2's mitochondrial mode. Reported
  inter-caller concordance for heteroplasmic variants is low, so such a
  comparison would be diagnostic rather than validating.
- **The T>A error mechanism is unresolved.** It is controlled by a threshold,
  not explained.
- **The ablation's upper bound is pending.** The refitted arm understates the
  NUMT contribution; the fixed-model arm bounds it from the other side and is
  not yet reflected here.
- **Single tissue, population and platform.** The thresholds reported here
  should be re-derived rather than transferred to other settings.
