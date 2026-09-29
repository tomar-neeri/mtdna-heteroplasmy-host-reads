# Mitochondrial heteroplasmy from shotgun metagenomic host reads

Pipeline and results for *Repurposing Shotgun Metagenomic Human Host Reads for
Mitochondrial Heteroplasmy Detection: A Proof-of-Concept Study*.

- **Version:** v1.6 (haplogroup assignment revised; see *Version history*)
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


---

## Principal results

| | |
|---|---|
| Mitochondrial genome recovery | complete in all 47 samples; no uncovered positions |
| Depth | mean 433×, median 226×, range 9–3,029× |
| NUMT contamination | 3.22% of recruited reads relocate to nuclear loci; removing them eliminated 96% of heteroplasmy calls in validation |
| Microbial cross-mapping | <0.05% of candidate read pairs; negligible after competitive human alignment |
| Error model | pooled, substitution-specific beta-binomial; Q37.0 over 2,083,234 background observations |
| Heteroplasmy calls | 134 across 47 samples; 16 clear the 5% allele-fraction limit |
| Detection limit | statistical 1.9%; **empirical 5%** |

The central finding is that the statistical model alone does not give a
trustworthy detection limit. Below 5% allele fraction the call set is dominated
by a T>A error signature whose mechanism is unresolved: it is not sequence
context, not homopolymer-driven, and not positional. It is controlled by an
allele-fraction threshold rather than by the model. Heteroplasmy calls from
these data should not be interpreted below approximately 5%.

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
| `cohort.sh`, `manage.sh`, `run_container.sh` | run and maintenance wrappers |
| `rerun_haplogroups.sh` | applies the v1.6 haplogroup fix to existing outputs |
| `CONTAINERIZATION.md` | container build and usage notes |

### results/summary/

| File | Contents |
|---|---|
| `error_model.json` | fitted error rates for all 12 substitution classes |
| `site_background.tsv` | per-position background from the panel of normals |
| `variant_recurrence.tsv` | positions called in more than one sample |
| `region_distribution.tsv` | calls by mitochondrial region |
| `cohort_overview.tsv` | cohort-level summary counts |
| `versions.txt` | tool versions and reference checksums for this run |

Per-participant tables (per-sample haplogroups, marker bases, and annotated
variant calls) are **not** included in this deposit. A mitochondrial haplotype
is identifying: it is heritable, shared with maternal relatives, and matchable
against other datasets. Sequence data (FASTQ, BAM) are likewise not deposited
here and require controlled access. Requests for individual-level data should
be directed to the corresponding author and are subject to the original
consent and institutional approval.

---

## Reproducing the analysis

The pipeline runs offline in a pinned Docker container.

```bash
# build
docker build -t mtdna-pipeline:1.6.0 -f code/Dockerfile .

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
[Authors]. Repurposing shotgun metagenomic human host reads for mitochondrial
heteroplasmy detection: a proof-of-concept study. [Journal, year].
[Authors]. Pipeline and results (v1.6). Zenodo. https://doi.org/[DOI]
```

---

## Known limitations

- **No external validation.** Public metagenomes are routinely host-scrubbed
  before deposit (NCBI's Human Read Removal Tool), so the host reads this
  approach depends on are generally unavailable in open archives.
- **No orthogonal confirmation** of individual calls, and no technical
  replicate libraries.
- **No benchmark** against mutserve or Mutect2's mitochondrial mode.
- **The T>A error mechanism is unresolved.** It is controlled by a threshold,
  not explained.
- **Single tissue, population and platform.** The thresholds reported here
  should be re-derived rather than transferred to other settings.
