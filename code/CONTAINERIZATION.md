# Containerized pipeline — build, run, reproduce

Turns the interactive workflow into a single reproducible command. Three bugs
found during development are now structurally prevented rather than documented.

---

## What is in the box

| File | Role |
|---|---|
| `Dockerfile` | Pinned image: tools, Java binaries, rCRS baked in |
| `envs/mtdna.yaml` | Exact conda versions — no version ranges anywhere |
| `Snakefile` | The three-pass DAG with validation gates |
| `config/config.yaml` | Every parameter. Deposit this with results. |

---

## Build

```bash
mkdir -p refs && cp /path/to/rCRS.fasta refs/
docker build -t mtdna-pipeline:1.0.0 .
```

The rCRS is baked in deliberately. It is 16.6 kb, immutable, and every
coordinate you report depends on it — a mismatched reference silently shifts
every position. hg38 stays outside the image because it is 3 GB and changes
with assembly patches.

For HPC without Docker:

```bash
apptainer build mtdna.sif docker-daemon://mtdna-pipeline:1.0.0
```

---

## Prepare the human reference (once, outside the container)

```bash
mkdir -p refs && cd refs
aria2c -x16 https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
docker run --rm -v $PWD:/refs mtdna-pipeline:1.0.0 bwa index /refs/hg38.fa
md5sum hg38.fa > hg38.fa.md5      # record for the methods section
```

`bwa index` takes ~2 hours and ~4.5 GB RAM. One time.

---

## Run

```bash
docker run --rm \
  -v /path/to/fastqs:/data/raw:ro \
  -v /path/to/refs:/refs:ro \
  -v /fast/scratch:/scratch \
  -v $PWD:/work -w /work \
  --memory 12g --cpus 4 \
  mtdna-pipeline:1.0.0 \
  snakemake -c4 --configfile config/config.yaml
```

**`/scratch` must be a fast local filesystem.** Alignment does heavy random
I/O; a network mount or a WSL2 `/mnt/d` 9p mount is 10–20× slower and will
dominate runtime.

Useful invocations:

```bash
snakemake -c4 -n                    # dry run: what would execute
snakemake --dag | dot -Tsvg > dag.svg
snakemake -c4 --until verify_human  # stop after Pass 2a
snakemake -c4 --rerun-incomplete    # resume after interruption
```

---

## Bugs this prevents

Each cost real debugging time. All three now fail loudly rather than silently.

**1 · R1/R2 desynchronisation.** `samtools fastq` cannot pair reads from a
coordinate-sorted BAM. Four reads drifted out of sync, `bwa mem` emitted
`[W::bseq_read] the 1st file has fewer sequences`, and produced an empty BAM
that looked like a legitimate result. The `candidate_fastq` rule name-sorts
first, routes singletons to `/dev/null`, and **asserts the counts match**.

**2 · Silent OOM.** `bwa mem` against hg38 needs ~6 GB. Run alongside
`bwa index`, it was killed, and stderr went to `/dev/null`. Every rule now
writes stderr to `logs/`, and `verify_human` **fails if the output has zero
reads**.

**3 · Contig-name mismatch.** `chrM` / `MT` / `NC_012920.1` all appear in the
wild. Now a single `chrm_contig` config key, used everywhere.

---

## Critical: the working directory must be on a native Linux filesystem

Snakemake writes lock and state files into `.snakemake/` continuously to track
running jobs. On a Windows mount under WSL2 (`/mnt/d`, 9p/drvfs) those writes
fail with `Operation not permitted`, Snakemake loses track of the job, and
**reports failure while the tool is still running successfully**.

Observed symptom: `trim` reported an error at 18 s; the fastp log showed it ran
for 87 s and completed normally, writing all outputs.

```bash
mkdir -p ~/mtdna-work          # ext4
```

Mount read-only inputs and references from anywhere — sequential reads over 9p
are fine. Only `/work` has this constraint.

Use `run_container.sh`, which checks the filesystem type before launching and
refuses to start if the working directory is on 9p.

---

## Reproducibility guarantees

| Layer | Mechanism |
|---|---|
| Tool versions | Exact pins in `envs/mtdna.yaml`, no ranges |
| Java tools | Release tags as build args |
| rCRS | Baked into image, length-asserted at build |
| hg38 | md5 recorded in `versions.txt` at runtime |
| Parameters | `config.yaml`, copied verbatim into `versions.txt` |
| Provenance | `results/summary/versions.txt` per run |

Pin the image by digest in any publication:

```bash
docker inspect --format='{{index .RepoDigests 0}}' mtdna-pipeline:1.0.0
```

That digest, not the tag, is what makes the run reproducible — tags can be
re-pushed.

---

## Outputs

```
results/
├── chrM/       <sample>.chrM.bam        final, verified
├── qc/         fastp reports, coverage, numt_destinations, kraken
├── logs/       per-rule stderr
└── summary/
    ├── coverage_metrics.tsv
    ├── filtering_stats.tsv     ← NUMT removal per sample
    └── versions.txt            ← provenance
```

`filtering_stats.tsv` is the one to look at first. It reports reads retained,
reads relocated, and the top nuclear destination per sample — the Pass 2a
result that drove the headline finding.

---

## Pass 2b is off by default

On the validation sample, Pass 2a alone removed 96% of spurious heteroplasmy
calls (174 → 7) by relocating 2.7% of reads, mostly to chr1. Whether Kraken2
adds anything is an open question, and the config defaults to `false` so that
enabling it is a deliberate, reportable choice.

Test it properly: run a high-bacterial sample both ways and compare calls. If
2b changes nothing there either, the negative result is worth reporting —
it means alignment alone suffices for this data type.

```yaml
run_kraken: true
kraken_db: "/refs/k2_standard_08gb"
```

---

## Publishing

1. Push the image to a registry and record the digest
2. Tag the repo to match the image version
3. Archive the tagged release on Zenodo for a DOI
4. Deposit `config.yaml` and `versions.txt` as supplementary files
5. Cite the hg38 URL and md5, not just "hg38"

A reader with the digest, the config, and the raw FASTQs reproduces your
numbers exactly. That is the bar.
