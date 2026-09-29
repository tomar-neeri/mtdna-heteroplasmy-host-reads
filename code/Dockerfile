# syntax=docker/dockerfile:1
# ---------------------------------------------------------------------------
# mtDNA heteroplasmy recovery from shotgun metagenomes
#
#   docker build -t mtdna-pipeline:1.0.0 .
#   docker run --rm -v $PWD:/work -w /work mtdna-pipeline:1.0.0 snakemake -c4
#
# Every version is pinned. Rebuilding this file in five years produces the
# same image, which is the entire point.
# ---------------------------------------------------------------------------
FROM mambaorg/micromamba:1.5.8-jammy

LABEL org.opencontainers.image.title="mtdna-pipeline"
LABEL org.opencontainers.image.description="Mitochondrial heteroplasmy recovery from host reads in shotgun metagenomic data"
LABEL org.opencontainers.image.licenses="MIT"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
      procps unzip curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# --- conda environment (pinned) --------------------------------------------
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/mtdna.yaml /tmp/mtdna.yaml
# patient downloads: a stalled mirror should retry, not abort the build
ENV MAMBA_NO_LOW_SPEED_LIMIT=1 MAMBA_REMOTE_MAX_RETRIES=10
RUN micromamba install -y -n base -f /tmp/mtdna.yaml \
    && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH

# --- Java tools (pinned release tags) --------------------------------------
# Not in bioconda at usable versions; fetched from GitHub releases.
ARG MUTSERVE_VER=2.0.0-rc13
ARG HAPLOGREP_VER=3.2.1
ARG HAPLOCHECK_VER=1.3.3

USER root
RUN mkdir -p /opt/tools && cd /opt/tools \
 && curl -sSfL "https://github.com/seppinho/mutserve/releases/download/v${MUTSERVE_VER}/mutserve.zip" -o m.zip \
 && unzip -oq m.zip -d mutserve && rm m.zip \
 && curl -sSfL "https://github.com/genepi/haplogrep3/releases/download/v${HAPLOGREP_VER}/haplogrep3-${HAPLOGREP_VER}-linux.zip" -o h.zip \
 && unzip -oq h.zip -d haplogrep3 && rm h.zip \
 && curl -sSfL "https://github.com/genepi/haplocheck/releases/download/v${HAPLOCHECK_VER}/haplocheck.zip" -o c.zip \
 && unzip -oq c.zip -d haplocheck && rm c.zip \
 && chmod -R a+rX /opt/tools && find /opt/tools -type f -name "haplo*" -o -name "mutserve" | xargs -r chmod a+x

ENV MUTSERVE_JAR=/opt/tools/mutserve/mutserve.jar
ENV HAPLOGREP=/opt/tools/haplogrep3/haplogrep3
ENV HAPLOCHECK=/opt/tools/haplocheck/haplocheck

# --- reference: rCRS shipped in-image, hg38 mounted at runtime -------------
# rCRS is 16.6 kb and immutable, so baking it in guarantees coordinates match.
# hg38 is 3 GB and stays outside the image.
COPY --chown=root:root refs/rCRS.fasta /opt/refs/rCRS.fasta
RUN samtools faidx /opt/refs/rCRS.fasta \
 && bwa index /opt/refs/rCRS.fasta \
 && python3 -c "\
from pathlib import Path; \
s=''.join(Path('/opt/refs/rCRS.fasta').read_text().splitlines()[1:]); \
assert len(s)==16569, f'rCRS length {len(s)}'; \
sh=s[8000:]+s[:8000]; \
Path('/opt/refs/rCRS_shifted.fasta').write_text('>chrM_shifted\n'+'\n'.join(sh[i:i+60] for i in range(0,len(sh),60))+'\n')" \
 && samtools faidx /opt/refs/rCRS_shifted.fasta \
 && bwa index /opt/refs/rCRS_shifted.fasta

ENV RCRS=/opt/refs/rCRS.fasta
ENV RCRS_SHIFTED=/opt/refs/rCRS_shifted.fasta

# analysis scripts (pooled error model, calling, kraken filter)
COPY --chown=root:root scripts/ /opt/scripts/
RUN chmod -R a+rX /opt/scripts \
 && test -f /opt/scripts/collect_variants.py \
 && test -f /opt/scripts/site_background.py \
 && test -f /opt/scripts/haplogroup.py \
 && test -f /opt/scripts/lineage_matrix.py

# provenance — MUST run as root, before the USER switch, or /opt is unwritable
RUN micromamba list -n base > /opt/package_versions.txt && chmod a+r /opt/package_versions.txt

USER $MAMBA_USER
WORKDIR /work

# version LAST: bumping it then rebuilds only this layer, not the conda env
LABEL org.opencontainers.image.version="1.5.1"

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["snakemake", "--help"]
