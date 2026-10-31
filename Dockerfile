FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337 AS genomes

# GENCODE basic annotations, the files SnapATAC2 2.10.0's genome objects (snap.genome.GRCh38 /
# GRCm39) download lazily; fetched here with the sha256 SnapATAC2 pins for them, so a run needs
# no network. Chromosome sizes are SnapATAC2's own tables (primary chromosomes).
ADD --checksum=sha256:b82a655bdb736ca0e463a8f5d00242bedf10fa88ce9d651a017f135c7c4e9285 \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.basic.annotation.gff3.gz \
    /opt/genomes/GRCh38/annotation.gff3.gz
ADD --checksum=sha256:6f433e2676e26569a678ce78b37e94a64ddd50a09479e433ad6f75e37dc82e48 \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.basic.annotation.gff3.gz \
    /opt/genomes/GRCm39/annotation.gff3.gz

FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /processor

# Pinned Python dependencies, binary wheels only, installed without resolving (requirements.txt is
# the complete set; see its header for what is deliberately left out).
COPY requirements.txt /processor/requirements.txt
RUN pip install --only-binary=:all: --no-deps -r /processor/requirements.txt \
    && find /usr/local/lib/python3.12/site-packages -name '*.pyc' -delete \
    && find /usr/local/lib/python3.12/site-packages -depth -type d \( -name tests -o -name test \) -exec rm -rf {} + \
    && rm -rf /usr/local/lib/python3.12/site-packages/pyarrow/include

COPY --from=genomes /opt/genomes /opt/genomes
RUN python -c "import snapatac2 as s; [open(f'/opt/genomes/{n}/chrom_sizes.tsv', 'w').write(''.join(f'{c}\t{l}\n' for c, l in getattr(s.genome, n).chrom_sizes.items())) for n in ('GRCh38', 'GRCm39')]" \
    && chmod -R a+rX /opt/genomes

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    GENOMES_DIR=/opt/genomes \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba \
    HOME=/tmp

# Matplotlib's font cache built now; the cache folders it leaves in /tmp are made writable for uid 1000
# (the platform does not remount /tmp).
RUN python -c "import matplotlib.pyplot, snapatac2" \
    && mkdir -p /tmp/matplotlib /tmp/numba && chmod -R a+rwX /tmp/matplotlib /tmp/numba

# Runs as uid 1000 on the platform; nothing here needs root.
USER 1000

ENTRYPOINT ["python3", "-m", "processor.entry"]
