# The runtime image Nextflow processes use. It carries the riborescue command and nothing else: the
# Pixi environment is for development, and a process that needs a tool the package does not install
# should say so rather than inherit one from a shared image.
FROM python:3.12-slim-bookworm

# cyvcf2 builds against htslib; the compiler is not kept in the runtime layer.
RUN apt-get update \
    && apt-get install --no-install-recommends --yes build-essential zlib1g-dev libbz2-dev \
        liblzma-dev libcurl4-openssl-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/riborescue
COPY pyproject.toml README.md ./
COPY src ./src

RUN pip install --no-cache-dir . \
    && apt-get purge --yes build-essential \
    && apt-get autoremove --yes

# A non-root user, because Nextflow bind-mounts the work directory and files written as root are a
# nuisance to clean up afterwards.
RUN useradd --create-home --uid 1000 riborescue
USER riborescue
WORKDIR /work

ENTRYPOINT []
CMD ["riborescue", "--help"]
