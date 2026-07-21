# conda-forge builds linux-64 against glibc, so this base cannot be musl. The tag must match the
# pixi that wrote pixi.lock, which older releases cannot read.
FROM ghcr.io/prefix-dev/pixi:0.73.0-bookworm-slim

# Nextflow runs ps to collect task metrics.
RUN apt-get update \
    && apt-get install --no-install-recommends --yes procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY pyproject.toml pixi.lock README.md ./
COPY src ./src
RUN pixi install --locked --environment runtime

# Thread count decides the order a matrix product is summed in, so pinning it keeps predictions
# identical to the environment the lockfile describes.
ENV PATH=/src/.pixi/envs/runtime/bin:$PATH \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

WORKDIR /work
CMD ["riborescue", "--help"]
