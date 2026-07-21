# The runtime image Nextflow processes use.
#
# It is built from pixi.lock, so the environment inside the container is the one the tests ran
# against rather than a fresh resolution of the same constraints. The runtime environment shares a
# solve group with the development one, which fixes both to identical package versions. The pixi that
# builds it must be the one that wrote the lockfile; an older release cannot read its format.
#
# The base is glibc rather than musl: conda-forge builds linux-64 against glibc, so the environment's
# own Python cannot start on Alpine.
FROM ghcr.io/prefix-dev/pixi:0.73.0-bookworm-slim

# Nextflow shells out to ps to collect the metrics its trace and report files carry, and a task fails
# outright without it.
RUN apt-get update \
    && apt-get install --no-install-recommends --yes procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY pyproject.toml pixi.lock README.md ./
COPY src ./src

# --locked refuses to proceed if the manifest and the lockfile disagree, so an image can never be
# built from an environment nobody resolved.
RUN pixi install --locked --environment runtime

# The environment needs no activation beyond its bin directory; nothing in it reads the variables an
# activation script would set.
#
# One thread per process, for two reasons that agree. These design matrices are small enough that
# threaded linear algebra spends longer coordinating than computing — a model fit takes 0.53s on one
# thread against 2.39s on sixteen — and the thread count decides the order a matrix product is summed
# in, which moves the last bits of every prediction. Pinning it makes a container run byte-identical
# to the environment it was tested in.
ENV PATH=/src/.pixi/envs/runtime/bin:$PATH \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

WORKDIR /work
CMD ["riborescue", "--help"]
