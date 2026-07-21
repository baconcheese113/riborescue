# The runtime image Nextflow processes use.
#
# It is built from pixi.lock, so the environment inside the container is the one the tests ran
# against rather than a fresh resolution of the same constraints. The runtime environment shares a
# solve group with the development one, which fixes both to identical package versions.
#
# The base is glibc rather than musl: conda-forge builds linux-64 against glibc, so an Alpine image
# cannot run these packages at all. The base contributes about a tenth of the image either way.
#
# The builder's pixi must be the one that wrote the lockfile; an older one cannot read its format.

FROM ghcr.io/prefix-dev/pixi:0.73.0-bookworm-slim AS build

WORKDIR /src
COPY pyproject.toml pixi.lock README.md ./
COPY src ./src

# --locked refuses to proceed if the manifest and the lockfile disagree, so an image can never be
# built from an environment nobody resolved.
RUN pixi install --locked --environment runtime \
    && pixi shell-hook --environment runtime --shell bash > /activate.sh \
    && echo 'exec "$@"' >> /activate.sh

FROM debian:bookworm-slim AS runtime

# Nextflow shells out to ps to collect the metrics its trace and report files carry, and a task fails
# outright without it.
RUN apt-get update \
    && apt-get install --no-install-recommends --yes procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /src/.pixi/envs/runtime /src/.pixi/envs/runtime
COPY --from=build /activate.sh /activate.sh

# One thread per process, for two reasons that happen to agree. The design matrices are small enough
# that threaded linear algebra spends longer coordinating than computing — a model fit takes 0.53s on
# one thread against 2.39s on sixteen — and the number of threads decides the order a matrix product
# is summed in, which moves the last bits of every prediction. Pinning it makes a container run
# byte-identical to the environment it was tested in.
ENV OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

# Nextflow bind-mounts its work directory, and files written as root are a nuisance to clean up.
RUN useradd --create-home --uid 1000 riborescue
USER riborescue
WORKDIR /work

ENTRYPOINT ["/bin/bash", "/activate.sh"]
CMD ["riborescue", "--help"]
