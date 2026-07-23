#!/usr/bin/env bash
# Refuse to start work that will not fit.
#
# Alignment writes a transcriptome BAM, a genome BAM and a published copy of each, so a dataset
# costs a few hundred gigabytes and the filesystem holding it also holds the container images, the
# STAR indexes and everything already computed. A run that fills the disk half way through leaves
# truncated BAMs that look finished, so the reserve is checked before the run rather than after.

set -euo pipefail

reserve_gb=${1:-200}
path=${2:-.}

free_gb=$(df -BG --output=avail "$path" | tail -1 | tr -dc '0-9')
if [ "$free_gb" -lt "$reserve_gb" ]; then
    echo "only ${free_gb} GB free under $path, below the ${reserve_gb} GB reserve" >&2
    echo "free space or lower the reserve before starting" >&2
    exit 1
fi
echo "${free_gb} GB free, above the ${reserve_gb} GB reserve"
