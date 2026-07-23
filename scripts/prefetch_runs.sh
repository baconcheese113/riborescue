#!/usr/bin/env bash
# Fetch the FASTQ a samplesheet names, several streams at once.
#
# The archive throttles a single connection to about 1.5 MB/s regardless of the link, so fetching
# thirty libraries one after another takes hours of doing nothing. Streams are independent, so
# several at once multiply the rate.
#
# This only puts bytes on disk. `riborescue stage-runs` remains what decides a file is the declared
# one: a download lands beside the target and is moved into place only once its checksum matches, so
# an interrupted fetch leaves nothing that could be mistaken for a complete file.

set -euo pipefail

samplesheet=${1:-pipeline/assets/riboseq_samples.tsv}
streams=${2:-6}
dest=${RIBORESCUE_DATA:-data}/fastq
mkdir -p "$dest"

fetch() {
    local url=$1 md5=$2 name
    name=$(basename "$url")
    local target="$dest/$name" partial="$dest/.$name.partial"

    if [ -f "$target" ] && [ "$(md5sum "$target" | cut -c1-32)" = "$md5" ]; then
        echo "have $name"
        return 0
    fi
    # Resumed rather than restarted: a partial file from an interrupted run is most of the work.
    curl -sSfL -C - "$url" -o "$partial"
    if [ "$(md5sum "$partial" | cut -c1-32)" != "$md5" ]; then
        echo "checksum mismatch for $name" >&2
        rm -f "$partial"
        return 1
    fi
    mv "$partial" "$target"
    echo "fetched $name"
}
export -f fetch
export dest

awk -F'\t' '
    NR==1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    { print $h["fastq_1_url"], $h["fastq_1_md5"]
      if ($h["fastq_2_url"] != "") print $h["fastq_2_url"], $h["fastq_2_md5"] }
' "$samplesheet" | xargs -P "$streams" -n 2 bash -c 'fetch "$0" "$1"'
