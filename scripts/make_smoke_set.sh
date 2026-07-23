#!/usr/bin/env bash
# A small deterministic slice of every aligned library, for developing the workflow against.
#
# The validation workflow was built against whole libraries, which meant every wiring mistake cost
# half an hour to find. This subsamples each transcriptome alignment to a fixed fraction with a fixed
# seed, so the same reads come out every time and the whole chain runs in about a minute.
#
# Its numbers are mechanically useful and scientifically meaningless: a thousandth of a library
# cannot support a readthrough claim, and the calibration thresholds that a real dataset must clear
# are deliberately left where they are, so a smoke run is expected to report itself inconclusive.

set -euo pipefail

fraction=${1:-0.002}
alignments=${2:-results/reads/alignments}
out=${3:-results/smoke/alignments}
seed=1

command -v samtools >/dev/null || { echo "samtools is not on PATH" >&2; exit 1; }
mkdir -p "$out"

count=0
for bam in "$alignments"/*.Aligned.toTranscriptome.out.bam; do
    [ -e "$bam" ] || { echo "no transcriptome alignments under $alignments" >&2; exit 1; }
    name=$(basename "$bam")
    target="$out/$name"
    if [ -s "$target" ]; then
        echo "have $name"
        count=$((count + 1))
        continue
    fi
    # -s takes seed.fraction as one number, so the same slice is reproduced on every machine.
    samtools view -@ 2 -b -s "${seed}${fraction#0}" "$bam" >"$target.tmp"
    mv "$target.tmp" "$target"
    echo "sliced $name to $(du -h "$target" | cut -f1)"
    count=$((count + 1))
done
echo "$count libraries under $out at fraction $fraction"
