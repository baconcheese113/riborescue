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

# Records taken from the head of each alignment. STAR writes the transcriptome BAM in read order
# rather than by coordinate, so the head is a fair slice across transcripts, and taking it costs a
# read of the first few hundred megabytes instead of the whole file. A fraction would be tidier and
# would turn a one-minute smoke run into a twenty-minute one.
records=${1:-1500000}
alignments=${2:-results/reads/alignments}
out=${3:-results/smoke/alignments}
pattern=${SMOKE_PATTERN:-*}

command -v samtools >/dev/null || { echo "samtools is not on PATH" >&2; exit 1; }
mkdir -p "$out"

count=0
for bam in "$alignments"/$pattern.Aligned.toTranscriptome.out.bam; do
    [ -e "$bam" ] || { echo "no transcriptome alignments under $alignments" >&2; exit 1; }
    name=$(basename "$bam")
    target="$out/$name"
    if [ -s "$target" ]; then
        echo "have $name"
        count=$((count + 1))
        continue
    fi
    # `head` closes the pipe once it has enough, which reaches samtools as SIGPIPE. Under pipefail
    # that reads as a failed slice, so this one pipeline is exempted rather than the whole script.
    (
        set +o pipefail
        { samtools view -H "$bam"; samtools view "$bam" | head -n "$records"; } \
            | samtools view -b -o "$target.tmp" -
    )
    mv "$target.tmp" "$target"
    echo "sliced $name to $(du -h "$target" | cut -f1)"
    count=$((count + 1))
done
echo "$count libraries under $out, $records records each"
