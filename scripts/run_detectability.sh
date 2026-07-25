#!/usr/bin/env bash
# The detectability arm of ADR-0019: how large an effect a design could have resolved.
#
# Only for a dataset whose libraries fail calibration on depth. Both length windows are run, and
# both are reported, in the order this script fixes rather than one chosen after the numbers exist.
# The contrasts resting on a library that failed on validity are put to the command too, so that its
# refusal and the reason behind it land on the log beside the arm that ran; nothing is computed for
# them, which is why they are there.
#
#   run_detectability.sh <dataset> <treated> <control> [refused-arm:its-control ...]

set -euo pipefail

dataset=${1:-gse179274_fibroblast}
treated=${2:-suptrna_tyr}
control=${3:-egfp}
[ "$#" -ge 3 ] && shift 3 || shift "$#"

samplesheet=${SAMPLESHEET:-results/staged_runs.tsv}
counts=${COUNTS:-results/psite/$dataset/readthrough_counts.tsv}
gtf=${GTF:-data/gencode/gencode.v50.primary_assembly.annotation.gtf.gz}
manifest=${MANIFEST:-results/psite/$dataset/calibration.json}
reference=${REFERENCE:-results/readthrough/gse144140/g418_vs_dmso.unpaired}
outdir=${DETECTABILITY_OUTDIR:-results/exploratory/detectability/$dataset}

[ -f "$counts" ] || { echo "no counts at $counts — calibrate $dataset first" >&2; exit 1; }
[ -f "$manifest" ] || { echo "no calibration manifest at $manifest" >&2; exit 1; }
mkdir -p "$outdir"

# run <label> <reference> <treated> <control> [extra riborescue arguments...]
run() {
    riborescue detectability "$counts" --gtf "$gtf" --samplesheet "$samplesheet" \
        --dataset "$dataset" --manifest "$manifest" --reference "$2" \
        --treated "$3" --control "$4" "${@:5}" \
        --out "$outdir/${3}_vs_${4}.${1}.tsv"
}

echo "=== $treated against $control, selected set ==="
run selected "$reference.tsv" "$treated" "$control"

echo
echo "=== $treated against $control, published window ==="
run published "$reference.published.tsv" "$treated" "$control" --published-lengths 28 35

for pair in "$@"; do
    echo
    echo "=== ${pair%%:*} against ${pair##*:}, guard check — expected to be refused ==="
    run selected "$reference.tsv" "${pair%%:*}" "${pair##*:}" || true
done
