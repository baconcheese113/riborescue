#!/usr/bin/env bash
# Every treated arm of one dataset against its control arm, by both estimators.
#
# Which estimator is primary is fixed per dataset in ADR-0010 and not chosen from the answers, so
# both are run every time and both are written. A contrast is never selected afterwards either:
# the arms come from the samplesheet, so an arm that disappoints still appears in the results.

set -euo pipefail

usage() {
    echo "usage: $0 <dataset> <control-arm> [samplesheet] [counts] [gtf]" >&2
    exit 2
}

dataset=${1:-}
control=${2:-}
samplesheet=${3:-results/staged_runs.tsv}
counts=${4:-results/psite/$dataset/readthrough_counts.tsv}
gtf=${5:-data/gencode/gencode.v50.primary_assembly.annotation.gtf.gz}
[ -n "$dataset" ] && [ -n "$control" ] || usage
[ -f "$counts" ] || { echo "no counts at $counts — calibrate $dataset first" >&2; exit 1; }

outdir="results/readthrough/$dataset"
mkdir -p "$outdir"

arms=$(awk -F'\t' -v want="$dataset" -v control="$control" '
    NR==1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    $h["assay"] != "riboseq" || $h["dataset"] != want { next }
    $h["treatment"] == control { next }
    !seen[$h["treatment"]]++ { print $h["treatment"] }' "$samplesheet")
[ -n "$arms" ] || { echo "$samplesheet gives $dataset no arm to contrast with $control" >&2; exit 1; }

for arm in $arms; do
    for estimator in unpaired paired; do
        if [ "$estimator" = paired ]; then paired_flag=(--paired); else paired_flag=(); fi
        echo "=== $arm against $control, $estimator ==="
        riborescue readthrough "$counts" --gtf "$gtf" --samplesheet "$samplesheet" \
            --dataset "$dataset" --treated "$arm" --control "$control" "${paired_flag[@]}" \
            --out "$outdir/${arm}_vs_${control}.${estimator}.tsv"
    done
done
