#!/usr/bin/env bash
# The contrasts of one dataset, each by the estimators its design admits.
#
# Which estimator applies is fixed per dataset in ADR-0010, before any of its numbers exist, so it
# is named on the command line rather than chosen from the answers. A dataset whose arms are
# separate library preparations gets `unpaired` only: matching replicate numbers across arms are not
# a block, and running `both` there would invent a pairing the design does not have.
#
# The arms are named too, because the control differs by contrast — a suppressor tRNA is tested
# against the construct it was delivered in, not against untreated cells.

set -euo pipefail

usage() {
    echo "usage: $0 <dataset> <control-arm> <unpaired|paired|both> [arm...]" >&2
    exit 2
}

dataset=${1:-}
control=${2:-}
estimators=${3:-}
[ -n "$dataset" ] && [ -n "$control" ] && [ -n "$estimators" ] || usage
shift 3

samplesheet=${SAMPLESHEET:-results/staged_runs.tsv}
counts=${COUNTS:-results/psite/$dataset/readthrough_counts.tsv}
gtf=${GTF:-data/gencode/gencode.v50.primary_assembly.annotation.gtf.gz}
manifest=${MANIFEST:-results/psite/$dataset/calibration.json}
[ -f "$counts" ] || { echo "no counts at $counts — calibrate $dataset first" >&2; exit 1; }
[ -f "$manifest" ] || { echo "no calibration manifest at $manifest" >&2; exit 1; }

case "$estimators" in
    unpaired|paired) run=("$estimators") ;;
    both) run=(unpaired paired) ;;
    *) usage ;;
esac

# Named arms are used as given. Otherwise every other arm of the dataset is contrasted with the
# control, so an arm cannot be quietly left out once its numbers are known.
if [ "$#" -gt 0 ]; then
    arms=("$@")
else
    mapfile -t arms < <(awk -F'\t' -v want="$dataset" -v control="$control" '
        NR==1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
        $h["assay"] != "riboseq" || $h["dataset"] != want || $h["treatment"] == control { next }
        !seen[$h["treatment"]]++ { print $h["treatment"] }' "$samplesheet")
fi
[ "${#arms[@]}" -gt 0 ] || { echo "$dataset has no arm to contrast with $control" >&2; exit 1; }

outdir="results/readthrough/$dataset"
mkdir -p "$outdir"

for arm in "${arms[@]}"; do
    for estimator in "${run[@]}"; do
        if [ "$estimator" = paired ]; then paired_flag=(--paired); else paired_flag=(); fi
        echo "=== $arm against $control, $estimator ==="
        riborescue readthrough "$counts" --gtf "$gtf" --samplesheet "$samplesheet" \
            --dataset "$dataset" --manifest "$manifest" \
            --treated "$arm" --control "$control" "${paired_flag[@]}" \
            --out "$outdir/${arm}_vs_${control}.${estimator}.tsv"
    done
done
