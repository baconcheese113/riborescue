#!/usr/bin/env bash
# The codon occupancy tables of one dataset's vehicle libraries, per ADR-0020.
#
# Four tables, because two conventions are declared rather than chosen: the A site the features are
# built from and the P site beside it, over the manifest's calibrated lengths and over the published
# 28-35 nt window. All four are reported, and a disagreement between them is a finding about which
# ribosome state carries the signal rather than a reason to prefer one.
#
# Only the named arm contributes. The command checks that against the samplesheet and refuses a
# treated library, so passing the wrong file list is an error rather than a silent circularity.

set -euo pipefail

dataset=${1:-gse144140}
control=${2:-dmso}

samplesheet=${SAMPLESHEET:-results/staged_runs.tsv}
psite=${PSITE_DIR:-results/psite/$dataset}
reference=${REFERENCE_DIR:-results/psite/_reference/gencode.v50.primary_assembly.annotation}
transcripts=${TRANSCRIPTS:-data/gencode/gencode.v50.transcripts.fa.gz}
outdir=${KINETICS_OUTDIR:-results/kinetics/$dataset}

mapfile -t counts < <(awk -F'\t' -v want="$dataset" -v arm="$control" '
    NR==1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    $h["assay"] == "riboseq" && $h["dataset"] == want && $h["treatment"] == arm { print $h["sample"] }
    ' "$samplesheet" | sed "s|^|$psite/|; s|$|.codon.tsv|")
[ "${#counts[@]}" -gt 0 ] || { echo "$samplesheet names no $control libraries in $dataset" >&2; exit 1; }
for file in "${counts[@]}"; do
    [ -f "$file" ] || { echo "no codon counts at $file — calibrate $dataset first" >&2; exit 1; }
done

mkdir -p "$outdir"

# build <site> <window-label> [extra riborescue arguments...]
build() {
    local site=$1 label=$2
    shift 2
    echo "=== $site site, $label window ==="
    riborescue codon-occupancy "${counts[@]}" \
        --annotation "$reference/annotation.tsv" --transcripts "$transcripts" \
        --dataset "$dataset" --manifest "$psite/calibration.json" \
        --samplesheet "$samplesheet" --control-arm "$control" --site "$site" "$@" \
        --out "$outdir/codon_occupancy.$site.$label.tsv"
}

build a selected
echo
build p selected
echo
build a published --published-lengths 28 35
echo
build p published --published-lengths 28 35
