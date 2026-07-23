#!/usr/bin/env bash
# One dataset, from alignments to both analysis arms, in a single resumable command.
#
# Calibration surveys the whole length range once. The length set is then chosen from that survey by
# ADR-0011, and both arms come out of the same per-transcript counts by summing different lengths —
# so each alignment is read once, not once per arm.
#
# Every stage is skippable when its output is already there, and the calibration resumes per library,
# so re-running after an interruption continues rather than starting over.

set -euo pipefail

usage() {
    echo "usage: $0 <dataset> <control-arm> <unpaired|paired|both> [alignments-dir]" >&2
    echo "       OUTDIR=... to write elsewhere (used by the smoke run)" >&2
    exit 2
}

dataset=${1:-}
control=${2:-}
estimators=${3:-}
alignments=${4:-results/reads/alignments}
[ -n "$dataset" ] && [ -n "$control" ] && [ -n "$estimators" ] || usage

samplesheet=${SAMPLESHEET:-results/staged_runs.tsv}
outdir=${OUTDIR:-results/psite/$dataset}
manifest="$outdir/calibration.json"
published=${PUBLISHED_LENGTHS:-28 35}

echo "== calibrating $dataset over the survey range =="
OUTDIR="$outdir" bash scripts/run_psite_all.sh "$dataset" "$samplesheet" "$alignments"

echo
echo "== choosing the length set =="
riborescue select-lengths "$dataset" \
    --frames "$outdir/frame_by_length.tsv" \
    --offsets "$outdir/psite_offsets.tsv" \
    --script-md5 "$outdir/script.md5" \
    --out "$manifest"

echo
echo "== primary analysis, on the selected set =="
COUNTS="$outdir/readthrough_counts.tsv" MANIFEST="$manifest" \
    bash scripts/run_readthrough.sh "$dataset" "$control" "$estimators"

echo
echo "== sensitivity analysis, on the published window ($published nt) =="
# Reported whatever the primary showed, and it cannot rescue a primary that failed: the same
# manifest gates it, so a dataset whose libraries did not calibrate has no path through here either.
COUNTS="$outdir/readthrough_counts.tsv" MANIFEST="$manifest" \
    PUBLISHED="$published" OUTDIR_SUFFIX=".published" \
    bash scripts/run_readthrough.sh "$dataset" "$control" "$estimators"

echo
echo "== $dataset complete =="
