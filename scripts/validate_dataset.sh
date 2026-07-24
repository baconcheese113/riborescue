#!/usr/bin/env bash
# One dataset, from alignments to both analysis arms, in a single resumable command.
#
# Calibration surveys the whole length range once. The length set is then chosen from that survey by
# ADR-0011, and both arms come out of the same per-transcript counts by summing different lengths —
# so each alignment is read once, not once per arm.
#
# The stages live in two Pixi environments: calibration runs in `psite`, which pins the older R that
# riboWaltz needs, and the analysis runs in the default environment that holds the riborescue CLI.
# Each stage is invoked through `pixi run` so this can be called from a plain shell, and each is
# skippable when its output is already there, so re-running after an interruption continues.

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

export SAMPLESHEET=${SAMPLESHEET:-results/staged_runs.tsv}
outdir=${OUTDIR:-results/psite/$dataset}
manifest="$outdir/calibration.json"
published=${PUBLISHED_LENGTHS:-28 35}

export OUTDIR="$outdir"

echo "== calibrating $dataset over the survey range =="
pixi run -e psite bash scripts/run_psite_all.sh "$dataset" "$SAMPLESHEET" "$alignments"

echo
echo "== choosing the length set =="
pixi run riborescue select-lengths "$dataset" \
    --frames "$outdir/frame_by_length.tsv" \
    --offsets "$outdir/psite_offsets.tsv" \
    --script-md5 "$outdir/script.md5" \
    --out "$manifest"

echo
echo "== primary analysis, on the selected set =="
COUNTS="$outdir/readthrough_counts.tsv" MANIFEST="$manifest" \
    pixi run bash scripts/run_readthrough.sh "$dataset" "$control" "$estimators"

echo
echo "== sensitivity analysis, on the full-length window ($published nt) =="
# Reported whatever the primary showed, and it cannot rescue a primary that failed: the same
# manifest gates it, so a dataset whose libraries did not calibrate has no path through here either.
COUNTS="$outdir/readthrough_counts.tsv" MANIFEST="$manifest" \
    PUBLISHED="$published" OUTDIR_SUFFIX=".published" \
    pixi run bash scripts/run_readthrough.sh "$dataset" "$control" "$estimators"

echo
echo "== $dataset complete =="
echo "  primary   results/readthrough/$dataset/*.unpaired.tsv / *.paired.tsv"
echo "  sensitivity results/readthrough/$dataset/*.published.tsv"
