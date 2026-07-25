#!/usr/bin/env bash
# The kinetics head-to-head of ADR-0020, over every split in its evaluation grid.
#
# Three splits, each fitting the four models plain and the three kinetic ones under each shuffle.
# All of them run: which one carries the claim is fixed in the ADR — grouped by gene — and running
# only that one would leave no way to see whether a gain is capacity rather than signal.
#
# Drugs run in parallel because the fits are independent and the saturated comparator is slow. The
# width is bounded rather than unlimited: each saturated fit holds a design matrix of a few hundred
# megabytes, and six at once is more memory than the gain in wall clock is worth.

set -euo pipefail

table=${1:-results/kinetics/gse144140/codon_occupancy.a.published.tsv}
site=${2:-a}

oracle=${ORACLE:-tests/fixtures/oracle}
outdir=${HEAD_TO_HEAD_OUTDIR:-results/kinetics/head_to_head}
width=${WIDTH:-3}

# Each fit is already a dense linear algebra problem, so an unbounded BLAS would run one drug across
# every core and leave the others waiting. One thread each keeps the parallelism where it is useful.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}

[ -f "$table" ] || { echo "no codon table at $table — run 'pixi run codon-occupancy' first" >&2; exit 1; }
mkdir -p "$outdir"

mapfile -t drugs < <(
    find "$oracle" -name 'features_*.tsv.gz' -printf '%f\n' | sed 's/^features_//; s/\.tsv\.gz$//' | sort
)
[ "${#drugs[@]}" -gt 0 ] || { echo "$oracle holds no per-drug features" >&2; exit 1; }

for config in published_random_cv grouped_by_gene grouped_by_sequence_cluster; do
    echo "=== $config, $site site ==="
    running=0
    for drug in "${drugs[@]}"; do
        riborescue head-to-head --codon-table "$table" --oracle "$oracle" --site "$site" \
            --config "$config" --drug "$drug" \
            --out "$outdir/$config.$site.$drug.tsv" >"$outdir/$config.$site.$drug.log" 2>&1 &
        running=$((running + 1))
        if [ "$running" -ge "$width" ]; then wait -n; running=$((running - 1)); fi
    done
    wait
    for drug in "${drugs[@]}"; do
        tail -n +2 "$outdir/$config.$site.$drug.log" | grep -E '^ ' || true
    done
done

# One table per artefact, gathered by header-aware concatenation rather than by a wildcard read: a
# file left by an earlier site or an abandoned config would otherwise be merged into these.
gather() {
    local suffix=$1 destination=$2
    local first=1
    : >"$destination"
    for config in published_random_cv grouped_by_gene grouped_by_sequence_cluster; do
        for drug in "${drugs[@]}"; do
            local part="$outdir/$config.$site.$drug$suffix.tsv"
            [ -f "$part" ] || { echo "missing $part; refusing to gather a partial grid" >&2; exit 1; }
            if [ "$first" = 1 ]; then cat "$part" >>"$destination"; first=0
            else tail -n +2 "$part" >>"$destination"; fi
        done
    done
}

echo
gather "" "$outdir/rounds.$site.tsv"
gather "_intervals" "$outdir/intervals.$site.tsv"
echo "gathered into $outdir/rounds.$site.tsv and $outdir/intervals.$site.tsv"
