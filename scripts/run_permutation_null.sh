#!/usr/bin/env bash
# The corrected shuffle null of ADR-0021: many permutations, synchronised across drugs.
#
# The work is one kinetic refit per drug per round per permutation, and it divides cleanly by
# permutation index. Each shard draws its own disjoint range, so no two shards can draw the same
# mapping and the shards concatenate into one null without deduplication.
#
# One BLAS thread each. The fits are small dense problems; an unbounded BLAS would have one shard
# take every core while the rest queued.

set -euo pipefail

table=${1:-tests/fixtures/kinetics/codon_occupancy.tsv}
permutations=${PERMUTATIONS:-999}
shards=${SHARDS:-14}
outdir=${NULL_OUTDIR:-results/kinetics/permutation_null}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

[ -f "$table" ] || { echo "no codon table at $table" >&2; exit 1; }
mkdir -p "$outdir"

# The last shard absorbs the remainder, so the shards cover exactly `permutations` draws with none
# repeated and none skipped.
per_shard=$((permutations / shards))
remainder=$((permutations % shards))
offset=0
for shard in $(seq 0 $((shards - 1))); do
    count=$per_shard
    [ "$shard" -lt "$remainder" ] && count=$((count + 1))
    [ "$count" -gt 0 ] || continue
    # Only the first shard writes the observed gains; the rest would recompute the same numbers.
    if [ "$shard" = 0 ]; then observed=(--observed); else observed=(--no-observed); fi
    riborescue permutation-null --codon-table "$table" \
        --permutations "$count" --offset "$offset" "${observed[@]}" \
        --out "$outdir/null.$shard.tsv" >"$outdir/null.$shard.log" 2>&1 &
    offset=$((offset + count))
done
echo "launched $shards shards covering $offset permutations of each shuffle family"
wait

head -1 "$outdir/null.0.tsv" >"$outdir/null.tsv"
for shard in $(seq 0 $((shards - 1))); do
    [ -f "$outdir/null.$shard.tsv" ] || continue
    tail -n +2 "$outdir/null.$shard.tsv" >>"$outdir/null.tsv"
done
echo "gathered $(($(wc -l <"$outdir/null.tsv") - 1)) rows into $outdir/null.tsv"
