#!/usr/bin/env bash
# Calibrate every footprint library, then merge the per-sample tables.
#
# One R process per library, sequentially. A transcriptome alignment records a footprint once per
# isoform it fits, so a whole library does not sit in memory beside its neighbours; an exited
# process is the only reliable way to return that memory. Peak resident memory is recorded per
# library rather than assumed from the smallest one.

set -euo pipefail

gtf=${1:-data/gencode/gencode.v50.primary_assembly.annotation.gtf.gz}
outdir=${2:-results/reads/qc/psite}
mkdir -p "$outdir"

# The extension windows come from `riborescue extensions`, which runs in the default environment.
if [ ! -f "$outdir/extensions.tsv" ]; then
    echo "missing $outdir/extensions.tsv — run 'pixi run extensions' first" >&2
    exit 1
fi
: >"$outdir/resources.tsv"

for bam in results/reads/alignments/*.toTranscriptome.out.bam; do
    sample=$(basename "$bam" .Aligned.toTranscriptome.out.bam)
    /usr/bin/time -f "${sample}\t%M\t%e" -a -o "$outdir/resources.tsv" \
        Rscript scripts/run_psite.R --bam "$bam" --sample "$sample" --gtf "$gtf" \
            --outdir "$outdir"
done

Rscript scripts/run_psite.R --combine --outdir "$outdir"

printf 'sample\tpeak_rss_kb\tseconds\n' | cat - "$outdir/resources.tsv" >"$outdir/resources.tmp"
mv "$outdir/resources.tmp" "$outdir/resources.tsv"
