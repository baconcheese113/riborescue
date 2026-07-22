#!/usr/bin/env bash
# Calibrate one dataset's footprint libraries, then merge that dataset's tables.
#
# One R process per library, sequentially. A transcriptome alignment records a footprint once per
# isoform it fits, so a whole library does not sit in memory beside its neighbours; an exited
# process is the only reliable way to return that memory. Peak resident memory is recorded per
# library rather than assumed from the smallest one.
#
# A dataset owns its own directory and the combiner reads only that directory. Libraries from
# different experiments are never merged by a wildcard, which is how tables from an abandoned run
# reached a combined table once already.

set -euo pipefail

usage() {
    echo "usage: $0 <dataset> <samplesheet> [cell-line] [alignments-dir] [gtf]" >&2
    exit 2
}

dataset=${1:-}
samplesheet=${2:-}
cell_line=${3:-}
alignments=${4:-results/reads/alignments}
gtf=${5:-data/gencode/gencode.v50.primary_assembly.annotation.gtf.gz}
[ -n "$dataset" ] && [ -n "$samplesheet" ] || usage
[ -f "$samplesheet" ] || { echo "no samplesheet at $samplesheet" >&2; exit 1; }

outdir="results/psite/$dataset"
# The annotation and the extension windows depend on the reference alone, so datasets sharing a
# reference share them rather than each spending four minutes rebuilding the same table.
cache="results/psite/_reference/$(basename "$gtf" .gtf.gz)"
mkdir -p "$outdir" "$cache"

if [ ! -f "$cache/extensions.tsv" ]; then
    echo "missing $cache/extensions.tsv — run 'pixi run extensions' first" >&2
    exit 1
fi

# The libraries of this dataset, named by its samplesheet rather than by whatever the alignment
# directory happens to hold. One samplesheet can describe more than one experiment, so the cell line
# narrows it: a shared sheet must not let another experiment's libraries into this dataset.
samples=$(awk -F'\t' -v want="$cell_line" '
    NR==1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    $h["assay"] != "riboseq" { next }
    want != "" && $h["cell_line"] != want { next }
    { print $h["sample"] }' "$samplesheet")
[ -n "$samples" ] || { echo "$samplesheet names no footprint libraries" >&2; exit 1; }

# A previous run of this dataset against a wider samplesheet leaves its per-library tables behind,
# and the combiner would gather those too — the same mixing the per-dataset directory prevents, one
# level further down.
rm -f "$outdir"/*.offsets.tsv "$outdir"/*.frame.tsv "$outdir"/*.region.tsv \
      "$outdir"/*.shift.tsv "$outdir"/*.metaprofile.tsv "$outdir"/*.readthrough.tsv
: >"$outdir/resources.tsv"

# Every library runs the same bytes. Each library is a fresh R process that re-reads the script, so
# an edit part-way through a run is read half-written by whichever process opens it at that moment —
# which happened, and surfaced as a parse error in the middle of a healthy run. The snapshot is
# read-only and removed on exit, and its checksum is recorded beside the results.
snapshot=$(mktemp -d)
trap 'chmod -R u+w "$snapshot"; rm -rf "$snapshot"' EXIT
cp scripts/run_psite.R "$snapshot/run_psite.R"
chmod a-w "$snapshot/run_psite.R"
md5sum "$snapshot/run_psite.R" | awk '{print $1}' >"$outdir/script.md5"

for sample in $samples; do
    bam="$alignments/$sample.Aligned.toTranscriptome.out.bam"
    [ -f "$bam" ] || { echo "no transcriptome alignment for $sample at $bam" >&2; exit 1; }
    /usr/bin/time -f "${sample}\t%M\t%e" -a -o "$outdir/resources.tsv" \
        Rscript "$snapshot/run_psite.R" --bam "$bam" --sample "$sample" --gtf "$gtf" \
            --outdir "$outdir" --cache "$cache"
done

Rscript "$snapshot/run_psite.R" --combine --outdir "$outdir" --evidence "docs/figures/$dataset"

printf 'sample\tpeak_rss_kb\tseconds\n' | cat - "$outdir/resources.tsv" >"$outdir/resources.tmp"
mv "$outdir/resources.tmp" "$outdir/resources.tsv"
