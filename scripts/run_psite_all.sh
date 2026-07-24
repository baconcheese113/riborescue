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
    echo "usage: [LENGTHS=a:b] $0 <dataset> <samplesheet> [alignments-dir] [gtf]" >&2
    exit 2
}

dataset=${1:-}
samplesheet=${2:-}
alignments=${3:-results/reads/alignments}
gtf=${4:-data/gencode/gencode.v50.primary_assembly.annotation.gtf.gz}
# The survey pass keeps everything ADR-0011 looks at; the second keeps what it chose. Passing the
# set explicitly is what makes the two passes distinguishable in the resources table beside them.
lengths=${LENGTHS:-18:40}
[ -n "$dataset" ] && [ -n "$samplesheet" ] || usage
[ -f "$samplesheet" ] || { echo "no samplesheet at $samplesheet" >&2; exit 1; }

outdir=${OUTDIR:-results/psite/$dataset}
# The annotation and the extension windows depend on the reference alone, so datasets sharing a
# reference share them rather than each spending four minutes rebuilding the same table.
cache="results/psite/_reference/$(basename "$gtf" .gtf.gz)"
mkdir -p "$outdir" "$cache"

if [ ! -f "$cache/extensions.tsv" ]; then
    echo "missing $cache/extensions.tsv — run 'pixi run extensions' first" >&2
    exit 1
fi

# The libraries of this dataset, named by its samplesheet rather than by whatever the alignment
# directory happens to hold. The samplesheet describes every experiment, and each library declares
# which one it belongs to — the cell line cannot narrow it, because two studies profile HEK293T.
samples=$(awk -F'\t' -v want="$dataset" '
    NR==1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    $h["assay"] != "riboseq" || $h["dataset"] != want { next }
    { print $h["sample"] }' "$samplesheet")
[ -n "$samples" ] || { echo "$samplesheet names no footprint libraries" >&2; exit 1; }

# A library whose tables are all present was calibrated by an identical script and is not repeated:
# a run stopped part-way through nine libraries continues rather than starting over. Tables belonging
# to no library of this dataset are removed, because the combiner gathers by wildcard and a table
# left by a wider samplesheet would be merged into these results.
keep=$(mktemp)
for sample in $samples; do echo "$sample" >>"$keep"; done
for existing in "$outdir"/*.readthrough.tsv; do
    [ -e "$existing" ] || continue
    name=$(basename "$existing" .readthrough.tsv)
    grep -qxF "$name" "$keep" || rm -f "$outdir/$name".*.tsv "$outdir/$name".*.png
done
[ -f "$outdir/resources.tsv" ] || : >"$outdir/resources.tsv"

# Every library runs the same bytes. Each library is a fresh R process that re-reads the script, so
# an edit part-way through a run is read half-written by whichever process opens it at that moment —
# which happened, and surfaced as a parse error in the middle of a healthy run. The snapshot is
# read-only and removed on exit, and its checksum is recorded beside the results.
snapshot=$(mktemp -d)
trap 'chmod -R u+w "$snapshot"; rm -rf "$snapshot" "$keep"' EXIT
cp scripts/run_psite.R "$snapshot/run_psite.R"
chmod a-w "$snapshot/run_psite.R"
checksum=$(md5sum "$snapshot/run_psite.R" | awk '{print $1}')
if [ -f "$outdir/script.md5" ] && [ "$(cat "$outdir/script.md5")" != "$checksum" ]; then
    echo "the calibration script changed; recalibrating every library of $dataset" >&2
    rm -f "$outdir"/*.tsv "$outdir"/*.done "$outdir"/*.png
    : >"$outdir/resources.tsv"
fi
echo "$checksum" >"$outdir/script.md5"

# What a library's calibration depends on. A marker is honoured only when every one of these is
# unchanged: the script bytes, the footprint lengths, the reference the extension windows came from,
# the annotation, and the alignment itself. A re-aligned BAM or a rebuilt extension table changes the
# fingerprint, so a stale marker cannot let a superseded input be skipped.
ext_md5=$(md5sum "$cache/extensions.tsv" | awk '{print $1}')
gtf_md5=$(md5sum "$gtf" | awk '{print $1}')
fingerprint_of() {
    local bam=$1
    printf '%s|%s|%s|%s|%s\n' \
        "$lengths" "$checksum" "$ext_md5" "$gtf_md5" "$(stat -c %s "$bam")"
}

done_count=0
for sample in $samples; do
    bam="$alignments/$sample.Aligned.toTranscriptome.out.bam"
    [ -f "$bam" ] || { echo "no transcriptome alignment for $sample at $bam" >&2; exit 1; }
    marker="$outdir/$sample.done"
    fingerprint=$(fingerprint_of "$bam")
    if [ -s "$marker" ] && [ "$(cat "$marker")" = "$fingerprint" ]; then
        echo "$sample: already calibrated, inputs unchanged"
        done_count=$((done_count + 1))
        continue
    fi
    rm -f "$outdir/$sample".*.tsv "$marker"
    /usr/bin/time -f "${sample}\t%M\t%e" -a -o "$outdir/resources.tsv" \
        Rscript "$snapshot/run_psite.R" --bam "$bam" --sample "$sample" --gtf "$gtf" \
            --outdir "$outdir" --cache "$cache" --lengths "$lengths"
    # Written by the driver, only after R exits 0 and every table is on disk. It is the authority on
    # whether a library is done, so a run interrupted between two of its six tables is not mistaken
    # for a finished one.
    printf '%s' "$fingerprint" >"$marker"
done
echo "$done_count of $(echo "$samples" | wc -w) libraries were already calibrated"

# The combiner gathers by wildcard, so it must not run on a set that is short a library or holds one
# calibrated against a since-changed input. Every sample of this dataset has to carry a marker whose
# fingerprint matches the run that just finished.
for sample in $samples; do
    bam="$alignments/$sample.Aligned.toTranscriptome.out.bam"
    marker="$outdir/$sample.done"
    if [ ! -s "$marker" ] || [ "$(cat "$marker")" != "$(fingerprint_of "$bam")" ]; then
        echo "$sample has no current calibration; refusing to combine a partial set" >&2
        exit 1
    fi
done

Rscript "$snapshot/run_psite.R" --combine --outdir "$outdir" --evidence "docs/figures/$dataset"

printf 'sample\tpeak_rss_kb\tseconds\n' | cat - "$outdir/resources.tsv" >"$outdir/resources.tmp"
mv "$outdir/resources.tmp" "$outdir/resources.tsv"
