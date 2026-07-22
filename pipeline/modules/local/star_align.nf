process STAR_ALIGN {
    tag "${meta.sample}"
    label 'process_align'
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'
    publishDir "${params.outdir}/alignments", mode: 'copy', pattern: '*.bam*'
    publishDir "${params.outdir}/qc/alignment", mode: 'copy', pattern: '*.Log.final.out'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path('*.sorted.bam')            , emit: bam
    tuple val(meta), path('*.toTranscriptome.out.bam'), optional: true, emit: transcriptome
    path '*.Log.final.out'                           , emit: log

    script:
    // Footprints are shorter than the mismatch defaults assume, and a spliced footprint is not a
    // meaningful call at 30 nt, so mismatches are bounded as a fraction of the read. Footprints also
    // get a transcriptome-coordinate alignment, which is what riboWaltz reads to place the P-site
    // without the exon-junction arithmetic a genome coordinate would need.
    def footprint = meta.assay == 'riboseq'
    def tuning = footprint
        ? '--outFilterMismatchNoverLmax 0.07 --alignIntronMax 1000000 --outFilterMultimapNmax 20'
        : '--outFilterMultimapNmax 20'
    def quant = footprint ? 'GeneCounts TranscriptomeSAM' : 'GeneCounts'
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn ${reads} \\
        --readFilesCommand gzip -cd \\
        --outFileNamePrefix ${meta.sample}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS nM \\
        --quantMode ${quant} \\
        ${tuning}
    mv ${meta.sample}.Aligned.sortedByCoord.out.bam ${meta.sample}.sorted.bam
    """
}
