process DEPLETE_CONTAMINANTS {
    tag "${meta.sample}"
    label 'process_reads'
    publishDir "${params.outdir}/qc/depletion", mode: 'copy', pattern: '*.Log.final.out'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path('*.depleted*.fastq.gz'), emit: reads
    path '*.Log.final.out'                       , emit: log

    script:
    // What fails to match the structural RNA is what carries the coding signal, so the reads kept
    // here are the ones STAR could not align.
    def rename = meta.layout == 'paired'
        ? "mv ${meta.sample}.Unmapped.out.mate1 ${meta.sample}.depleted_1.fastq \
           && mv ${meta.sample}.Unmapped.out.mate2 ${meta.sample}.depleted_2.fastq"
        : "mv ${meta.sample}.Unmapped.out.mate1 ${meta.sample}.depleted_1.fastq"
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn ${reads} \\
        --readFilesCommand gzip -cd \\
        --outFileNamePrefix ${meta.sample}. \\
        --outSAMtype None \\
        --outReadsUnmapped Fastx
    ${rename}
    gzip ${meta.sample}.depleted_*.fastq
    """
}
