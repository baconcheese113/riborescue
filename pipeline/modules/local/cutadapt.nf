process CUTADAPT {
    tag "${meta.sample}"
    label 'process_reads'
    container 'quay.io/biocontainers/cutadapt:5.2--py312hfabe715_2'
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: '*.fastq.gz'
    publishDir "${params.outdir}/qc/trimming", mode: 'copy', pattern: '*.json'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trimmed.fastq.gz'), emit: reads
    path "${meta.sample}.cutadapt.json"        , emit: report

    script:
    // The second adapter and output only exist for a paired library, so both halves of the
    // invocation are assembled rather than branched on downstream.
    def paired = meta.layout == 'paired'
    def second = paired ? "-A ${meta.adapter_3p_2} -p ${meta.sample}_2.trimmed.fastq.gz" : ''
    def minimum = meta.assay == 'riboseq' ? params.min_length_riboseq : params.min_length_rnaseq
    """
    cutadapt \\
        --cores ${task.cpus} \\
        --cut ${meta.cut_5p} \\
        --adapter ${meta.adapter_3p} \\
        --overlap ${meta.adapter_overlap} \\
        --quality-cutoff ${params.quality_cutoff} \\
        --minimum-length ${minimum} \\
        --json ${meta.sample}.cutadapt.json \\
        ${second} \\
        --output ${meta.sample}_1.trimmed.fastq.gz \\
        ${reads}
    """
}
