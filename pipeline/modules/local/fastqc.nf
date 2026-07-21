process FASTQC {
    tag "${meta.sample} ${qc_stage}"
    label 'process_reads'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir path: { "${params.outdir}/qc/${qc_stage}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)
    val qc_stage

    output:
    path '*_fastqc.zip' , emit: zip
    path '*_fastqc.html', emit: html

    script:
    """
    fastqc --threads ${task.cpus} --quiet --outdir . ${reads}
    """
}
