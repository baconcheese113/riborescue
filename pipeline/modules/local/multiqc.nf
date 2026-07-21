process MULTIQC {
    label 'process_reads'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path reports, stageAs: 'reports/*'

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data'       , emit: data

    script:
    """
    multiqc --no-ansi --force --outdir . reports
    """
}
