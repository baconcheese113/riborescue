process TRIM_SUMMARY {
    label 'process_single'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path reports
    path samplesheet

    output:
    path 'trim_summary.tsv', emit: summary

    script:
    """
    riborescue trim-summary ${reports} --samplesheet ${samplesheet} --out trim_summary.tsv
    """
}
