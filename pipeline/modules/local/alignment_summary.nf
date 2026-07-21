process ALIGNMENT_SUMMARY {
    label 'process_single'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path logs

    output:
    path 'alignment_summary.tsv', emit: summary

    script:
    """
    riborescue alignment-summary ${logs} --out alignment_summary.tsv
    """
}
