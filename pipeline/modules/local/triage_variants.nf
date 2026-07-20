process TRIAGE_VARIANTS {
    tag "${variants.name}"
    label 'process_single'
    publishDir "${params.outdir}/triage", mode: 'copy'

    input:
    path variants

    output:
    path 'triaged.tsv', emit: triaged

    script:
    """
    riborescue triage-table ${variants} --out triaged.tsv
    """
}
