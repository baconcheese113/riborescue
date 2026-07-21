process TRNA_COVERAGE {
    tag "${contexts.name}"
    label 'process_single'
    publishDir "${params.outdir}/scores", mode: 'copy'

    input:
    path contexts

    output:
    path 'trna_coverage.tsv', emit: coverage

    script:
    """
    riborescue trna-coverage ${contexts} --out trna_coverage.tsv
    """
}
