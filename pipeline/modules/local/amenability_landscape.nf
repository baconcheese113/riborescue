process AMENABILITY_LANDSCAPE {
    tag "${contexts.name}"
    label 'process_annotation'
    publishDir "${params.outdir}/scores", mode: 'copy'

    input:
    path contexts
    path scores

    output:
    path 'amenability_landscape.tsv', emit: landscape
    path 'landscape_summary.tsv', emit: summary

    script:
    """
    riborescue landscape ${contexts} ${scores} \\
        --out amenability_landscape.tsv \\
        --summary landscape_summary.tsv
    """
}
