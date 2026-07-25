process SCORE_VARIANTS {
    tag "${contexts.name}"
    label 'process_annotation'
    publishDir "${params.outdir}/scores", mode: 'copy'

    input:
    path contexts
    path training
    path held_out

    output:
    path 'variant_therapy_scores.tsv', emit: scored

    script:
    // The held-out rounds are staged beside the training tables; the scorer reads each therapy's
    // measured error from its sibling and refuses to score without it.
    """
    riborescue score ${contexts} ${training} --out variant_therapy_scores.tsv
    """
}
