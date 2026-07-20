process VALIDATE_LABELS {
    tag "${labels.name}"
    label 'process_single'
    publishDir "${params.outdir}/labels", mode: 'copy'

    input:
    path labels

    output:
    path labels, emit: labels
    path 'labels.log', emit: log

    script:
    """
    riborescue validate-labels ${labels} | tee labels.log
    """
}
