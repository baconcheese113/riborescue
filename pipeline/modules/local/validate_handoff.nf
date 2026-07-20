process VALIDATE_HANDOFF {
    tag "${manifest.name}"
    label 'process_single'
    publishDir "${params.outdir}/handoff", mode: 'copy'

    input:
    path manifest
    path results_root, stageAs: 'upstream'

    output:
    path manifest, emit: manifest
    path 'handoff.log', emit: log

    script:
    """
    riborescue validate-handoff ${manifest} --check-files --results-root upstream | tee handoff.log
    """
}
