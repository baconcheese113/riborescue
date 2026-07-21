process VARIANT_CONTEXTS {
    tag "${variants.name}"
    label 'process_annotation'
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    path variants
    path annotation
    path transcripts
    path proteins

    output:
    path 'clinvar_contexts.tsv', emit: contexts

    script:
    """
    riborescue contexts ${variants} \\
        --annotation ${annotation} \\
        --transcripts ${transcripts} \\
        --proteins ${proteins} \\
        --out clinvar_contexts.tsv
    """
}
