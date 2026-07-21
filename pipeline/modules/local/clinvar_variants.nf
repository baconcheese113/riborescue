process CLINVAR_VARIANTS {
    tag "${vcf.name}"
    label 'process_single'
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    path vcf

    output:
    path 'clinvar_nonsense.tsv', emit: variants

    script:
    """
    riborescue clinvar ${vcf} --out clinvar_nonsense.tsv
    """
}
