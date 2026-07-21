process CONTAMINANT_INDEX {
    label 'process_reads'
    storeDir params.contaminant_index_store

    input:
    path transcripts
    path rdna

    output:
    path 'contaminant_index', emit: index

    script:
    // STAR scales its suffix array index to the reference; the default overshoots one this small.
    """
    riborescue contaminants ${transcripts} --include ${rdna} --out contaminants.fa
    mkdir contaminant_index
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir contaminant_index \\
        --genomeFastaFiles contaminants.fa \\
        --genomeSAindexNbases ${params.contaminant_sa_index_nbases}
    rm contaminants.fa
    """
}
