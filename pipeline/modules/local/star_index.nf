process STAR_INDEX {
    label 'process_index'
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'
    storeDir params.star_index_store

    input:
    path genome
    path gtf

    output:
    path 'star_index', emit: index

    script:
    // A dense human index needs more memory than this machine has. Sampling every second suffix
    // array entry halves that, at the cost of slower alignment.
    """
    gzip -cd ${genome} > genome.fa
    gzip -cd ${gtf} > annotation.gtf
    mkdir star_index
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles genome.fa \\
        --sjdbGTFfile annotation.gtf \\
        --sjdbOverhang ${params.sjdb_overhang} \\
        --genomeSAsparseD ${params.genome_sa_sparse} \\
        --limitGenomeGenerateRAM ${task.memory.toBytes()}
    rm genome.fa annotation.gtf
    """
}
