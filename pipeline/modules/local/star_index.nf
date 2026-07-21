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
    // Sampling every second suffix array entry halves the memory an alignment holds, which buys
    // concurrency across samples at the cost of a slower search within each.
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
