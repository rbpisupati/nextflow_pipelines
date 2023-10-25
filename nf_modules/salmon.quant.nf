nextflow.enable.dsl=2


params.args_salmon = false
params.args_salmon_prep = false


workflow SALMON {
    take:
    genome
    raw_reads
    outputDir

    main:
    index = SALMON_PREPARE ( genome.fasta, genome.transcript )


        
}



process SALMON_PREPARE {
    tag "$genome_id"
    label "process_medium"

    conda "bioconda::salmon=1.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.9.0--h7e5ed60_1' :
        'quay.io/biocontainers/salmon:1.9.0--h7e5ed60_1' }"

    input:
        tuple val(genome_id), path(fasta)
        path(transcript_fasta)

    output:
        path "${genome_id}_salmon", emit: index

    
    script:
        def get_decoy_ids = "grep '^>' $fasta | cut -d ' ' -f 1 > decoys.txt"
        def gentrome      = "gentrome.fa"
        if (fasta.endsWith('.gz')) {
            get_decoy_ids = "grep '^>' <(gunzip -c $fasta) | cut -d ' ' -f 1 > decoys.txt"
            gentrome      = "gentrome.fa.gz"
        }


        def args = params.args_salmon_prep ?: ''
        """
        $get_decoy_ids
        sed -i.bak -e 's/>//g' decoys.txt
        cat $transcript_fasta $fasta > $gentrome

        salmon index \\
            --threads $task.cpus \\
            -t $gentrome \\
            -d decoys.txt \\
            -i ${genome_id}_salmon
            $args
        """

}


process SALMON_QUANT {
    tag "$name"
    label "process_medium"

    conda "bioconda::salmon=1.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.9.0--h7e5ed60_1' :
        'quay.io/biocontainers/salmon:1.9.0--h7e5ed60_1' }"

    input:
        tuple val(name), path(bam)
        // path  salmon_index
        // path  gtf
        path  transcript_fasta
        val(outputdir)

    output:
        tuple val(name), path("${name}.salmon*") , emit: results


    publishDir "$outputdir", mode: "copy", overwrite: true
    
    script:
    def args = params.args_salmon ?: ''
    args = (args.contains('libType')) ? args : args + ' --libType=IU '
    // -i $salmon_index \\
    """
    salmon quant \\
        $args \\
        --threads $task.cpus \\
        -t $transcript_fasta \\
        -a $bam \\
        -o ${name}.salmon
    """
    

}