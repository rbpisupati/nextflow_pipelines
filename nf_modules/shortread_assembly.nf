nextflow.enable.dsl=2

// This module is for short read assembly using SPAdes


process SPADES {
    tag "${sample_name}"
    label 'process_high'

    conda "bioconda::spades=3.11.0"
    container "https://depot.galaxyproject.org/singularity/spades:3.11.0--py27_0"

    publishDir "${outdir}/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(shortreads)
    val outdir
    val args

    output:
    tuple val(sample_name), path("*scaffolds.fasta")                   , emit: assembly
    tuple val(sample_name), path("*contigs.fasta")                     , emit: contigs
    tuple val(sample_name), path("*paths")                             , emit: paths
    tuple val(sample_name), path("*spades.log")                        , emit: log
    tuple val(sample_name), path("params.txt")                        , emit: params

    script:
    // def args    = task.ext.args ?: ''
    // memory
    args = args + " --memory ${task.memory.toGiga()} "
    //cpus
    args = args + " --threads $task.cpus "
    
    """
    spades.py \\
        -o ./ \\
        $args \\
        -1 ${shortreads[0]} \\
        -2 ${shortreads[1]}

    """

}