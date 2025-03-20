


process DRAGONFLYE {
    tag ${sample_name}

    conda "bioconda::dragonflye=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dragonflye:1.1.2--hdfd78af_0' :
        'biocontainers/dragonflye:1.1.2--hdfd78af_0' }"

    publishDir "${params.outdir}/dragonflye_assembly/${sample_name}"

    input:
    tuple val(sample_name), path(shortreads), path(longreads)

    output:
    tuple val(meta), path("*.fa")                                               , emit: contigs
    tuple val(meta), path("dragonflye.log")                                     , emit: log
    tuple val(meta), path("{flye,miniasm,raven}.fasta")                         , emit: raw_contigs
    tuple val(meta), path("{flye,miniasm,raven}-unpolished.gfa"), optional:true , emit: gfa
    tuple val(meta), path("flye-info.txt")                      , optional:true , emit: txt

    script:
    def args    = task.ext.args ?: ''

    // memory
    args = args + " --ram $task.memory.toGiga() "
    //cpus
    args = args + " --cpus $task.cpus "
    // shortreads polishing
    args = shortreads ? args + "--R1 ${shortreads[0]} --R2 ${shortreads[1]}" : args

    """
    dragonflye \\
        $args \\
        --force \\
        --prefix ${sample_name} \\
        --outdir ./ \\
        --reads ${longreads}

    """

}