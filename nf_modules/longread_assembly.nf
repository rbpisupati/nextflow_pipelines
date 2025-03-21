nextflow.enable.dsl=2

// Multiple assemblies are already present
// 1. Canu: https://github.com/marbl/canu
// single molecule sequence assembly for large and small genomes



process CANU {
    tag ${sample_name}

    conda "bioconda::canu=2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/canu:2.2--ha47f30e_0':
        'biocontainers/canu:2.2--ha47f30e_0' }"


    publishDir "${params.outdir}/canu_assembly/${sample_name}"
    
    input:
    tuple val(sample_name), path(reads)
    val mode
    val genome_size

    output:
    tuple val(sample_name), path("*.report")                   , emit: report
    tuple val(sample_name), path("*.contigs.fasta.gz")         , emit: assembly                , optional: true
    tuple val(sample_name), path("*.unassembled.fasta.gz")     , emit: contigs
    tuple val(sample_name), path("*.correctedReads.fasta.gz")	, emit: corrected_reads         , optional: true
    tuple val(sample_name), path("*.trimmedReads.fasta.gz")	, emit: corrected_trimmed_reads , optional: true
    tuple val(sample_name), path("*.contigs.layout")           , emit: metadata                , optional: true
    tuple val(sample_name), path("*.contigs.layout.readToTig") , emit: contig_position         , optional: true
    tuple val(sample_name), path("*.contigs.layout.tigInfo")   , emit: contig_info             , optional: true

    script:
    def args = task.ext.args ?: ''
    def valid_mode = ["-pacbio", "-nanopore", "-pacbio-hifi"]
    if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Canu. Options: ${valid_mode.join(', ')}" }

    """
    canu \\
        -p ${sample_name} \\
        $mode \\
        -d ${params.outdir}/canu_assembly/${sample_name} \\
        genomeSize=${genome_size} \\
        $reads

    gzip *.fasta
    """
}


// 2. Flye
// Dragonflye is a hybrid assembler that uses both short and long reads to assemble genomes.
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


// 6. Miniasm
// 7. Miniasm2
// 8. Wtdbg2: Redbean: A fuzzy Bruijn graph approach to long noisy reads assembly

// Falcon : https://github.com/PacificBiosciences

// https://nf-co.re/bacass/2.4.0/parameters/



