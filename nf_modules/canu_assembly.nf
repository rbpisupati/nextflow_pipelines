// Multiple assemblies are already present
// 1. Canu: https://github.com/marbl/canu
// single molecule sequence assembly for large and small genomes

// 2. Flye
// 6. Miniasm
// 7. Miniasm2
// 8. Wtdbg2: Redbean: A fuzzy Bruijn graph approach to long noisy reads assembly

// Falcon : https://github.com/PacificBiosciences

// https://nf-co.re/bacass/2.4.0/parameters/



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