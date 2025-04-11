nextflow.enable.dsl=2

// Multiple assemblies are already present
// 1. Canu: https://github.com/marbl/canu
// single molecule sequence assembly for large and small genomes
process CANU {
    tag "${sample_name}"
    label 'process_high'

    conda "bioconda::canu=2.2"
    container "https://depot.galaxyproject.org/singularity/canu:2.2--ha47f30e_0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/canu:2.2--ha47f30e_0':
    //     'biocontainers/canu:2.2--ha47f30e_0' }"


    publishDir "${outdir}", mode: 'copy'
    
    input:
    tuple val(sample_name), path(reads)
    val outdir
    val args

    output:
    tuple val(sample_name), path("$sample_name/*.report")                   , emit: report
    tuple val(sample_name), path("$sample_name/*.contigs.fasta.gz")         , emit: assembly                , optional: true
    tuple val(sample_name), path("$sample_name/*.unassembled.fasta.gz")     , emit: contigs
    tuple val(sample_name), path("$sample_name/*.correctedReads.fasta.gz")	, emit: corrected_reads         , optional: true
    tuple val(sample_name), path("$sample_name/*.trimmedReads.fasta.gz")	, emit: corrected_trimmed_reads , optional: true
    tuple val(sample_name), path("$sample_name/*.contigs.layout")           , emit: metadata                , optional: true
    tuple val(sample_name), path("$sample_name/*.contigs.layout.readToTig") , emit: contig_position         , optional: true
    tuple val(sample_name), path("$sample_name/*.contigs.layout.tigInfo")   , emit: contig_info             , optional: true

    script:
    // def args = task.ext.args ?: ''
    // def valid_mode = ["-pacbio", "-nanopore", "-pacbio-hifi"]
    // if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Canu. Options: ${valid_mode.join(', ')}" }
    // genomeSize=${genome_size} \\

    """
    canu \\
        -p ${sample_name} \\
        $args \\
        -d ${sample_name} \\
        $reads

    cd ${sample_name} && gzip *.fasta
    """
}


// 2. Flye
// Dragonflye is a hybrid assembler that uses both short and long reads to assemble genomes.
process DRAGONFLYE {
    tag "${sample_name}"
    label 'process_high'

    conda "bioconda::dragonflye=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dragonflye:1.1.2--hdfd78af_0' :
        'biocontainers/dragonflye:1.1.2--hdfd78af_0' }"

    publishDir "${params.outdir}/assembly_dragonflye/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(shortreads), path(longreads)

    output:
    tuple val(sample_name), path("*.fa")                                               , emit: contigs
    tuple val(sample_name), path("dragonflye.log")                                     , emit: log
    tuple val(sample_name), path("{flye,miniasm,raven}.fasta")                         , emit: raw_contigs
    tuple val(sample_name), path("{flye,miniasm,raven}-unpolished.gfa"), optional:true , emit: gfa
    tuple val(sample_name), path("flye-info.txt")                      , optional:true , emit: txt

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


// 3. Circlator
// Tool for assembling circular genomes using long reads
// Sanger Labs
// Run Circlator
process CIRCLATOR {
    tag "${sample_name}"
    label 'process_high'

    conda "bioconda::circlator=1.5.2--py35_1"
    container "quay.io/biocontainers/circlator:1.5.2--py35_1"

    publishDir "${params.outdir}/assembly_circlator/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(shortreads), path(longreads)

    output:
    tuple val(sample_name), path("*.fa")                                               , emit: contigs

    script:
    def args    = task.ext.args ?: ''

    // memory
    args = args + " --ram $task.memory.toGiga() "
    //cpus
    args = args + " --threads $task.cpus "
    // --merge_min_id ${params.merge_min_id} \
    //     --merge_breaklen ${params.merge_breaklen} \
    //     --b2r_length_cutoff ${params.b2r_length_cutoff} \
    //     --merge_min_length_merge ${params.merge_min_length_merge} \
    //     --merge_reassemble_end ${params.merge_reassemble_end} \
    //     --merge_ref_end ${params.merge_ref_end} \
    """
    circlator all \
        $args \
        ${fasta} \
        $longreads \
        ${sample_name}

    """

}



// 6. Miniasm
// 7. Miniasm2
// 8. Wtdbg2: Redbean: A fuzzy Bruijn graph approach to long noisy reads assembly

// Falcon : https://github.com/PacificBiosciences

// https://nf-co.re/bacass/2.4.0/parameters/



