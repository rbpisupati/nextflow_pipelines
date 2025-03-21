nextflow.enable.dsl=2

// Process to annotate long reads using Prokka
// Mainly used for bacterial genomes and prokaryotic genomes
process PROKKA {
    tag "$assembly_name"
    label 'process_low'

    conda "bioconda::prokka=1.14.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl5321hdfd78af_4' :
        'biocontainers/prokka:1.14.6--pl5321hdfd78af_4' }"

    publishDir "${outdir}/${assembly_name}", mode: 'copy'

    input:
    tuple val(assembly_name), path(fasta)
    val(outdir)
    // path proteins
    // path prodigal_tf

    output:
    tuple val(assembly_name), path("${assembly_name}/*.gff"), emit: gff
    tuple val(assembly_name), path("${assembly_name}/*.gbk"), emit: gbk
    tuple val(assembly_name), path("${assembly_name}/*.fna"), emit: fna
    tuple val(assembly_name), path("${assembly_name}/*.faa"), emit: faa
    tuple val(assembly_name), path("${assembly_name}/*.ffn"), emit: ffn
    tuple val(assembly_name), path("${assembly_name}/*.sqn"), emit: sqn
    tuple val(assembly_name), path("${assembly_name}/*.fsa"), emit: fsa
    tuple val(assembly_name), path("${assembly_name}/*.tbl"), emit: tbl
    tuple val(assembly_name), path("${assembly_name}/*.err"), emit: err
    tuple val(assembly_name), path("${assembly_name}/*.log"), emit: log
    tuple val(assembly_name), path("${assembly_name}/*.txt"), emit: txt
    tuple val(assembly_name), path("${assembly_name}/*.tsv"), emit: tsv
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    // args = proteins ? args + " --proteins ${proteins[0]} "  : args
    // args = prodigal_tf ? args + " --prodigaltf ${prodigal_tf[0]} " : args
    """
    prokka \\
        $args \\
        --cpus $task.cpus \\
        --prefix $assembly_name \\
        $fasta

    """
}
