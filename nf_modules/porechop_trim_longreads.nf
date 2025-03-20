
// Process to trim adapters from long reads using Porechop
process PORECHOP_PORECHOP {
    tag "$sample_name"
    label 'process_medium'

    conda "bioconda::porechop=0.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39h7cff6ad_2' :
        'biocontainers/porechop:0.2.4--py39h7cff6ad_2' }"

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("*.fastq.gz"), emit: reads
    tuple val(sample_name), path("*.log")     , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_name}"
    """
    porechop \\
        -i $reads \\
        -t $task.cpus \\
        $args \\
        -o ${sample_name}.fastq.gz \\
        > ${sample_name}.log
    
    """

}
