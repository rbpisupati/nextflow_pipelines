
// SNIFFLES to call structural variants from long-read data
process SNIFFLES {
    tag "$sample_name"
    label 'process_medium'

    conda "bioconda::sniffles=2.4"
    // container 'https://depot.galaxyproject.org/singularity/sniffles:2.4--pyhdfd78af_0'
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sniffles:2.4--pyhdfd78af_0' :
        'biocontainers/sniffles:2.4--pyhdfd78af_0' }"

    input:
    tuple val(sample_name), path(input_bam), path(input_bam_index)
    tuple val(genome_name), val(genome_id), path(reference_fasta)
    val output_dir
    val sniffles_args
    // val(vcf_output)
    // val(snf_output)

    publishDir "${output_dir}", mode: 'copy', overwrite: true

    output:
    tuple val(sample_name), path("${sample_name}.vcf.gz"), path("${sample_name}.vcf.gz.tbi"),       emit: vcf, optional: true
    tuple val(sample_name), path("*.snf"),                                                          emit: snf, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    def args = sniffles_args ? "${sniffles_args}" : ''
    
    def prefix = task.ext.prefix ?: "${sample_name}"
    // def tandem_repeats = tandem_file ? "--tandem-repeats ${tandem_file}" : ''

    """
    sniffles \\
        --input $input_bam \\
        --reference ${reference_fasta} \\
        -t $task.cpus \\
        --vcf ${sample_name}.vcf.gz \\
        --snf ${sample_name}.snf

    """
    
}