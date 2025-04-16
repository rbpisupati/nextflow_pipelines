nextflow.enable.dsl=2

// Process to remove duplicate IDs in the fastq files
process SEQKIT_RMDUP {
    tag "$sample_name"
    label 'process_medium'

    conda "bioconda::seqkit=0.7.1"
    container 'https://depot.galaxyproject.org/singularity/seqkit:0.7.1--0'

    publishDir "${outdir}/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)
    val outdir

    output:
    tuple val(sample_name), path("*.fastq.gz"), emit: reads
    tuple val(sample_name), path("*.log")     , emit: log

    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_name}"
    """
    zcat $reads | seqkit rmdup  -i -o ${sample_name}.rmdup.fastq.gz -d ${sample_name}.duplicated.fastq.gz -D duplicated.detail.txt
    """

}

// Process to trim adapters from long reads using Porechop
process PORECHOP {
    tag "$sample_name"
    label 'process_medium'

    conda "bioconda::porechop=0.2.4"
    container 'https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39h7cff6ad_2'
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39h7cff6ad_2' :
    //     'biocontainers/porechop:0.2.4--py39h7cff6ad_2' }"

    publishDir "${outdir}/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)
    val outdir

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


// PycoQC report for long reads
process PYCOQC {
    tag "$sample_name"
    label 'process_medium'

    conda "bioconda::pycoqc=2.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pycoqc:2.5.2--py_0' :
        'biocontainers/pycoqc:2.5.2--py_0' }"

    input:
    tuple val(sample_name), path(fast5)

    output:
    tuple val(sample_name), path("*.html"), emit: html
    tuple val(sample_name), path("*.json"), emit: json
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def run_summary = file("${fast5}/sequencing_summary.txt").exists() ? "cp ${fast5}/sequencing_summary.txt ./sequencing_summary.txt" : "Fast5_to_seq_summary -f $fast5 -t ${task.cpus} -s './sequencing_summary.txt' --verbose_level 2"
    def barcode_me  = file("${fast5}/barcoding_sequencing.txt").exists() ? "-b ${fast5}/barcoding_sequencing.txt" : ''

    """
    $run_summary

    pycoQC \\
        $args \\
        -f "sequencing_summary.txt" \\
        $barcode_me \\
        -o ${sample_name}.html \\
        -j ${sample_name}.json

    """
}
