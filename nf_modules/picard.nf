nextflow.enable.dsl=2


process PICARD_DEDUP {
    tag "$name"

    conda (params.enable_conda ? 'bioconda::picard=2.25.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.25.7--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.25.7--hdfd78af_0' }"

    // publishDir "${outdir}", mode: 'copy'
    
    input:
    tuple val(name), path(bam), path(bam_idx)
    // val(outdir)

    output:
    tuple val(name), path("${name}.dedup.bam"), emit: dedup_bam
    tuple val(name), path("${name}.metrics.txt"), emit: metrics

    script:
    """
    picard MarkDuplicates\
        I=$bam O=${name}.dedup.bam METRICS_FILE=${name}.metrics.txt
    """

}

process PICARD_ADDINFO {
    tag "$name"

    conda (params.enable_conda ? 'bioconda::picard=2.25.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.25.7--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.25.7--hdfd78af_0' }"
    
    publishDir "${outdir}", mode: 'copy'
    
    input:
    tuple val(name), path(bam)
    val(outdir)

    output:
    tuple val(name), path("${name}.mkdup.bam"), emit: mkdup_bam

  
    script:
    """
    picard AddOrReplaceReadGroups\
        I=${bam} O=${name}.mkdup.bam\
        ID=$name LB=$name PL=illumina PU=none SM=$name

    """
}
