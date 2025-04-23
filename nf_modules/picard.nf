nextflow.enable.dsl=2

process PICARD_COLLECTMULTIPLEMETRICS {
    tag "$name"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::picard=2.25.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.25.7--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.25.7--hdfd78af_0' }"

    publishDir "${outputdir}/${name}", mode: 'copy'

    input:
    tuple val(name), path(bam), path(bam_idx)
    tuple val(id), val(genome_id), path(ref_fasta)
    val(outputdir)

    output:
    tuple val(name), path("*_metrics"),     emit: metrics
    tuple val(name), path("*pdf"),          emit: pdf, optional: true
    path("*_metrics"),                      emit: report

    script:
    avail_mem = (task.memory.mega*0.8).intValue()
    """
    picard \\
        -Xmx${avail_mem}M \\
        CollectMultipleMetrics \\
        --REFERENCE_SEQUENCE ${ref_fasta} \\
        --INPUT $bam \\
        --OUTPUT ${name}.CollectMultipleMetrics \\

    """
}

process PICARD_DEDUP {
    tag "$name"

    conda (params.enable_conda ? 'bioconda::picard=2.25.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.25.7--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.25.7--hdfd78af_0' }"

    publishDir "$outputdir",
		mode: "copy", overwrite: true,
		saveAs: {filename ->
			if( params.save_intermediate ) filename
			else null
		}
    
    input:
    tuple val(name), path(bam), path(bam_idx)
    val(outputdir)

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
    
    publishDir "$outputdir",
		mode: "copy", overwrite: true,
		saveAs: {filename ->
			if( params.save_intermediate ) filename
			else null
		}
    
    input:
    tuple val(name), path(bam)
    val(outputdir)

    output:
    tuple val(name), path("${name}.mkdup.bam"), emit: mkdup_bam

  
    script:
    """
    picard AddOrReplaceReadGroups\
        I=${bam} O=${name}.mkdup.bam\
        ID=$name LB=$name PL=illumina PU=none SM=$name

    """
}
