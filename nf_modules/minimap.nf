

process MINIMAP2_ALIGN {
    tag "$sample_name"
    label 'process_high'

    conda "bioconda::minimap2=2.29 bioconda::samtools=1.21"
    // container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/66/66dc96eff11ab80dfd5c044e9b3425f52d818847b9c074794cf0c02bfa781661/data' :
        'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c' }"

    publishDir "${outdir}/$sample_name",
		mode: "copy", overwrite: true,
		saveAs: {filename ->
			if( params.save_intermediate ) filename
			else null
		}

    input:
    tuple val(genome_id), val(genome_name), path(ref_genome)
    tuple val(sample_name), path(reads)
    val outdir
    val minimap2_args

    output:
    tuple val(sample_name), path("*.sam")           , emit: sam
    path("*.log")                                   , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = minimap2_args
    def prefix = task.ext.prefix ?: "${sample_name}"
    """
    minimap2 $args $ref_genome $reads > ${sample_name}.sam 2> ${sample_name}.log
        
    """
}