nextflow.enable.dsl=2

params.samtools_sort_args = ""
params.samtools_index_options = ""

process SAMTOOLS_SORT {
	tag "$name"
	label 'bigMem'
	publishDir "$outputdir", 
			mode: "copy", overwrite: true

	conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
        'quay.io/biocontainers/samtools:1.14--hb421002_0' }"


	input:
		tuple val(name), path(bam)
		val(outputdir)
	
	output:
		tuple val(name), path("*_sorted.bam"), path("*_sorted.bam.bai"), 	emit: bam

	script:
	"""
	samtools sort $params.samtools_sort_args $bam -o ${name}_sorted.bam
	samtools index $params.samtools_index_options ${name}_sorted.bam
	"""

}