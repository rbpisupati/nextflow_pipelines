nextflow.enable.dsl=2
params.samtools_sort_args = [:]

process SAMTOOLS_VIEW {
	tag "$name"
	label 'bigMem'

	// publishDir "$outputdir", mode: "copy", overwrite: true

	conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
        'quay.io/biocontainers/samtools:1.14--hb421002_0' }"


	input:
		tuple val(name), path(sam)
		val(outputdir)
	
	output:
		tuple val(name), path("${name}.bam"), 	emit: bam

	script:
	"""
	samtools view -b -o ${name}.bam -S $sam
	"""
}

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
		tuple val(name), path("*sorted.bam"), path("*sorted.bam.bai"), 	emit: bam

	script:
	def avail_mem = task.memory ? ((task.memory.toGiga() - 6) / task.cpus).trunc() : false
	def sort_mem = avail_mem && avail_mem > 2 ? "-m ${avail_mem}G" : ''
	output_name = bam.baseName.replaceAll(/.deduplicated/,"")
	"""
	samtools sort $params.samtools_sort_args $bam $sort_mem -o ${output_name}.sorted.bam
	samtools index ${output_name}.sorted.bam
	"""

}