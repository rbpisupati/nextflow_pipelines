nextflow.enable.dsl=2
params.save_intermediate = false
params.samtools_sort_args = ''


process SAMTOOLS_INDEX {
	tag "$name"
	label 'process_medium'

	publishDir "$outputdir", mode: "copy", overwrite: true

	conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
        'quay.io/biocontainers/samtools:1.14--hb421002_0' }"


	input:
		tuple val(name), path(bam)
		val(outputdir)
	
	output:
		tuple val(name), path("*bam"), path("*bam.bai"), 	emit: bam

	script:
	"""
	samtools index $bam
	"""
}

process SAMTOOLS_VIEW {
	tag "$name"
	label 'bigMem'

	// publishDir "$outputdir", mode: "copy", overwrite: true
	publishDir "$outputdir",
		mode: "copy", overwrite: true,
		saveAs: {filename ->
			if( params.save_intermediate ) filename
			else null
		}

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
	
	publishDir "$outputdir/$name", mode: "copy", overwrite: true
	// publishDir "$outputdir",
	// 	mode: "copy", overwrite: true,
	// 	saveAs: {filename ->
	// 		if( params.save_intermediate ) filename
	// 		else null
	// 	}

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

	sort_args = params.samtools_sort_args
    sort_args = (sort_args.contains('-m')) ? sort_args : sort_args + sort_mem

	sort_args = sort_args + " -o ${name}.sorted.bam "
	
	// output_name = bam.baseName.replaceAll(/.deduplicated/,"")
	"""
	samtools sort $sort_args $bam
	samtools index ${name}.sorted.bam
	"""

}



process SAMTOOLS_STATS{	
    
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
        'quay.io/biocontainers/samtools:1.14--hb421002_0' }"

	input:
		tuple val(name), path( bam ), path(bam_idx)
		val (outputdir)
		val (verbose)

	output:
		path "*txt",        emit: stats

	publishDir "${outputdir}/${name}", mode: "copy"
	
    script:
		if (verbose){
			println ("[MODULE] SAMTOOLS IDXSTAT " )
		}
		"""
		samtools idxstat $bam > ${name}_idxstat.txt
		samtools flagstat ${bam} > ${name}_flagstat_report.txt
		samtools stats ${bam} > ${name}_stats_report.txt
    	"""
}