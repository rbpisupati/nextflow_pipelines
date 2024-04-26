nextflow.enable.dsl=2
params.save_intermediate = false

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
		path "*txt",        emit: idxstat

	publishDir "$outputdir",
		mode: "copy", overwrite: true,
		saveAs: {filename ->
			if( params.save_intermediate ) filename
			else null
		}
	
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