nextflow.enable.dsl=2

process SAMTOOLS_IDXSTAT{	
    
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
        'quay.io/biocontainers/samtools:1.14--hb421002_0' }"

	input:
		tuple val(name), path( bam )
		val (outputdir)
		val (verbose)

	output:
		path "*txt",        emit: idxstat

	publishDir "$outputdir",
		mode: "copy", overwrite: true
	
    script:
		if (verbose){
			println ("[MODULE] SAMTOOLS IDXSTAT " )
		}
		"""
		samtools sort $bam -o ${name}_sorted.bam
		samtools idxstat ${name}_sorted.bam > ${name}_idxstat.txt
    	"""
}