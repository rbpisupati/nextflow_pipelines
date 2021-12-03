nextflow.enable.dsl=2

process MULTIQC {
	
	label 'quadCore'
	conda (params.enable_conda ? 'bioconda::multiqc=1.11' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0' }"

	// dynamic directive
	memory { 20.GB * task.attempt }  
	maxRetries 3

    input:
	    path (file)
		val (outputdir)
		val (multiqc_args)
		val (verbose)

	output:
	    path "*html",       emit: html
		// path "*stats.txt", emit: stats 

	publishDir "$outputdir",
		mode: "copy", overwrite: true

    script:
		
		if (verbose){
			println ("[MODULE] MULTIQC ARGS: " + multiqc_args)
		}

		"""
		multiqc $multiqc_args -x work .
		"""

}