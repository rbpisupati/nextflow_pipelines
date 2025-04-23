nextflow.enable.dsl=2

process MULTIQC {
	
	label 'label_medium'

	conda (params.enable_conda ? 'bioconda::multiqc=1.11' : null)
	// container 'https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0' }"

	// dynamic directive
	memory { 20.GB * task.attempt }  
	maxRetries 3

    input:
	    path (file)
		val ( prefix )
		val (outputdir)
		val (multiqc_args)
		val (verbose)

	output:
	    path "*html",       emit: html
		// path "*stats.txt", emit: stats 

	publishDir "$outputdir",
		mode: "copy", overwrite: true,
		saveAs: {filename -> "${prefix}_${filename}" }

    script:
		
		if (verbose){
			println ("[MODULE] MULTIQC ARGS: " + multiqc_args)
		}

		"""
		multiqc $multiqc_args -x work .
		"""

}