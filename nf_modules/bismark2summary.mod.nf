nextflow.enable.dsl=2

process BISMARK2SUMMARY {

	conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.23.0--0' :
        'quay.io/biocontainers/bismark:0.23.0--0' }"
    
	input:
	    file (file)
		val (outputdir)
		val (bismark2summary_args)
		val (verbose)

	output:
	    path "*html",       emit: html
		path "*txt",        emit: report 

	publishDir "$outputdir",
		mode: "copy", overwrite: true

    script:
		// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
		// This is only a temporary workaround until Paolo has fixed the Nextflow bug.
		// https://github.com/nextflow-io/nextflow/issues/1519
		bismark2summary_args = bismark2summary_args.replaceAll(/'/,"")
		if (verbose){
			println ("[MODULE] BISMARK2SUMMARY ARGS: " + bismark2summary_args)
		}

		"""
		bismark2summary
		"""

}