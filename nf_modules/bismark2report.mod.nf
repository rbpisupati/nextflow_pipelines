nextflow.enable.dsl=2

process BISMARK2REPORT {

	conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.23.0--0' :
        'quay.io/biocontainers/bismark:0.23.0--0' }"
	
    input:
	    file (file)
		val (outputdir)
		val (bismark2report_args)
		val (verbose)

	output:
	    path "*html",       emit: html
		
	publishDir "$outputdir",
		mode: "copy", overwrite: true

    script:
		if (verbose){
			println ("[MODULE] BISMARK2REPORT ARGS: " + bismark2report_args)
		}

		"""
		bismark2report
		"""

}