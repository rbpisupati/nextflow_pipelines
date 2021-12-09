nextflow.enable.dsl=2

process BISMARK_DEDUPLICATION {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.23.0--0' :
        'quay.io/biocontainers/bismark:0.23.0--0' }"

	// dynamic directive to increase memory as required
	cpus = 1
	memory { 20.GB * task.attempt }  
  	maxRetries 5
  	
    input:
	    tuple val(name), path(bam)
		val (outputdir)
		val (deduplicate_bismark_args)
		val (verbose)

	output:
		path "*report.txt", emit: report
		tuple val(name), path ("*bam"),        emit: bam

	publishDir "$outputdir",
		mode: "copy", overwrite: true

    script:
		if (verbose){
			println ("[MODULE] BISMARK DEDUPLICATION ARGS: " + deduplicate_bismark_args)
		}

		// Options we add are
		deduplication_options = deduplicate_bismark_args
		deduplication_options += " --bam "

		"""
		deduplicate_bismark ${deduplication_options} ${bam}
		"""

}