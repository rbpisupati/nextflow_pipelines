nextflow.enable.dsl=2
params.nogroup = false

process FASTQC {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'quay.io/biocontainers/fastqc:0.11.9--0' }"

	input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (fastqc_args)
		val (verbose)

	output:
	    tuple val(name), path ("*fastqc*"), emit: all
		path "*.zip",  emit: report
	
	publishDir "$outputdir",
		mode: "copy", overwrite: true

	script:

		if (params.nogroup){
			// println ("ADDING --nogroup: " + fastqc_args)
			fastqc_args += " --nogroup "
		}
		
		if (verbose){
			println ("[MODULE] FASTQC ARGS: "+ fastqc_args)
		}

		"""
		fastqc $fastqc_args -q -t 2 ${reads}
		"""
}