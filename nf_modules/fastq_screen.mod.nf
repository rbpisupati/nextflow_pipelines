nextflow.enable.dsl=2
params.bisulfite = ''
params.single_end = false
params.fastq_screen_subset = 500000

process FASTQ_SCREEN {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	// label 'hugeMem'
	
	label 'multiCore'

	conda ( "/users/rahul.pisupati/.conda/envs/fastq_screen/" )
	// conda (params.enable_conda ? "bioconda::fastq-screen=0.13.0" : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
	// 	'https://depot.galaxyproject.org/singularity/fastq-screen%3A0.13.0--pl526_1' :
    //     'quay.io/biocontainers/fastq-screen:0.13.0' }"
	
	memory { 30.GB * task.attempt }  
  	maxRetries 3
  	
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (fastq_screen_args)
		val (verbose)

	output:
	    path "*png",  emit: png
	    path "*html", emit: html
		path "*txt",  emit: report

	publishDir "$outputdir",
		mode: "copy", overwrite: true

    script:
		fastq_screen_args += " --threads $task.cpus "
		fastq_screen_args += " --subset ${params.fastq_screen_subset} "
		if (verbose){
			println ("[MODULE] FASTQ SCREEN ARGS: "+ fastq_screen_args)
		}

		if (params.single_end){
			// TODO: Add single-end parameter
		}
		else{
			// for paired-end files we only use Read 1 (as Read 2 tends to show the exact same thing)
			if (reads instanceof List) {
				reads = reads[0]
			}
		}
		if (params.bisulfite){
			// println("Setting --bisulfite")
			fastq_screen_args += " --bisulfite "
			// println (fastq_screen_args)
		}	

	"""
	fastq_screen $fastq_screen_args $reads
	"""

}