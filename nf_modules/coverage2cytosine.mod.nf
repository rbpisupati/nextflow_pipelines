nextflow.enable.dsl=2

params.singlecell = false
params.rrbs       = false
params.verbose    = false
params.pbat       = false
params.nome       = false

genome = params.genome["bismark"]

process COVERAGE2CYTOSINE {
	tag "$coverage_file" // Adds name to job submission instead of (1), (2) etc.

	conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.23.0--0' :
        'quay.io/biocontainers/bismark:0.23.0--0' }"

	// dynamic directive
	memory { 20.GB * task.attempt }  
	maxRetries 5

    input:
	    path(coverage_file)
		val (outputdir)
		val (coverage2cytosine_args)
		val (verbose)

	output:
	    path "*{report.txt.gz,report.txt}", emit: report
		path "*{.cov.gz,.cov}",             emit: coverage
		path "*cytosine_context_summary.txt", optional: true, emit: summary
	
	publishDir "$outputdir",
		mode: "copy", overwrite: true
    
	script:
		
		// removing the file extension from the input file name 
		// (https://www.nextflow.io/docs/latest/script.html#removing-part-of-a-string)
		outfile_basename = coverage_file.toString()  // Important to convert nextflow.processor.TaskPath object to String first
		outfile_basename = (outfile_basename - ~/.bismark.cov.gz$/)
		outfile_basename = (outfile_basename - ~/.cov.gz$/)
		outfile_basename = (outfile_basename - ~/.cov$/)

		if (verbose){
			println ("[MODULE] BISMARK COVERAGE2CYTOSINE ARGS: " + coverage2cytosine_args)
			println ("Bismark Genome is: " + genome)
		}

		// Options we add are
		cov2cyt_options = coverage2cytosine_args + " --gzip "
		
		if (params.nome){
			if (verbose){
				println ("NOMe-seq outfile basename: $outfile_basename")
			}
			cov2cyt_options += " --nome"
		}

		
		if (verbose){
			println ("Now running command: coverage2cytosine --genome $genome $cov2cyt_options --output ${outfile_basename} $coverage_file ")
		}

		"""
		coverage2cytosine --genome $genome $cov2cyt_options --output ${outfile_basename} $coverage_file
		"""
		
		
}