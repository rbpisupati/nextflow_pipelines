nextflow.enable.dsl=2

params.singlecell = false
params.rrbs       = false
params.verbose    = false
params.pbat       = false
params.nonCG      = true
params.save_intermediate = false
params.bismark_methylation_extractor_args = " "

process BISMARK_METHYLATION_EXTRACTOR {
	label 'bigMem'          // 20G
	label 'quadCore'        // 4 cores
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.23.0--0' :
        'quay.io/biocontainers/bismark:0.23.0--0' }"

    input:
	    tuple val (name), path(bam), val(genome_name), path (genome)
		val (outputdir)
		val (verbose)

	output:
	    tuple val (name), path ("CpG*"),        		emit: context_files_CG
		tuple val(name), path("*CX_report.txt.gz"),		emit: cx_report
		path "CH*",                             		emit: context_files_nonCG
		path "*report.txt",                     		emit: report
		path "*M-bias.txt",                     		emit: mbias
		path "*cov.gz",                         		emit: coverage
	
	publishDir "$outputdir", mode: "copy", overwrite: true,
		saveAs: {filename ->
			if( filename.indexOf("CX_report.txt.gz") > 0 ) "cx_report/$filename"
			else if( filename.indexOf("report.txt") > 0 ) "report/$filename"
			else if( filename.indexOf("M-bias.txt") > 0 && params.save_intermediate ) "bias/$filename"
			else if( filename.indexOf("cov.gz") > 0 && params.save_intermediate ) "coverage/$filename"
			// else if( filename.indexOf("fq.gz" ) > 0) "unmapped/$filename"
			// else if( params.save_intermediate ) filename
			else null
		}

	script:
		
		if (verbose){
			println ("[MODULE] BISMARK METHYLATION EXTRACTOR ARGS: " + params.bismark_methylation_extractor_args)
		}


		// Options we add are
		methXtract_options = params.bismark_methylation_extractor_args + " --gzip --genome_folder $genome --cytosine_report --bedGraph  "
		
		if (params.singlecell){
			// println ("FLAG SINGLE CELL SPECIFIED: PROCESSING ACCORDINGLY")
		}

		if (params.nonCG){
			if (verbose){
				println ("FLAG nonCG specified: adding flag --CX ")
			}
			methXtract_options +=  " --CX "
		}

		isPE = isPairedEnd(bam)
		if (isPE){
			// not perform any ignoring behaviour for RRBS or single-cell libraries
			if (!params.rrbs && !params.singlecell && !params.pbat){
				// default ignore parameters for paired-end libraries
				methXtract_options +=  " --ignore_r2 2 "
			}
		}
		else{
			// println("File seems to be single-end")
		}

		// println ("Now running command: bismark_methylation_extractor -parallel ${cores} ${methXtract_options} ${bam}")
		"""
		bismark_methylation_extractor --buffer 10G -parallel $task.cpus ${methXtract_options} ${bam}
		"""

}


def isPairedEnd(bamfile) {

	// need to transform the nextflow.processor.TaskPath object to String
	bamfile = bamfile.toString()
	if (params.verbose){
		println ("Processing file: " + bamfile)
	}
	
	if (bamfile =~ /_pe/){
		if (params.verbose){
			println ("File is paired-end!")
		}
		return true
	}
	else{
	 	if (params.verbose){
			 println ("File is single-end")
		 }
		return false
	}
}