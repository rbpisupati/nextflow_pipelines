nextflow.enable.dsl=2
params.local = ''

process BOWTIE2 {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	conda (params.enable_conda ? 'bioconda::bowtie2=2.4.2 bioconda::samtools=1.11 conda-forge::pigz=2.3.4' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0' :
        'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0' }"

	label 'bigMem'
	label 'multiCore'
		
    input:
		tuple val(genome_id), path(index)
	    tuple val(name), path(reads)
		val (outputdir)
		val (bowtie2_args)
		val (verbose)

	output:
	    tuple val(name), path ("*sam"),        emit: sam
		path "*stats.txt", emit: stats 

	publishDir "$outputdir",
		mode: "copy", overwrite: true,
		saveAs: {filename ->
			if( params.save_intermediate ) filename
			else null
		}

	script:
		if (verbose){
			println ("[MODULE] BOWTIE2 ARGS: " + bowtie2_args)
		}

		cores = 8
		readString = ""

		// Options we add are
		bowtie_options = bowtie2_args
		// bowtie_options +=  " --no-unal " // We don't need unaligned reads in the BAM file
		
		if (params.local == '--local'){
			// println ("Adding option: " + params.local )
			bowtie_options += " ${params.local} " 
		}

		if (reads instanceof List) {
			readString = "-1 " + reads[0] + " -2 " + reads[1]
			// bowtie_options += " --no-discordant --no-mixed " // just output properly paired reads
		}
		else {
			readString = "-U " + reads
		}


		bowtie_name = name + "_" + genome_id

		"""
		bowtie2 -x ${index}/${index} -p ${cores} ${bowtie_options} ${readString} > ${bowtie_name}.sam  2>${bowtie_name}_bowtie2_stats.txt
		"""

}