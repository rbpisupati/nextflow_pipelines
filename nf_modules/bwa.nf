nextflow.enable.dsl=2
params.local = ''

process BWA_MEM {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	conda (params.enable_conda ? 'bioconda::bwa=0.7.18' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa%3A0.7.8--hed695b0_5' :
        'quay.io/biocontainers/bwa:0.7.17--he4a0461_11' }"

	label 'bigMem'
	label 'quadCore'
		
    input:
		tuple val(genome_id), path(index)
	    tuple val(name), path(reads)
		val (outputdir)
		val (bwa_args)
		val (verbose)

	output:
	    tuple val(name), path ("*sam"),        emit: sam 

    // we do not need to save sam files
	//  publishDir "$outputdir",
	// 	mode: "copy", overwrite: true
	publishDir "$outputdir",
		mode: "copy", overwrite: true,
		saveAs: {filename ->
			if( params.save_intermediate ) filename
			else null
		}

	script:
		if (verbose){
			println ("[MODULE] BWA ARGS: " + bwa_args)
		}

		readString = ""

		// Options we add are
		bwa_options = bwa_args

		// if (reads instanceof List) {
		// 	readString = "-1 " + reads[0] + " -2 " + reads[1]
		// 	bowtie_options += " --no-discordant --no-mixed " // just output properly paired reads
		// }
		// else {
		// 	readString = "-U " + reads
		// }

		"""
        bwa mem $bwa_options -t ${task.cpus} $genome_id/genome.fa $reads > ${name}.sam
        """
        // bowtie2 -x ${index}/${index} -p ${cores} ${bowtie_options} ${readString}  2>${bowtie_name}_bowtie2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${bowtie_name}_bowtie2.bam
		// 

}