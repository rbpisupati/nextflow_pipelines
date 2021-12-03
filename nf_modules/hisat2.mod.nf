nextflow.enable.dsl=2

process HISAT2 {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	label 'bigMem'
	label 'multiCore'

	conda (params.enable_conda ? "bioconda::hisat2=2.2.0 bioconda::samtools=1.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0' :
        'quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0' }"

    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (hisat2_args)
		val (verbose)

	output:
	    path "*bam",       emit: bam
		path "*stats.txt", emit: stats

	publishDir "$outputdir",
		mode: "copy", overwrite: true

    script:
	
		if (verbose){
			println ("[MODULE] HISAT2 ARGS: " + hisat2_args)
		}
	
		cores = 8
		readString = ""
		hisat_options = hisat2_args

		// Options we add are
		hisat_options = hisat_options + " --no-unal --no-softclip --new-summary"

		if (reads instanceof List) {
			readString = "-1 "+reads[0]+" -2 "+reads[1]
			hisat_options = hisat_options + " --no-mixed --no-discordant"
		}
		else {
			readString = "-U "+reads
		}
		index = params.genome["hisat2"]
		
		// TODO: need to add a check if the splice-site infile is present or not, and leave out this parameter otherwise 
		splices = " --known-splicesite-infile " + params.genome["hisat2_splices"]
		hisat_name = name + "_" + params.genome["name"]

		"""
		hisat2 -p ${cores} ${hisat_options} -x ${index} ${splices} ${readString}  2>${hisat_name}_hisat2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${hisat_name}_hisat2.bam
		"""

}