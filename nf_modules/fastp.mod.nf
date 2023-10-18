nextflow.enable.dsl=2
params.nogroup = false
params.single_end = [:]

process FASTP {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	conda (params.enable_conda ? "bioconda::fastp=0.23.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'biocontainers/fastp:0.23.4--h5f740d0_0' }"

	input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (fastp_args)
		val (verbose)

	output:
	    tuple val(name), path ("*fastq.gz"), emit: reads
		path "*.json",  emit: report
		path "*.html",  emit: html
	
	publishDir "$outputdir",
		mode: "copy", overwrite: true

	script:
		
		if (verbose){
			println ("[MODULE] FASTP ARGS: "+ fastp_args)
		}
		def sh_reads = "--in1 ${reads[0]} --out1 ${name}_1.fastp.fastq.gz"
		if (!params.single_end){
			sh_reads = " --in2 ${reads[1]} --out2 ${name}_2.fastp.fastq.gz "
		}

		"""
		fastp \\
            $sh_reads \\
            --json ${name}.fastp.json \\
            --html ${name}.fastp.html \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $args \\
            2> ${name}.fastp.log

		"""
}