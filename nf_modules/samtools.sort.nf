nextflow.enable.dsl=2

params.samtools_index_options = ""
params.samtools_sort_args = ""

process SAMTOOLS_SORT {
	tag "$name"
	label "bigMem"

	conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
        'quay.io/biocontainers/samtools:1.14--hb421002_0' }"

	input:
		tuple val(name), path(bam)
		val(outputdir)
	
	output:
		tuple val(name), path "*_sorted.bam", path "*_sorted.bam.bai", 	emit: bam

	publishDir "$outputdir", mode: "copy"

	script:
	"""
	samtools sort $samtools_sort_args $bam -o ${name}_sorted.bam
	samtools index $samtools_index_options ${name}_sorted.bam
	"""
}

// process SAMTOOLS_SORT {

// 	tag "$name" // Adds name to job submission instead of (1), (2) etc.
// 	label 'bigMem' // 20GB

// 	conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
//         'quay.io/biocontainers/samtools:1.14--hb421002_0' }"

// 	input:
// 		tuple val(name), path( bam )
// 		val (outputdir)
// 		val (verbose)

// 	output:
// 		tuple val(name), path "${name}_sorted.bam", path "${name}_sorted.bam.bai",        emit: bam

// 	publishDir "$outputdir",
// 		mode: "copy", overwrite: true

	
//     script:		
// 		if (verbose){
// 			println ("[MODULE] SAMTOOLS SORT ARGS: " + samtools_sort_args)
// 		}
// 		"""
// 		samtools sort $samtools_sort_args $bam -o ${name}_sorted.bam
// 		samtools index $samtools_index_options ${name}_sorted.bam
//     	"""	
// }

