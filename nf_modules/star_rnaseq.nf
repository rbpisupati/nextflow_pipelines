nextflow.enable.dsl=2

params.seq_platform='illumina'
params.seq_center='ncbi'



process STAR_ALIGN {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	label 'bigMem'
	label 'multiCore'

	conda "bioconda::star=2.6.1d bioconda::samtools=1.10 conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

    input:
        tuple val(name), path(reads), val(genome_name), path (genome)
        val (outputdir)

	output:
        tuple val(prefix), path('*d.out.bam')         , emit: bam
        tuple val(prefix), path('*Log.final.out')     , emit: log_final
        tuple val(prefix), path('*Log.out')           , emit: log_out
        tuple val(prefix), path('*Log.progress.out')  , emit: log_progress

        tuple val(prefix), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
        tuple val(prefix), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
        tuple val(prefix), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
        tuple val(prefix), path('*.tab')                   , optional:true, emit: tab
        tuple val(prefix), path('*.out.junction')          , optional:true, emit: junction
        tuple val(prefix), path('*.out.sam')               , optional:true, emit: sam


	publishDir "$outputdir", mode: "copy", overwrite: true,
		saveAs: {filename ->
			if( filename.indexOf("bam") > 0 ) "aligned_bam/$filename"
			else if( filename.indexOf("out" ) > 0) "log_files/$filename"
            else if( filename.indexOf("tab" ) > 0) "tab/$filename"
			else null
        }

    script:
    star_args = params.star_args
    star_args = (star_args.contains('--outSAMtype')) ? star_args : star_args + ' --outSAMtype BAM SortedByCoordinate '
    star_args = star_args + "--outSAMattrRGline ID:$name 'CN:$params.seq_center' 'SM:$name' $params.seq_platform "
    prefix = "${name}.${genome_name}"
    """
    STAR \\
        --genomeDir $genome \\
        --readFilesIn $reads  \\
        --readFilesCommand zcat \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix ${prefix}. \\
        $star_args
    """
}