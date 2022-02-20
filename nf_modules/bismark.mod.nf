nextflow.enable.dsl=2

// parameters passed in by specialised pipelines
params.singlecell = ''
params.pbat = false
params.unmapped = false
params.read_identity = ''
params.save_intermediate = false

process BISMARK {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.
	label 'process_high'

	conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.23.0--0' :
        'quay.io/biocontainers/bismark:0.23.0--0' }"

	
		
    input:
	    tuple val(name), path(reads), val(genome_name), path (genome)
		val (outputdir)
		val (bismark_args)
		val (verbose)

	output:
	    tuple val(name), path ("${bismark_name}*bam"),       					emit: bam
		path "${bismark_name}*_report.txt",										emit: report
		// we always pass back the original name so we can use .join() later on, e.g. for bismark2bedGraph
		tuple val(name), path ("*unmapped_reads_1.fq.gz"), optional: true,		emit: unmapped1
		tuple val(name), path ("*unmapped_reads_2.fq.gz"), optional: true,		emit: unmapped2

	
	publishDir "$outputdir", mode: "copy", overwrite: true,
		saveAs: {filename ->
			if( filename.indexOf("report.txt") > 0 ) "reports/$filename"
			// else if( filename.indexOf("fq.gz" ) > 0) "unmapped/$filename"
			else if( params.save_intermediate ) filename
			else null
		}


    script:
		readString = ""

		if (verbose){
			println ("[MODULE] BISMARK ARGS: " + bismark_args)
		}

		// Options we add are
		bismark_options = bismark_args
		if (params.singlecell){
			bismark_options += " --non_directional "
		}
		else{
		
		}
		
		unmapped_1_name = ''
		unmapped_2_name = ''
		
		if (params.unmapped){
			bismark_options += " --unmapped "
			unmapped_1_name = name + "_unmapped_R1"
			unmapped_2_name = name + "_unmapped_R2"
		}

		if (params.pbat){
			bismark_options += " --pbat "
		}

		if (reads instanceof List) {
			readString = "-1 "+reads[0]+" -2 "+reads[1]
		}
		else {
			readString = reads
		}


		unmapped_name = ''	
			// add Genome build and aligner to output name
		if (params.read_identity == "1" || params.read_identity == "2"){
			// println ("FILENAME IS: $reads")
			if (params.read_identity == "1"){
				unmapped_name = name + "_unmapped_R1"
			}
			else{
				unmapped_name = name + "_unmapped_R2"
			}

			if (bismark_args =~ /-hisat/){ // if HISAT2 was given on the command line
				bismark_name = unmapped_name + "_" + genome_name + "_bismark_hisat2"
			}
			else{ // default is Bowtie 2
				bismark_name = unmapped_name + "_" + genome_name + "_bismark_bt2"
			}
		}
		else{
			if (bismark_args =~ /-hisat/){ // if HISAT2 was given on the command line
				bismark_name = name + "_" + genome_name + "_bismark_hisat2"
			}
			else{ // default is Bowtie 2
				bismark_name = name + "_" + genome_name + "_bismark_bt2"
			}
		}	
		// println ("Output basename: $bismark_name")
		// --basename $bismark_name
		"""
		bismark  \
			--parallel $task.cpus \
			--genome $genome \
			$bismark_options $readString
		[ -f *_pe.bam ] && mv *_pe.bam ${bismark_name}_pe.bam
		[ -f *_se.bam ] && mv *_se.bam ${bismark_name}_se.bam
		[ -f *_PE_report.txt ] && mv *_PE_report.txt ${bismark_name}_PE_report.txt
		[ -f *_SE_report.txt ] && mv *_SE_report.txt ${bismark_name}_SE_report.bam
		"""
}