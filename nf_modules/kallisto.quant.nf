nextflow.enable.dsl=2


params.args_kallisto = false
params.args_kallisto_prep = false
params.single_end = [:]


workflow KALLISTO {
    take:
    ch_fasta_file
    ch_transcript_file
    ch_gtf_file
    ch_raw_reads
    outputDir

    main:
    prepare_idx = KALLISTO_INDEX( ch_transcript_file, "$outputDir/index_transcript" )
    
    // ch_input_files_to_map = TRIM_GALORE.out.reads.combine( genome.star )
    kallisto_out = KALLISTO_QUANT(ch_raw_reads, prepare_idx.index.collect(), ch_gtf_file.collect(), "$outputDir"   )


    emit:
    index       = prepare_idx.index   //
    quant       = kallisto_out.abundance_h5

}


process KALLISTO_INDEX {
    tag "$fasta"
    label 'process_medium'
    storeDir "${outputDir}"

    conda "bioconda::kallisto=0.48.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.48.0--h15996b6_2':
        'quay.io/biocontainers/kallisto:0.48.0--h15996b6_2' }"
    
    input:
    path fasta
    val(outputDir)

    output:
    path "kallisto" , emit: index

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = params.args_kallisto_prep ?: ''
    """
    kallisto \\
        index \\
        $args \\
        -i kallisto \\
        $fasta
    """
}


process KALLISTO_QUANT {
    tag "$name"
    label 'process_medium'

    conda "bioconda::kallisto=0.48.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.48.0--h15996b6_2':
        'quay.io/biocontainers/kallisto:0.48.0--h15996b6_2' }"

    input:
    tuple val(name), path(reads)
    path index
    path gtf
    val(outputDir)
    // path chromosomes

    output:
    tuple val(name), path("abundance.tsv"), emit: abundance
    tuple val(name), path("abundance.h5") , emit: abundance_h5
    tuple val(name), path("run_info.json"), emit: run_info
    tuple val(name), path("*.log.txt")    , emit: log

    publishDir "$outputDir", mode: "copy", overwrite: true,
        saveAs: {filename ->
			if( filename.indexOf("tsv") > 0 ) "counts/${name}.abundance.tsv"
			else if( filename.indexOf("h5" ) > 0) "counts_h5/${name}.abundance.h5"
            else if( filename.indexOf("json" ) > 0) "run_info/${name}.run_info.json"
            else if( filename.indexOf("txt" ) > 0) "log/${filename}"
			else null
        }


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = params.args_kallisto ?: ''
    def gtf_input = gtf ? "--gtf ${gtf}" : ''
    
    args = (args.contains('--bootstrap')) ? args : args + ' --bootstrap-samples 100 '

    if (params.single_end){
        args = args + ' --single --fragment-length 500 --sd 10 '
    }
    // def single = meta.single_end ? "--single  ${task.ext.fragment_len} --sd ${task.ext.sd}" : ""        
    
    // def chromosomes_input = chromosomes ? "--chromosomes ${chromosomes}" : ''
    // ${chromosomes_input} \\
    """
    kallisto quant \\
        --threads ${task.cpus} \\
        --index ${index} \\
        ${gtf_input} \\
        ${args} \\
        -o . \\
        ${reads} 2> >(tee -a ${name}.log.txt >&2)

    """

}