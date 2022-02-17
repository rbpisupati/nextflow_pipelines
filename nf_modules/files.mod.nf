#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.file_ext = [:]
params.single_end = [:]

workflow INPUT_FILES {
    take:
    fileList

    main:
    if (params.file_ext == "fastq"){
        file_ch = Channel.fromFilePairs( fileList )
    } else if (params.file_ext == "bam"){
        ch_input_bams = Channel
                            .fromPath( fileList )
                            .map{ row -> [ file(row).baseName, file(row, checkIfExists: true) ] }
        file_ch = BAMtoFastq ( ch_input_bams ).fastq_ch
    } else if ( params.file_ext == "sra" ){
        ch_input_sra = Channel
                            .fromPath( fileList )
                            .map{ row -> [ file(row).baseName, file(row, checkIfExists: true) ] }
        file_ch = SRAtoFastq ( ch_input_sra ).fastq_ch
    }

    emit:
    file_ch
}



process BAMtoFastq {
    tag "${name}"
    // storeDir "${workflow.workDir}/rawreads"

    conda (params.enable_conda ? 'bioconda::picard=2.25.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.25.7--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.25.7--hdfd78af_0' }"

    input:
    tuple val(name), path(input_reads)

    output:
    tuple val(name), path(output), emit: fastq_ch

    script:
    def c_single_end = params.single_end ? "FASTQ=${name}.fastq.gz" : "FASTQ=${name}_1.fastq.gz SECOND_END_FASTQ=${name}_2.fastq.gz"
    output = params.single_end ? '*.fastq.gz' : '*_{1,2}.fastq.gz'
    """
    picard SamToFastq I=$input_reads $c_single_end VALIDATION_STRINGENCY=LENIENT
    """
}

process SRAtoFastq {
    tag "${name}"
    // storeDir "${workflow.workDir}/rawreads"

    module = ['build-env/2020', 'sra-toolkit/2.9.6-1-centos_linux64']
    // conda (params.enable_conda ? 'bioconda::sra-tools=2.11.0' : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0' :
    //     'quay.io/biocontainers/sra-tools:2.11.0--pl5262h314213e_0' }"

    input:
    tuple val(name), path(input_reads)

    output:
    tuple val(name), path(output), emit: fastq_ch

    script:
    def c_single_end = params.single_end ? "" : "--split-files"
    output = params.single_end ? '*.fastq.gz' : '*_{1,2}.fastq.gz'
    """
    fastq-dump --gzip $c_single_end $input_reads
    """
}