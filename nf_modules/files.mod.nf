#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.file_ext = [:]
params.single_end = [:]

include { BISMARK_GENOMEPREPARATION }   from './genomes.mod.nf'
include { BOWTIE2_BUILD }   from './genomes.mod.nf'

workflow INPUT_FILES {
    take:
    fileList

    main:
    if (params.file_ext == "fastq"){
        if (params.single_end){
            file_ch = Channel
                        .fromPath( fileList )
                        .map{ row ->  [ file(row).baseName.replace(".fastq", "").replace(".fq", ""), file(row, checkIfExists: true) ]  }

        } else {
            file_ch = Channel.fromFilePairs( fileList )
        }
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

    // module = ['build-env/2020', 'sra-toolkit/2.9.6-1-centos_linux64']
    conda (params.enable_conda ? 'bioconda::sra-tools=2.11.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0' :
        'quay.io/biocontainers/sra-tools:2.11.0--pl5262h314213e_0' }"

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


workflow INPUT_SAMPLESHEET {
    take:
    fileList
    outputDir

    main:

    // ch_input_genomes = Channel
    //     .fromPath( fileList )
    //     .ifEmpty { exit 1, "Provide a tab separated table indicating the read pairs and reference fasta."}
    //     .splitCsv(sep: ',')
    //     .map{ row -> [row[0], row[1], file(row[2])] }
    //     // sample_id, reference_id, reference fasta

    ch_input_files = Channel
        .fromPath( fileList )
        .splitCsv(sep: ',')
        .map{ row -> [row[0], file(row[1]), file(row[2])] }
        .groupTuple(by: 0)
        // sample_id, read_1, read_2

    
    // need to merge files with same sample ID
    ch_input_files_merged = mergeInputFastq( ch_input_files.filter{ it[1].size() > 1 }, outputDir )

    ch_input_files_final = ch_input_files_merged.concat( ch_input_files.filter{ it[1].size() == 1 }.map{ it -> [it[0], it[1][0], it[2][0] ] } )
        .map{ it -> [ it[0], [it[1], it[2]] ] }
    
    // Prepare genomes
    // bowtie_indices = BOWTIE2_BUILD (ch_input_genomes, outputDir)
    // bismark_indices = BISMARK_GENOMEPREPARATION( ch_input_genomes, outputDir )

    emit:
    reads       = ch_input_files_final
    // fasta       = ch_input_genomes
    // bowtie2     = bowtie_indices
}

process mergeInputFastq {
    tag { "${sample_id}" }

    publishDir "$outputDir", mode: 'copy'
    // storeDir "$params.outdir"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(read_1), path(read_2)
    val(outputDir)

    output:
    tuple val(sample_id), path("${sample_id}*_1.fq.gz"), path("${sample_id}*_2.fq.gz")

    script:
    // def file_ext = ${read_1[0]}.getExtension()
    // if( file_ext == 'gz')
    """
    cat ${read_1[0]} ${read_1[1]} > ${sample_id}_merged_1.fq.gz
    cat ${read_2[0]} ${read_2[1]} > ${sample_id}_merged_2.fq.gz
    """
    // gzip -c ${sample_id}_1.fq > ${sample_id}_1.fq.gz
    // else if(file_ext == 'fq' || file_ext == 'fastq')
    // """
    // cat ${read_1[0]} ${read_1[1]} > ${sample_id}_1.fq
    // gzip -c ${sample_id}_1.fq > ${sample_id}_1.fq.gz
    // cat ${read_2[0]} ${read_2[1]} > ${sample_id}_2.fq
    // gzip -c ${sample_id}_2.fq > ${sample_id}_2.fq.gz
    // """
    // else if(file_ext == 'bam')
    // """
    // samtools merge -n --threads ${task.cpus} ${sample_id}.bam ${read_1[0]} ${read_1[1]}
    // """
    // else
    // error "Invalid file extensions"
    
}