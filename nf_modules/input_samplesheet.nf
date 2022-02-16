#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.file_ext = [:]
params.single_end = [:]

include { BAMtoFastq }                  from './files.mod.nf'
include { SRAtoFastq }                  from './files.mod.nf'
include { BISMARK_GENOMEPREPARATION }   from './genomes.mod.nf'

workflow INPUT_SAMPLESHEET {
    take:
    fileList
    outputDir

    main:

    ch_input_genomes = Channel
        .fromPath( fileList )
        .ifEmpty { exit 1, "Provide a tab separated table indicating the read pairs and reference fasta."}
        .splitCsv(sep: '\t')
        .map{ row -> [row[0], file(row[1]).baseName, row[1]] }

    ch_input_files = Channel
        .fromPath( fileList )
        .splitCsv(sep: '\t')
        .map{ row -> [row[0], row[2]] }
    

    prepare_genomes = BISMARK_GENOMEPREPARATION( ch_input_genomes, outputDir )

    if (params.file_ext == "fastq"){
        file_ch = Channel.fromFilePairs( fileList )
    } else if (params.file_ext == "bam"){
        file_ch = BAMtoFastq ( ch_input_files ).fastq_ch
    } else if ( params.file_ext == "sra" ){
        file_ch = SRAtoFastq ( ch_input_files ).fastq_ch
    }

    emit:
    reads       = file_ch
    fasta       = ch_input_genomes
    bismark     = prepare_genomes.index
}