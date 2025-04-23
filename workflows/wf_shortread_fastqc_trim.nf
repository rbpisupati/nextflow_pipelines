
nextflow.enable.dsl=2

include { INPUT_FILES }         from '../nf_modules/files.mod.nf'
include { FASTQC }              from '../nf_modules/fastqc.mod.nf'
include { FASTQC as FASTQC2 }   from '../nf_modules/fastqc.mod.nf'
include { TRIM_GALORE }         from '../nf_modules/trim_galore.mod.nf'
include { MULTIQC }             from '../nf_modules/multiqc.mod.nf'

// Perform the required fastqc and trimming on short reads

workflow FASTQC_TRIM_MULTIQC {

    take:
    ch_input_files                 // channel: (sample_name, path(fastq))
    val_outdir                      // value: val('outdir')
    val_fastqc_args                 // value: val('args')
    val_trim_galore_args            // value: val('args')
    val_multiqc_args                // value: val('args')
    val_verbose                     // value: boolean

    main:
    
    // fastqc
    FASTQC                  (ch_input_files, "${val_outdir}/fastqc", val_fastqc_args, val_verbose)

    // trim_galore
    TRIM_GALORE             (ch_input_files, "${val_outdir}/trim_galore", val_trim_galore_args, val_verbose)

    // fastqc after trimming
    FASTQC2                 (TRIM_GALORE.out.reads, "${val_outdir}/fastqc_after_trim", val_fastqc_args, val_verbose)


    // multiqc
    multiqc_ch = FASTQC.out.report.mix(
            TRIM_GALORE.out.report,
            // FASTP.out.report,
            // FASTQ_SCREEN.out.report.ifEmpty([]),
            FASTQC2.out.report.ifEmpty([])
        ).collect()
    MULTIQC                 (multiqc_ch, "multiqc_rawdata", "${val_outdir}/multiqc", val_multiqc_args, val_verbose)


    emit:
    input                       =   ch_input_files
    trim_reads                  =   TRIM_GALORE.out.reads
    fastqc_report               =   FASTQC.out.report
    trim_report                 =   TRIM_GALORE.out.report
    fastqc_report_trimmed       =   FASTQC2.out.report
    multiqc_report              =   MULTIQC.out.html
    
}
