nextflow.enable.dsl=2

include { INPUT_FILES }             from '../nf_modules/files.mod.nf'
include { PORECHOP }                from '../nf_modules/longread_qc.nf'
include { PYCOQC }                  from '../nf_modules/longread_qc.nf'
include { NANOPLOT }                from '../nf_modules/longread_qc.nf'

include { MULTIQC }                 from '../nf_modules/multiqc.mod.nf'


// Perform the required fastqc and trimming on short reads

workflow PORECHOP_QC {

    take:
    ch_longreads                    // channel: (sample_name, path(fastq))
    val_outdir                      // value: val('outdir')

    main:
    
    // Trip long reads using Porechop
    PORECHOP                    (ch_longreads, "${val_outdir}/porechop")
    
    // Perform quality control on long reads using PycoQC for FAST5 files
    // PYCOQC                     (PORECHOP.out.reads, "${val_outdir}/pycoqc")

    // QC using Nanoplot
    NANOPLOT                    (PORECHOP.out.reads, "${val_outdir}/nanoplot")

    // Generate a MultiQC report 
    multiqc_ch = PORECHOP.out.log.mix(
            NANOPLOT.out.report,
        ).collect()
    MULTIQC                     (multiqc_ch, "multiqc_rawdata", "${val_outdir}/multiqc", ' ', false)
    

    emit:
    input                           =   ch_longreads
    porechop_reads                  =   PORECHOP.out.reads
    
}
