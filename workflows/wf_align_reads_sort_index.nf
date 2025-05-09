
nextflow.enable.dsl=2

include { PREPARE_GENOME }                                      from '../nf_modules/genomes.mod.nf'

include { BWA_MEM }                                             from '../nf_modules/bwa.nf'
include { BOWTIE2 }                                             from '../nf_modules/bowtie2.mod.nf'
include { MINIMAP2_ALIGN }                                      from '../nf_modules/minimap.nf'
include { SAMTOOLS_VIEW; SAMTOOLS_SORT; SAMTOOLS_STATS }        from '../nf_modules/samtools.nf' 
include { PICARD_ADDINFO; PICARD_COLLECTMULTIPLEMETRICS }       from '../nf_modules/picard.nf'
include { MULTIQC }                                             from '../nf_modules/multiqc.mod.nf'

// Align the reads to the reference genome
workflow BWA_ALIGN_SORT_INDEX {

    take:
    ch_input_reads                  // channel: (sample_name, path(fastq))
    val_input_genome                // value: (genome_ref, path('genome_fasta') )
    val_outdir                      // value: val('outdir')
    val_bwa_args                    // value: val('args')
    val_verbose                     // value: boolean

    main:
    // build genome index
    ch_input_genome = PREPARE_GENOME    ( val_input_genome, "${val_outdir}/genome_index", false, false, true )
    // .map{ row -> [file(row).baseName, file(row).baseName, file(row)] }

    // Align reads
    BWA_MEM                             ( ch_input_genome.bwa.collect(), ch_input_reads, "${val_outdir}/bam", val_bwa_args, val_verbose)
    SAMTOOLS_VIEW                       ( BWA_MEM.out.sam, "${val_outdir}/bam")

    // Add read group information and sorting bam file
    PICARD_ADDINFO                      ( SAMTOOLS_VIEW.out.bam, "${val_outdir}/bam")
    SAMTOOLS_SORT                       ( PICARD_ADDINFO.out.mkdup_bam, "${val_outdir}/sorted_bam")

    // Stats
    SAMTOOLS_STATS                      ( SAMTOOLS_SORT.out.bam, "${val_outdir}/aligned_stats", val_verbose)
    PICARD_COLLECTMULTIPLEMETRICS       ( SAMTOOLS_SORT.out.bam, ch_input_genome.fasta.collect(), "${val_outdir}/aligned_stats")
    
    // multiqc
    multiqc_ch = SAMTOOLS_STATS.out.stats.mix(
            PICARD_COLLECTMULTIPLEMETRICS.out.report.ifEmpty([])
        ).collect()
    MULTIQC                 (multiqc_ch, "multiqc_aligned", "${val_outdir}/multiqc", ' ', val_verbose)

    emit:
    ref                         = ch_input_genome.fasta                               //  channel: id, genome_name, path(fasta)
    sorted_bam                  = SAMTOOLS_SORT.out.bam                         //  channel: sample_name, path(bam), path(bam_idx)
    stats_report                = PICARD_COLLECTMULTIPLEMETRICS.out.metrics     //  channel: sample_name, path(stats)
    multiqc_report              = MULTIQC.out.html
    
}

workflow BOWTIE2_ALIGN_SORT_INDEX {

    take:
    ch_input_reads                  // channel: (sample_name, path(fastq))
    val_input_genome                // value: (genome_ref, path('genome_fasta') )
    val_outdir                      // value: val('outdir')
    val_bowtie2_args                    // value: val('args')
    val_verbose                     // value: boolean

    main:
    // build genome index for bowtie2
    ch_input_genome = PREPARE_GENOME    ( val_input_genome, "${val_outdir}/genome_index", false, true, false )

    // Align reads
    BOWTIE2                             (ch_input_genome.bowtie2.collect(), ch_input_reads, "${val_outdir}/bowtie2", val_bowtie2_args, val_verbose)
    SAMTOOLS_VIEW                       ( BOWTIE2.out.sam, "${val_outdir}/bam")

    // Add read group information and sorting bam file
    PICARD_ADDINFO                      ( SAMTOOLS_VIEW.out.bam, "${val_outdir}/bam")
    SAMTOOLS_SORT                       ( PICARD_ADDINFO.out.mkdup_bam, "${val_outdir}/sorted_bam")

    // Stats
    SAMTOOLS_STATS                      ( SAMTOOLS_SORT.out.bam, "${val_outdir}/aligned_stats", val_verbose)
    PICARD_COLLECTMULTIPLEMETRICS       ( SAMTOOLS_SORT.out.bam, ch_input_genome.fasta.collect(), "${val_outdir}/aligned_stats")
    
    // multiqc
    multiqc_ch = SAMTOOLS_STATS.out.stats.mix(
            BOWTIE2.out.stats.ifEmpty([]),
            PICARD_COLLECTMULTIPLEMETRICS.out.report.ifEmpty([])
        ).collect()
    MULTIQC                 (multiqc_ch, "multiqc_aligned", "${val_outdir}/multiqc", ' ', val_verbose)

    emit:
    ref                         = ch_input_genome.fasta                               //  channel: id, genome_name, path(fasta)
    sorted_bam                  = SAMTOOLS_SORT.out.bam                         //  channel: sample_name, path(bam), path(bam_idx)
    stats_report                = PICARD_COLLECTMULTIPLEMETRICS.out.metrics     //  channel: sample_name, path(stats)
    multiqc_report              = MULTIQC.out.html
    
}


workflow MINIMAP2_SORT_INDEX {

    take:
    ch_input_reads                  // channel: (sample_name, path(fastq))
    val_input_genome                // value: (genome_ref, path('genome_fasta') )
    val_outdir                      // value: val('outdir')
    val_minimap2_args               // value: val('args')
    val_verbose                     // value: boolean

    main:
    // build genome index for bowtie2
    ch_input_genome = PREPARE_GENOME    ( val_input_genome, "${val_outdir}/genome_index", false, false, false )

    // Align reads
    MINIMAP2_ALIGN                      (ch_input_genome.fasta.collect(), ch_input_reads, "${val_outdir}/minimap2", val_minimap2_args)
    SAMTOOLS_VIEW                       ( MINIMAP2_ALIGN.out.sam, "${val_outdir}/minimap2")

    // Add read group information and sorting bam file
    PICARD_ADDINFO                      ( SAMTOOLS_VIEW.out.bam, "${val_outdir}/minimap2")
    SAMTOOLS_SORT                       ( PICARD_ADDINFO.out.mkdup_bam, "${val_outdir}/sorted_bam")

    // Stats
    SAMTOOLS_STATS                      ( SAMTOOLS_SORT.out.bam, "${val_outdir}/aligned_stats", val_verbose)
    PICARD_COLLECTMULTIPLEMETRICS       ( SAMTOOLS_SORT.out.bam, ch_input_genome.fasta.collect(), "${val_outdir}/aligned_stats")
    
    // multiqc
    multiqc_ch = SAMTOOLS_STATS.out.stats.mix(
            MINIMAP2_ALIGN.out.log.ifEmpty([]),
            PICARD_COLLECTMULTIPLEMETRICS.out.report.ifEmpty([])
        ).collect()
    MULTIQC                 (multiqc_ch, "multiqc_aligned", "${val_outdir}/multiqc", ' ', val_verbose)

    emit:
    ref                         = ch_input_genome.fasta                               //  channel: id, genome_name, path(fasta)
    sorted_bam                  = SAMTOOLS_SORT.out.bam                         //  channel: sample_name, path(bam), path(bam_idx)
    stats_report                = PICARD_COLLECTMULTIPLEMETRICS.out.metrics     //  channel: sample_name, path(stats)
    multiqc_report              = MULTIQC.out.html

}