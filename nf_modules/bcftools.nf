nextflow.enable.dsl=2

params.bcftools_snpcall_args = [:]
params.bcftools_mpileup_args = [:]
params.bcftools_filter_args = [:]


params.min_base_quality = 30
params.min_snp_qual = 30


process SNPCALL_BCFTOOLS {

    tag "${name}"
    label 'bigMem'
	label 'quadCore'

    conda (params.enable_conda ? "bioconda::bcftools=1.18-0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools%3A1.18--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.18--h8b25389_0' }"


    publishDir "${outdir}", mode: "copy"

    input:
    tuple val(genome_id), path(index)
    tuple val(name), file(bam)
    val(outdir)
    
    output:
    tuple val(name), path("${name}.vcf.gz"), path("${name}.vcf.gz.tbi"), emit: vcf

    script:
    
    mpileup_args = params.bcftools_mpileup_args
    mpileup_args = (mpileup_args.contains('-a')) ? mpileup_args : mpileup_args + ' -a FORMAT/DP '
    mpileup_args = (mpileup_args.contains('-Q')) ? mpileup_args : mpileup_args + " -Q $params.min_base_quality "

    snpcall_args = params.bcftools_snpcall_args
    snpcall_args = (snpcall_args.contains('-C')) ? snpcall_args : snpcall_args + ' -c '
    
    filter_args = params.bcftools_filter_args
    filter_args = (filter_args.contains('-e')) ? filter_args : filter_args + " -e 'QUAL < $params.min_snp_qual' "

    // flagW=99,147 flagC=83,163 // set for paired end
    // known_sites_mpile_cmd = params.known_sites != false ? "-T ${params.known_sites}" : ''
    // known_sites_call_cmd = params.known_sites != false ? "-C alleles -m -T ${params.known_sites}" : '-m'
    """
    bcftools mpileup --threads ${task.cpus}\
    $mpileup_args -f ${genome_id}/genome.fa $bam | \
    bcftools call --threads ${task.cpus} \
    $snpcall_args \
    -O z -o ${name}.raw.vcf.gz

    bcftools filter $filter_args \
    -Oz ${name}.raw.vcf.gz | \
    bcftools view -V indels > ${name}.vcf


    bgzip ${name}.vcf && tabix ${name}.vcf.gz

    """
}