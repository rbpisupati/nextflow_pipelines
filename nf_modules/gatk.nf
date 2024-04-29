nextflow.enable.dsl=2

params.gatk_snpcall_args = [:]


process GATK4_HAPLOCALLER {
    tag "$name"
    label 'bigMem'
	label 'quadCore'

    conda (params.enable_conda ? 'bioconda::gatk4=4.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4%3A4.2.0.0--0':
        'quay.io/biocontainers/gatk4:4.2.0.0--0' }"


    publishDir "${outdir}", mode: 'copy'
    
    input:
    tuple val(genome_id), path(index)
    tuple val(name), path(bam)
    val(outdir)

    output:
    tuple val(name), path("${name}.g.vcf.gz"), emit: vcf
    // tuple val(name), path("${name}.metrics.txt"), emit: metrics

  
    script:
    gatk_options = params.gatk_snpcall_args

    gatk_options = (gatk_options.contains('ERC')) ? gatk_options : gatk_options + ' -ERC GVCF '
    // --tmp-dir ${params.tmpdir}\
    """
    gatk  HaplotypeCaller $gatk_options \
        -R $genome_id/genome.fa \
        -I $bam -O ${name}.g.vcf.gz
    """
}
