nextflow.enable.dsl=2

params.umeth = ""

process METHYLPY_ALLC {
    tag "$name"

    conda ( "/users/rahul.pisupati/.conda/envs/nf-core-methylseq-1.5/" )
    // conda (params.enable_conda ? "bioconda::methylpy=1.4.3" : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/methylpy:1.4.3--py37h41a55b7_0' :
    //     'quay.io/biocontainers/methylpy:1.4.3--py37h41a55b7_0' }"

    publishDir "${outputdir}", mode: 'copy',
        saveAs: {filename ->
            if (filename =~ '^allc' ) "methylpy/$filename"
            else if (filename =~ '^conversion' ) "info/$filename"
            else if (filename =~ '^log' ) "info/log.${name}.txt"
        }

    input:
    tuple val(name), path(bam)
    path (genome)
    val (outputdir)

    output:
    tuple val(name), path ("allc_*"),       emit: allc
    path("log.txt"),                        emit: methylpy_log
    path ("conversion_rate_*"),             emit: conv_rate
    

    script:
    """
    samtools sort $bam -o ${name}_sorted.bam
	samtools index ${name}_sorted.bam
    methylpy call-methylation-state \
    --input-file ${name}_sorted.bam  \
    --paired-end True \
    --sample $name \
    --ref-fasta $genome \
    --unmethylated-control $params.umeth \
    --num-procs ${task.cpus} > log.txt 2>&1
    cat log.txt | grep "non-conversion rate" > conversion_rate_${name}.txt
    """
}