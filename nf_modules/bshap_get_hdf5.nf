nextflow.enable.dsl=2

params.umeth = false // "ChrC"

process BISMARK_TO_HDF5 {
    tag "$name"
    label 'bigMem'

    conda = "$HOME/.conda/envs/py3_quant"
    // conda ( "/users/rahul.pisupati/.conda/envs/nf-core-methylseq-1.5/" )
    // conda (params.enable_conda ? "bioconda::methylpy=1.4.3" : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/methylpy:1.4.3--py37h41a55b7_0' :
    //     'quay.io/biocontainers/methylpy:1.4.3--py37h41a55b7_0' }"

    publishDir "${outputdir}", mode: 'copy',
        saveAs: {filename ->
            if( filename.indexOf("hdf5") > 0 ) "hdf5/$filename"
            else if( filename.indexOf("conv_rate.txt") > 0 ) "conv_rates/$filename"
            // else if (filename =~ '^log' ) "info/log.${name}.txt"
            else null
        }

    input:
    tuple val(name), path(cx_report), val(genome_name), path(genome)
    val (outputdir)

    output:
    tuple val(name), path ("*hdf5"),        emit: allc_hdf5
    path("*conv_rate.txt"), optional: true, emit: conv_rate
    

    script:
    """
    bshap allc_to_hdf5 \
        --input_file $cx_report \
        -t "bismark" \
        -f $genome --umeth_control $params.umeth \
        --output ${name}_${genome_name}.bismark

    """
}