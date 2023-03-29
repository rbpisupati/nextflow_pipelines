#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.genome_id = false
params.star_genome_args = ''


workflow PREPARE_GENOME {

    take:
    fasta_file
    outputDir

    main:

    ch_input_genome = Channel
        .fromPath( fasta_file )
        .map{ row -> [file(row).baseName, file(row).baseName, file(row)] }

    // ch_fasta_dir = ch_fasta.parent

    prepare_genome = BISMARK_GENOMEPREPARATION ( ch_input_genome, outputDir )

    emit:
    fasta            = ch_input_genome            //    path: genome.fasta
    // chrom_sizes      = ch_chrom_sizes      //    path: genome.sizes
    bismark          = prepare_genome.index
    // rsem_index       = ch_rsem_index       //    path: rsem/index/
    // hisat2_index     = ch_hisat2_index     //    path: hisat2/index/
    // salmon_index     = ch_salmon_index     //    path: salmon/index/

}

workflow PREPARE_GENOME_RNA {

    take:
    fasta_file
    gtf_file
    outputDir

    main:

    ch_input_genome = Channel
        .fromPath( fasta_file )
        .map{ row -> [params.genome_id ? "$params.genome_id" : file(row).baseName, file(row)] }

        

    ch_gtf_file = Channel.fromPath( gtf_file )
    // ch_fasta_dir = ch_fasta.parent

    prepare_genome = STAR_GENOMEPREPARATION ( ch_input_genome, ch_gtf_file, outputDir )

    emit:
    fasta       = ch_input_genome            //    path: genome.fasta
    gtf         = ch_gtf_file
    // chrom_sizes      = ch_chrom_sizes      //    path: genome.sizes
    star        = prepare_genome.index
    // rsem_index       = ch_rsem_index       //    path: rsem/index/
    // hisat2_index     = ch_hisat2_index     //    path: hisat2/index/
    // salmon_index     = ch_salmon_index     //    path: salmon/index/

}


process STAR_GENOMEPREPARATION {

    tag "$genome_id"
    label 'process_medium'
    storeDir "${outdir}"

    conda "bioconda::star=2.6.1d bioconda::samtools=1.10 conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"
        
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/star%3A2.6.1d--0' :
    //     'quay.io/biocontainers/star:2.6.1d' }"

    input:
        tuple val(genome_id), path(fasta)
        path(gtf)
        val(outdir)

    output:
        tuple val(genome_id), path("$genome_id"), emit: index

    script:
    star_genome_args = params.star_genome_args
    star_genome_args = star_genome_args + ''
    def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''

    """
    samtools faidx $fasta
    NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fasta}.fai`

    mkdir $genome_id

    STAR --runThreadN $task.cpus \
        --runMode genomeGenerate \
        --genomeDir $genome_id \
        --genomeFastaFiles $fasta \
        --sjdbGTFfile $gtf \
        --genomeSAindexNbases \$NUM_BASES \
        $memory \
    """

}

// name	ARS-UCD1.2
// species	Bos_taurus
// fasta	/bi/scratch/Genomes/Bos_taurus/ARS-UCD1.2/
// bismark	/bi/scratch/Genomes/Bos_taurus/ARS-UCD1.2/
// bowtie	/bi/scratch/Genomes/Bos_taurus/ARS-UCD1.2/ARS-UCD1.2
// bowtie2	/bi/scratch/Genomes/Bos_taurus/ARS-UCD1.2/ARS-UCD1.2
// hisat2	/bi/scratch/Genomes/Bos_taurus/ARS-UCD1.2/ARS-UCD1.2
// gtf	/bi/scratch/Genomes/Bos_taurus/ARS-UCD1.2/Bos_taurus.ARS-UCD1.2.98.gtf
// hisat2_splices	/bi/scratch/Genomes/Bos_taurus/ARS-UCD1.2/Bos_taurus.ARS-UCD1.2.98.hisat2_splices.txt

process BISMARK_GENOMEPREPARATION {
    tag "$name"
    label 'process_medium'
    storeDir "${outdir}"
        // mode: params.publish_dir_mode,
        // saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bismark:0.23.0--0"
    } else {
        container "quay.io/biocontainers/bismark:0.23.0--0"
    }

    input:
    tuple val(name), val(genome_id), path(fasta, stageAs: "BismarkIndex/*")
    val(outdir)

    output:
    tuple val(name), val(genome_id), path("$genome_id"), emit: index

    script:
    """
    bismark_genome_preparation BismarkIndex
    mv BismarkIndex/ $genome_id
    """
    // $options.args \\
    // echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//' > ${software}.version.txt
}

def getGenome(name) {

    // Find a file with the same name as the genome in our genomes.d directory

    scriptDir = workflow.projectDir
    
    // // die gracefully if the user specified an incorrect genome
    // def fileName = scriptDir.toString() + "/genomes.d/" + name + ".genome"
    // def testFile = new File(fileName)
    // if (!testFile.exists()) {
    //     println("\nFile >>$fileName<< does not exist. Listing available genomes...\n")
    //     listGenomes()
    // }   
    // else { 
    //     // println ("File $fileName exists.")
    // }

    // genomeFH = new File (fileName).newInputStream()

    // genomeValues = [:]  // initialising map. name is also part of each .genome file

    // genomeFH.eachLine {
    //     sections =  it.split("\\s+",2)
    //     genomeValues[sections[0]] = sections[1]
    // }

    // migrated to igenomes..
    return params.genomes[ params.genome ]

}


def listGenomes(){
    
    println ("These genomes are currently available to choose from:")
    println ("=====================================================")
    scriptDir = workflow.projectDir + "/genomes.d/"
    // println (scriptDir) // last slash is consumed
    allFiles = scriptDir.list()
    
    for( def file : allFiles.sort() ) {
        
        if( file =~ /.genome$/){

            genomeFH = new File(scriptDir.toString() + "/$file").newInputStream()
            name = file.replaceFirst(/.genome/, "")
        
            println (name)
            genomeFH.eachLine {
                if (params.verbose){
                    println ("\t$it")
                }
            }
        }
    }
    println ("\nTo see this list of available genomes with more detailed information about paths and indexes,\nplease re-run the command including '--list_genomes --verbose'\n\n")

    System.exit(1)
}

