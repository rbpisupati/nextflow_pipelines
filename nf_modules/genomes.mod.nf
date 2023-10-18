#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.genome_id = false
params.star_genome_args = ''
params.bowtie2_build_args = ''


workflow PREPARE_GENOME {

    take:
    fasta_file
    outputDir
    bismark
    bowtie2

    main:

    ch_input_genome = Channel
        .fromPath( fasta_file )
        .map{ row -> [file(row).baseName, file(row).baseName, file(row)] }

    // Generate this index only if given
    if( bowtie2 ){
        bowtie2_index = BOWTIE2_BUILD ( ch_input_genome, outputDir)
    } else {
        bowtie2_index = false
    }

    // ch_fasta_dir = ch_fasta.parent
    if( bismark ){
        bismark_index = BISMARK_GENOMEPREPARATION ( ch_input_genome, outputDir )
    } else {
        bismark_index = false
    }

    emit:
    fasta            = ch_input_genome            //    path: genome.fasta
    // chrom_sizes      = ch_chrom_sizes      //    path: genome.sizes
    bismark          = bismark_index
    bowtie2           = bowtie2_index
    // rsem_index       = ch_rsem_index       //    path: rsem/index/
    // hisat2_index     = ch_hisat2_index     //    path: hisat2/index/
    // salmon_index     = ch_salmon_index     //    path: salmon/index/

}

workflow PREPARE_GENOME_RNA {

    take:
    fasta_file
    gtf_file
    transcript_fasta
    outputDir

    main:

    ch_input_genome = Channel
        .fromPath( fasta_file )
        .map{ row -> [params.genome_id ? "$params.genome_id" : file(row).baseName, file(row)] }

    ch_gtf_file = Channel.fromPath( gtf_file )
    ch_transcript = Channel.fromPath( transcript_fasta )
    // ch_fasta_dir = ch_fasta.parent

    emit:
    fasta       = ch_input_genome            //    path: genome.fasta
    gtf         = ch_gtf_file
    transcript  = ch_transcript
    // chrom_sizes      = ch_chrom_sizes      //    path: genome.sizes
    // star        = prepare_genome.index
    // rsem_index       = ch_rsem_index       //    path: rsem/index/
    // hisat2_index     = ch_hisat2_index     //    path: hisat2/index/
    // salmon_index     = ch_salmon_index     //    path: salmon/index/

}

process BOWTIE2_BUILD {
    tag "$name"
    label 'process_medium'

    storeDir "${outdir}"
    
    conda "bioconda::bowtie2=2.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0' :
        'biocontainers/bowtie2:2.4.4--py39hbb4e92a_0' }"

    input:
    tuple val(name), val(genome_id), path(fasta, stageAs: "bowtie2/*")
    val(outdir)

    output:
    tuple val(name), path("$genome_id")

    script:
    """
    bowtie2-build $params.bowtie2_build_args --threads $task.cpus $fasta bowtie2/${genome_id}
    mv bowtie2 $genome_id
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

process DOWNLOAD_GENOMES{
    tag "$common_name"
    label 'process_medium'


    publishDir "${outdir}", mode: 'move'
        // saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-biomartr:1.0.4--r42h3342da4_0"
    } else {
        container "quay.io/biocontainers/r-biomartr:1.0.4--r42h3342da4_0"
    }

    input:
    tuple val(common_name), val(species_name)
    val(outdir)
    val(assembly)
    val(annotation)

    output:
    tuple val(common_name), path("$common_name"), emit: genome

    script:
    """
    #!/usr/bin/env Rscript
    library(biomartr)

    if($assembly){
        genome.refseq <- getGenome(
                        db = "refseq", 
                        organism = "$species_name", 
                        path = file.path("$common_name", 'assembly')
                    )

        file.rename(genome.refseq, file.path("$common_name", 'assembly','genome.fa.gz'))
    }
    
    if ($annotation){
        cds.refseq = getCDS(
                        db = "refseq", 
                        organism = "$species_name",
                        path = file.path("$common_name", 'cds')
        )
        file.rename(cds.refseq, file.path("$common_name",'cds','cds.fa.gz'))

        gff.refseq = getGFF(
                        organism = "$species_name",
                        path = file.path("$common_name", 'annotation'),
                        remove_annotation_outliers = TRUE
                    )

        file.rename(gff.refseq, file.path("$common_name", 'annotation','gff.gz'))
    }
    """
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

