#!/usr/bin/env nextflow

nextflow.enable.dsl=2 


/* 
 * Proof of concept Nextflow based trim(v3v4 primer and nextera backbone), QC, filter before Qimme2 import and process
 * 
 */ 

 
/*
 * Defines some parameters in order to specify 
 * read pairs by using the command line options
 */
//params.reads = "/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/hammer/reads/*R{1,2}.fastq.gz"
params.reads = "xxx"
params.fqsConfig = "/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/16s/fastq_screen.conf"
params.outdir = 'results'
log.info """
         Q2 prescreen  P I P E L I N E:16S    
         =============================
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
         .stripIndent()

 
/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */
Channel 
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .view() 
 
 

Channel                                                                         
    .fromFilePairs( params.reads )                                              
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }        
    .set { read_pairs } 



process fastqc {
    input:
    tuple val(pair_id), path(read) 
    output:
    path("fastqc_${pair_id}_logs") 
 
    script:

    """
    mkdir fastqc_${pair_id}_logs
    module load fastqc
    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${pair_id}.R1.fastq.gz
    """
}
process fqscreenQc {
    input:
    tuple val(pair_id), path(read) 
    output:
     path("${pair_id}.R1_screen.txt") 

    script:
 
    """
    module load bowtie2
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app
    \$APP/fastq_screen/v0.14.0/fastq_screen --conf ${params.fqsConfig} --threads ${task.cpus} ${pair_id}.R1.fastq.gz


    """
}

process fqscreenFilt {
    publishDir "${params.outdir}/fqscreenFilt", mode: 'copy', pattern: '*tagged_filter.fastq.gz'
    input:
    tuple val(pair_id), path(read)
    output:
    path("${pair_id}*tagged_filter.fastq.gz")

    script:

    """
    module load bowtie2
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app
    #filter mitochondria reads and mouse seq. hold on filtering plastid reads due to db congtamination of ncbi refseq
    \$APP/fastq_screen/v0.14.0/fastq_screen --conf ${params.fqsConfig} --threads ${task.cpus} --tag --filter --0---------00---- --subset 0 ${pair_id}.R1.fastq.gz
    \$APP/fastq_screen/v0.14.0/fastq_screen --conf ${params.fqsConfig} --threads ${task.cpus} --tag --filter --0---------00---- --subset 0 ${pair_id}.R2.fastq.gz

    """
}




process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    path mqcFiles

    output:
    file('multiqc_report.html') optional true

    script:
    """
    module load multiqc
    multiqc -v .
    """
}

workflow {
    fastqc(read_pairs)
    fqscreenQc(read_pairs)
    fqscreenFilt(read_pairs)
    multiqc(fastqc.out.mix(fqscreenQc.out).collect())
}




