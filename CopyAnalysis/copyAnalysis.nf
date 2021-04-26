#!/usr/bin/env nextflow

nextflow.enable.dsl=2 


/* 
 * Proof of concept Nextflow based Sentieon germline CNV calling pipeline
 * 
 */ 

 
/*
 * Defines some parameters in order to specify 
 * input bam, its index, reference fasta, panel of normal, exome bed file and ucsc annotation table  by using the command line options
 */

params.bam="bam"
params.bamIdx="bam.bai"
params.fasta = "/uufs/chpc.utah.edu/common/PE/hci-bioinformatics1/TNRunner/Indexes/B38IndexForBwa-0.7.17/hs38DH.fa"
params.avatarPon = '/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/qli/cnvBenchmark/avatarPon'
params.exomBed = '/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/qli/cnvBenchmark/hg38NimIdtMergedPad150bp.bed'
params.ucscAno = '/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/qli/db/hg38RefSeq8Aug2018_Merged.ucsc.gz'
params.outdir = 'results'

bam=file(params.bam)
bamName=bam.getName()
id=(bamName=~ /(\S+)_NormalDNA_Hg38_final.bam/)[0][1]

log.info """
         Sentieon germline CNV calling pipeline    
         =============================
         bam : ${params.bam}
         outdir: ${params.outdir}
         """
         .stripIndent()

 

process callCnvs {
    publishDir "${params.outdir}/sentieonOut", mode: 'copy'
    input:
    file(params.bamIdx) 
    output:
    path("${id}"), emit: cnvCall 
    path("${id}.t*") 
 
    script:

    """
    sentieon driver -t ${task.cpus} -r ${params.fasta} -i ${params.bam} --algo CNV --pon ${params.avatarPon} ${id}

    """
}

process makeBed {
    input:
    path("${id}") 
    output:
    path("${id}.intersect.bed") 
 
    script:

    """
    /uufs/chpc.utah.edu/common/HIPAA/u0762203/bin/stCnvCall2bed.pl ${id}
    bedtools intersect -a ${id}.bed -b ${params.exomBed} > ${id}.intersect.bed

    """
}

process annot {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path("${id}.intersect.bed") 
    output:
    path("${id}.cnvCallWithGenes.txt.gz")

    script:
 
    """
    java -jar /uufs/chpc.utah.edu/common/HIPAA/u0762203/USeq_9.2.8/Apps/AnnotateBedWithGenes -p 100 -g -b ${id}.intersect.bed -r ${id}.cnvCallWithGenes.txt.gz -u ${params.ucscAno} 

    """
}



workflow {
    callCnvs(params.bamIdx)
    makeBed(callCnvs.out.cnvCall)
    annot(makeBed.out)
}




