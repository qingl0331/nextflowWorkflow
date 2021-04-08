#!/usr/bin/env nextflow

nextflow.enable.dsl=2 


/* 
 * Proof of concept Nextflow based GATK4 joint variant calling pipeline
 * 
 */ 

 
/*
 * Defines some parameters in order to specify 
 * read pairs by using the command line options
 */
params.reads = "/scratch/general/lustre/u0762203/nextflowTest/testReads/run1/*_R{1,2}*.fastq.gz"
params.index = '/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/Bacteria/Ecoli/ecoli119/E.coli.119.GCF_005222045.1_ASM522204v1_genomic.bwa'
params.ref = '/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/Bacteria/Ecoli/ecoli119/E.coli.119.GCF_005222045.1_ASM522204v1_genomic.fna'
params.snpeffDb = 'ecoli119'
params.outdir = 'results'
log.info """
         WGS Variant Calling  P I P E L I N E:Bacteria    
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


process trim {
    
    input:
    tuple val(pair_id), path(read) 
    output:
    tuple val(pair_id), path("${pair_id}.trimmed_*.fastq.gz"), emit: fq 
    path("${pair_id}.cutadapt.txt") , emit: report
 
    script:

    """
    module load cutadapt
    cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o ${pair_id}.trimmed_1.fastq.gz -p ${pair_id}.trimmed_2.fastq.gz -j ${task.cpus} ${read} > ${pair_id}.cutadapt.txt

    """
}

process fastqc {
    input:
    tuple val(pair_id), path("${pair_id}.trimmed_*.fastq.gz") 
    output:
    path("fastqc_${pair_id}_logs") 
 
    script:

    """
    mkdir fastqc_${pair_id}_logs
    module load fastqc
    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${pair_id}.trimmed_*.fastq.gz
    """
}

process alignSortmkdup {
    publishDir "${params.outdir}/finalBam", mode: 'copy', pattern: '*mkdup.ba*'
    input:
    tuple val(pair_id), path("${pair_id}.trimmed_*.fastq.gz") 
    output:
     tuple val(pair_id), path("${pair_id}.mkdup.bam"), path("${pair_id}.mkdup.bai"), emit: bam
     path("${pair_id}.markduplicates.txt"), emit: log  

    script:
 
    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app
    \$APP/bwa/0.7.8/bwa mem -v 2 -t ${task.cpus} -R "@RG\tID:${pair_id}\tPL:illuminaNovaSeq\tLB:NexteraWGS\tSM:${pair_id}\tCN:HCI\tPU:191115" ${params.index} ${pair_id}.trimmed_1.fastq.gz ${pair_id}.trimmed_2.fastq.gz | \$APP/samtools/1.8/samtools view -b -o ${pair_id}.raw.bam -

    \$APP/samtools/1.8/samtools sort ${pair_id}.raw.bam -o ${pair_id}.sorted.bam
    \$APP/samtools/1.8/samtools index ${pair_id}.sorted.bam	
    java -jar \$APP/picard/2.18.26/picard.jar MarkDuplicates I=${pair_id}.sorted.bam O=${pair_id}.mkdup.bam M=${pair_id}.markduplicates.txt REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000 ASSUME_SORTED=true CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING USE_JDK_DEFLATER=true USE_JDK_INFLATER=true 

    """
}


process alignSum {
    input:
    tuple val(pair_id), path("${pair_id}.mkdup.bam"), path("${pair_id}.mkdup.bai")
    output:
    path("${pair_id}.alignmentSummaryMetrics.txt")

    script:

    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app

    java -jar \$APP/picard/2.18.26/picard.jar CollectAlignmentSummaryMetrics I=${pair_id}.mkdup.bam R=${params.ref} O=${pair_id}.alignmentSummaryMetrics.txt VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING 

    """
}


process alignCI {
    input:
     tuple val(pair_id), path("${pair_id}.mkdup.bam"), path("${pair_id}.mkdup.bai") 
    output:
    path("${pair_id}.insertSizeMetrics.txt")

    script:

    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app

    java -jar \$APP/picard/2.18.26/picard.jar CollectInsertSizeMetrics I=${pair_id}.mkdup.bam O=${pair_id}.insertSizeMetrics.txt H=${pair_id}.insert_size_histogram.pdf M=0.5 VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING 

    """
}


process alignCW {
    input:
     tuple val(pair_id), path("${pair_id}.mkdup.bam"), path("${pair_id}.mkdup.bai")
    output:
    path("${pair_id}.wgsMetrics.txt")
 
    script:

    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app

    java -jar \$APP/picard/2.18.26/picard.jar CollectWgsMetrics I=${pair_id}.mkdup.bam R=${params.ref} O=${pair_id}.wgsMetrics.txt MINIMUM_MAPPING_QUALITY=20 MINIMUM_BASE_QUALITY=20 VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING 

    """
}

process gVcf {
    input:
    tuple val(pair_id), path("${pair_id}.mkdup.bam"), path("${pair_id}.mkdup.bai")
    output:
    path("${pair_id}.g.vcf*")

    script:

    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app

    \$APP/gatk/4.1.4.1/gatk HaplotypeCaller -ERC GVCF -R ${params.ref} -I ${pair_id}.mkdup.bam -O ${pair_id}.g.vcf --ploidy 1 --max-reads-per-alignment-start 0 

    """
}

process cgfa {
    publishDir params.outdir, mode: 'copy' 
    input:
     path(gVcfs)
    output:
     path("combined.g.vcf")  
     path("ggvcf.vcf")  
     path("ggvcf.filtered.vcf")  
     path("filtered_final.ann.vcf")  
     path("snpEff.summary")  

    script:

    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app
    
    mkdir gvcfs
    cp ${gVcfs} gvcfs/
    bash -x /uufs/chpc.utah.edu/common/home/u0762203/bin/print.cg.sh ${params.ref} 
    bash -x cg.sh

    \$APP/gatk/4.1.4.1/gatk GenotypeGVCFs -R ${params.ref} -V ./combined.g.vcf -O ggvcf.vcf
    /uufs/chpc.utah.edu/common/home/u0762203/bin/filter.vcf.pl ggvcf.vcf gvcfs/ 
    java -jar \$APP/snpEff/snpEff.jar -v ${params.snpeffDb} -s snpEff.summary ggvcf.filtered.vcf > filtered_final.ann.vcf
 
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
    trim(read_pairs)
    fastqc(trim.out.fq)
    alignSortmkdup(trim.out.fq)
    alignSum(alignSortmkdup.out.bam)
    alignCI(alignSortmkdup.out.bam)                                            
    alignCW(alignSortmkdup.out.bam)                                            
    gVcf(alignSortmkdup.out.bam)  
    cgfa(gVcf.out.collect())
    multiqc(trim.out.report.mix(fastqc.out).mix(alignSortmkdup.out.log).mix(alignSum.out).mix(alignCI.out).mix(alignCW.out))
}




