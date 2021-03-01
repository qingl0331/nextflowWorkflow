#!/usr/bin/env nextflow
 
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
    .set { read_pairs } 
 
 
/*
 * Step 1. Trimming with cutadapt
 */
process trim {
    cpus 20
    tag "$pair_id"
    
    input:
    set pair_id, file(reads) from read_pairs
    output:
    set pair_id, file("${pair_id}.trimmed_*.fastq.gz") into tobwaMem
    set pair_id, file("${pair_id}.trimmed_*.fastq.gz") into tofqc
    file("${pair_id}.cutadapt.txt") into forMqc1 
 
    """
    module load cutadapt
    cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o ${pair_id}.trimmed_1.fastq.gz -p ${pair_id}.trimmed_2.fastq.gz -j ${task.cpus} $reads > ${pair_id}.cutadapt.txt

    """
}

/*
 * Step 2. QC with fastqc
 */
process fastqc {
       
    input:
     set pair_id, file("${pair_id}.trimmed_*.fastq.gz") from tofqc
    output:
     file("fastqc_${pair_id}_logs") into forMqc2
 
    """
    mkdir fastqc_${pair_id}_logs
    module load fastqc
    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${pair_id}.trimmed_*.fastq.gz
    """
}
/*
 * Step 3. bwa mem alignment & samtools sort & picard markDuplicates
 */
process alignSortmkdup {
    publishDir "${params.outdir}/finalBam", mode: 'copy', pattern: '*mkdup.ba*'
    input:
     set pair_id, file("${pair_id}.trimmed_*.fastq.gz") from tobwaMem
    output:
     tuple val(pair_id), file("${pair_id}.mkdup.bam"), file("${pair_id}.mkdup.bai") into toalignSum 
     tuple val(pair_id), file("${pair_id}.mkdup.bam"), file("${pair_id}.mkdup.bai") into toalignCI 
     tuple val(pair_id), file("${pair_id}.mkdup.bam"), file("${pair_id}.mkdup.bai") into toalignCW 
     tuple val(pair_id), file("${pair_id}.mkdup.bam"), file("${pair_id}.mkdup.bai") into togvcf 
     file("${pair_id}.markduplicates.txt") into forMqc3 
 
    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app
    \$APP/bwa/0.7.8/bwa mem -v 2 -t ${task.cpus} -R "@RG\tID:${pair_id}\tPL:illuminaNovaSeq\tLB:NexteraWGS\tSM:${pair_id}\tCN:HCI\tPU:191115" ${params.index} ${pair_id}.trimmed_1.fastq.gz ${pair_id}.trimmed_2.fastq.gz | \$APP/samtools/1.8/samtools view -b -o ${pair_id}.raw.bam -

    \$APP/samtools/1.8/samtools sort ${pair_id}.raw.bam -o ${pair_id}.sorted.bam
    \$APP/samtools/1.8/samtools index ${pair_id}.sorted.bam	
    java -jar \$APP/picard/2.18.26/picard.jar MarkDuplicates I=${pair_id}.sorted.bam O=${pair_id}.mkdup.bam M=${pair_id}.markduplicates.txt REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000 ASSUME_SORTED=true CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING USE_JDK_DEFLATER=true USE_JDK_INFLATER=true 

    """
}

/*
 * Step 4. AlignmentSummary with picard
 */
process alignSum {
    input:
    tuple val(pair_id), file("${pair_id}.mkdup.bam"), file("${pair_id}.mkdup.bai") from toalignSum
    output:
     file("${pair_id}.alignmentSummaryMetrics.txt") into forMqc4
    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app

    java -jar \$APP/picard/2.18.26/picard.jar CollectAlignmentSummaryMetrics I=${pair_id}.mkdup.bam R=${params.ref} O=${pair_id}.alignmentSummaryMetrics.txt VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING 

    """
}


/*
 * Step 5. CollectInsertSize with picard
 */
process alignCI {
    input:
     tuple val(pair_id), file("${pair_id}.mkdup.bam"), file("${pair_id}.mkdup.bai") from toalignCI
    output:
     file("${pair_id}.insertSizeMetrics.txt") into forMqc5
    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app

    java -jar \$APP/picard/2.18.26/picard.jar CollectInsertSizeMetrics I=${pair_id}.mkdup.bam O=${pair_id}.insertSizeMetrics.txt H=${pair_id}.insert_size_histogram.pdf M=0.5 VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING 

    """
}


/*
 * Step 6. CollectWGSMetrics with picard
 */
process alignCW {
    input:
     tuple val(pair_id), file("${pair_id}.mkdup.bam"), file("${pair_id}.mkdup.bai") from toalignCW
    output:
     file("${pair_id}.wgsMetrics.txt") into forMqc6
    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app

    java -jar \$APP/picard/2.18.26/picard.jar CollectWgsMetrics I=${pair_id}.mkdup.bam R=${params.ref} O=${pair_id}.wgsMetrics.txt MINIMUM_MAPPING_QUALITY=20 MINIMUM_BASE_QUALITY=20 VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING 

    """
}

/*
 * Step 7. Generate gVCF file
 */
process gVcf {
    input:
    tuple val(pair_id), file("${pair_id}.mkdup.bam"), file("${pair_id}.mkdup.bai") from togvcf
    output:
     file("${pair_id}.g.vcf*") into forMerge
    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app

    \$APP/gatk/4.1.4.1/gatk HaplotypeCaller -ERC GVCF -R ${params.ref} -I ${pair_id}.mkdup.bam -O ${pair_id}.g.vcf --ploidy 1 --max-reads-per-alignment-start 0 

    """
}
/*
 * Step 8. Combine, genotype gvcfs, filter and annotate
 */
process cgfa {
    publishDir params.outdir, mode: 'copy' 
    input:
     file(gVcfs) from forMerge.collect()
    output:
     file("combined.g.vcf")  
     file("ggvcf.vcf")  
     file("ggvcf.filtered.vcf")  
     file("filtered_final.ann.vcf")  
     file("snpEff.summary")  
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

/*
 * Step 9. QC summary with multiqc
 */


process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    file('./*') from forMqc6.mix(forMqc1,forMqc2,forMqc3,forMqc4,forMqc5).collect()

    output:
    file('multiqc_report.html') optional true

    script:
    """
    module load multiqc
    multiqc -v .
    """
}


workflow.onComplete { 
	println ( workflow.success ? "Done!Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
