#!/usr/bin/env nextflow
 
/* 
 * Proof of concept Nextflow based microbiome pipeline
 * 
 */ 

 
/*
 * Defines some parameters in order to specify 
 * read pairs by using the command line options
 */
params.reads = "/scratch/general/lustre/u0762203/nextflowTest/testReads/*_R{1,2}*.fastq.gz"
params.g38Dir = "/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/Human/GRCh38/bowtie2"
params.outdir = 'results'
log.info """
         MICROBIOME   P I P E L I N E:metagenome    
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
 * Step 1. remove optical duplicates from NovaSeq 
 */ 
process dedup {
    tag "$pair_id"
     
    input:
    set pair_id, file(reads) from read_pairs
 
    output:
    set pair_id, "${pair_id}.clump*.fq.gz" into totrim
    set pair_id, file("${pair_id}.clump1.fq.gz") into tofqc 
    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app
    \$APP/BBmap/v38.34/clumpify.sh in1=${reads[0]} in2=${reads[1]} out1=${pair_id}.clump1.fq.gz out2=${pair_id}.clump2.fq.gz dupedist=12000 dedupe=t optical=t
    """
}
/*
 * Step 2. Trimming & decontam with kneaddata
 */
process trim {
    module 'bioBakery/1.7:bowtie2:trimmomatic'
    cpus 18
    //cpus 4
    tag "$pair_id"
    
    input:
    set pair_id, file("${pair_id}.clump*.fq.gz") from totrim
 
    output:
    set pair_id, file("${pair_id}.fq") into tohumann2
    file("./kneaddata_out/${pair_id}.clump1_kneaddata.log") into forMqc
 
    """
    runbioBakery kneaddata --input ${pair_id}.clump1.fq.gz --input ${pair_id}.clump2.fq.gz -o kneaddata_out/ -db ${params.g38Dir}/GRCh38 -t ${task.cpus} --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output
    runbioBakery kneaddata_read_count_table --input kneaddata_out --output ${pair_id}.kneaddata_read_counts.human.txt 
    cat kneaddata_out/${pair_id}.clump1_kneaddata_paired_{1,2}.fastq > ${pair_id}.fq
    """
}

/*
 * Step 3. QC with fastqc
 */
process fastqc {
       
    input:
     set pair_id, file("${pair_id}.clump1.fq.gz") from tofqc
     
    output:
     file("fastqc_${pair_id}_logs") into forMqc2
 
    """
    mkdir fastqc_${pair_id}_logs
    module load fastqc
    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${pair_id}.clump1.fq.gz
    """
}
/*
 * Step 4. HumanN2 classification
 */
process humann2 {
    module 'bioBakery/1.7'
    publishDir "$baseDir/humann2_results", mode: 'copy'  
    input:
     set pair_id, file("${pair_id}.fq") from tohumann2
    output:
     file("./humann2_out/${pair_id}_humann2_temp/${pair_id}.log") into forMqc3
     file("./humann2_out/${pair_id}_humann2_temp/${pair_id}_metaphlan_bugs_list.tsv") into tomerge 
     file ("./humann2_out/${pair_id}_genefamilies.tsv")
     file ("./humann2_out/${pair_id}_pathabundance.tsv")
     file ("./humann2_out/${pair_id}_pathcoverage.tsv")
 
    """
    #echo ${pair_id}

    runbioBakery humann2 --threads ${task.cpus} --input ${pair_id}.fq --output humann2_out/ --protein-database /uufs/chpc.utah.edu/sys/installdir/bioBakery/databases/uniref/ --nucleotide-database=/uufs/chpc.utah.edu/sys/installdir/bioBakery/databases/chocophlan --prescreen-threshold 0.0001  
    #cp humann2_out/*temp/*metaphlan_bugs_list.tsv  humann2_out/
    #cp  humann2_out/*temp/*bowtie2.txt humann2_out/
    """
}
/*
 * Step 5. HumanN2 merge and functional analysis
 */
process merge {
    module 'bioBakery/1.7'
    publishDir params.outdir, mode: 'copy' 
    input:
     file(bugsListTsv) from tomerge.collect()
    output:
     file("humann2*tsv")  
     file("metaphlan2_merged.txt")  
    """
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_results/ --file_name pathabundance --output humann2_pathabundance.tsv  
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_results/ --file_name pathcoverage --output humann2_pathcoverage.tsv
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_results/ --file_name genefamilies --output humann2_genefamilies.tsv
    runbioBakery humann2_renorm_table --input humann2_pathabundance.tsv --units relab --output humann2_pathabundance_relab.tsv
    runbioBakery humann2_renorm_table --input humann2_genefamilies.tsv --units relab --output humann2_genefamilies_relab.tsv
    runbioBakery merge_metaphlan_tables.py ${bugsListTsv} > metaphlan2_merged.txt
    sed -i 's/_metaphlan_bugs_list//g' metaphlan2_merged.txt
    """
}

/*
 * Step 6. QC summary with multiqc
 */


process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    file('./*') from forMqc3.mix(forMqc,forMqc2).collect()

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
