#!/usr/bin/env nextflow
 
/* 
 * Proof of concept Nextflow in Nextflow based microbiome pipeline
 * 
 */ 

 
/*
 * Defines some parameters in order to specify 
 * read pairs by using the command line options
 */
params.reads = "/scratch/general/lustre/u0762203/microbiomePipeSingu/dnaReads/*_R{1,2}*.fastq.gz"
params.outdir = '/scratch/general/lustre/u0762203/microbiomePipeSingu/dnaResults'
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
 * Step 1. genome run individual level 
 */ 
process gri {
    tag "$pair_id"
     
    input:
    set pair_id, file(reads) from read_pairs
 
    output:
    set pair_id, file("./humann2_out.${pair_id}/${pair_id}_humann2_temp/${pair_id}_metaphlan_bugs_list.tsv") into tomerge 
    file("./kneaddata_out.${pair_id}/${pair_id}.clump1_kneaddata.log") into forMqc1
    file("./fastqc_${pair_id}_logs") into forMqc2
    file("./humann2_out.${pair_id}/${pair_id}_humann2_temp/${pair_id}.log") into forMqc3
    """
    bash -x ${baseDir}/mb.run.1.sh
    """
}
/*
 * Step 2. QC summary with multiqc
 */


process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    file('./*') from forMqc3.mix(forMqc1,forMqc2).collect()

    output:
    file('multiqc_report.html') optional true

    script:
    """
    module load multiqc
    multiqc -v .
    """
}

/*
 * Step 3. HumanN2 merge and functional analysis
 */
process merge {
    module 'bioBakery/1.7'
    publishDir params.outdir, mode: 'copy', pattern: 'humann2*tsv' || 'metaphlan2_merged.txt' 
    input:
     file(bugsListTsv) from tomerge.collect()
    output:
     file("humann2*tsv")  
     file("metaphlan2_merged.txt")  
    """
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_out/ --file_name pathabundance --output humann2_pathabundance.tsv  
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_out/ --file_name pathcoverage --output humann2_pathcoverage.tsv
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_out/ --file_name genefamilies --output humann2_genefamilies.tsv
    runbioBakery humann2_renorm_table --input humann2_pathabundance.tsv --units relab --output humann2_pathabundance_relab.tsv
    runbioBakery humann2_renorm_table --input humann2_genefamilies.tsv --units relab --output humann2_genefamilies_relab.tsv
    runbioBakery merge_metaphlan_tables.py ${bugsListTsv} > metaphlan2_merged.txt
    sed -i 's/_metaphlan_bugs_list//g' metaphlan2_merged.txt
    """
}


