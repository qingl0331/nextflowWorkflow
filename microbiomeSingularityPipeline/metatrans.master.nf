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
params.rnaReads = "/scratch/general/lustre/u0762203/microbiomePipeSingu/rnaReads/*_R{1,2}*.fastq.gz"
params.sampSheet="/scratch/general/lustre/u0762203/microbiomePipeSingu/run/SampleSequencingDone_Master_190409missing15411X45.txt"
params.outdir = '/scratch/general/lustre/u0762203/microbiomePipeSingu/dnaResults'
params.rnaOutdir = '/scratch/general/lustre/u0762203/microbiomePipeSingu/rnaResults'
log.info """
         MICROBIOME   P I P E L I N E:metagenome and matched transcriptome    
         =============================
         reads : ${params.reads}
         rnaReads : ${params.rnaReads}
         outdir: ${params.outdir}
         rnaOutdir: ${params.rnaOutdir}
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
 * Step 3. HumanN2 merge and functional analysis, prepare for matched transcriptome normalization
 */
process merge {
    module 'bioBakery/1.7'
    publishDir "$baseDir/genome_normalization", mode: 'copy' , pattern: '*norm*tsv' 
    publishDir params.outdir, mode: 'copy', pattern: 'humann2*tsv' || 'metaphlan2_merged.txt' 
    input:
     file(bugsListTsv) from tomerge.collect()
    output:
     file("humann2*tsv")  
     file("*norm*tsv") into normTsv 
     file("metaphlan2_merged.txt")  
    """
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_out/ --file_name pathabundance --output humann2_pathabundance.tsv  
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_out/ --file_name pathcoverage --output humann2_pathcoverage.tsv
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_out/ --file_name genefamilies --output humann2_genefamilies.tsv
    runbioBakery humann2_renorm_table --input humann2_pathabundance.tsv --units relab --output humann2_pathabundance_relab.tsv
    runbioBakery humann2_renorm_table --input humann2_genefamilies.tsv --units relab --output humann2_genefamilies_relab.tsv
    runbioBakery merge_metaphlan_tables.py ${bugsListTsv} > metaphlan2_merged.txt
    sed -i 's/_metaphlan_bugs_list//g' metaphlan2_merged.txt
    changePartFileName.master.nf.g2t4metaphlan.pl ${baseDir}/metaphlan2_out ${params.sampSheet} 
    changePartFileName.master.nf.g2t4abundGenefam.pl ${baseDir}/humann2_out ${params.sampSheet} 
    """
}

/*
 *For transcriptome part:  Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.rnaReads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairsRna }


/*
 * Step 4. transcriptome run individual level:HumanN2 classification for paired transcriptome 
 */ 
process tri {
    tag "$pair_idRna"
     publishDir "$baseDir/humann2_resultsRna", mode: 'copy'
    input:
    set pair_idRna, file(reads) from read_pairsRna
   file(nt) from normTsv.collect() //just to make sure transcriptome run start after genome run is finished
 
    output:
     file("./${pair_idRna}.cutadapt.txt") into forMqc4
     file("./fastqc_${pair_idRna}_logs") into forMqc5
     file("./humann2_out.${pair_idRna}/${pair_idRna}.no_rRNA_humann2_temp/${pair_idRna}.no_rRNA.log") into forMqc6 //forMqc6 rf -> forMqc7, just to make sure transcriptome merge is after transcriptome individual run
     file ("./${pair_idRna}.rRNA.log") into forMqc8 
     file ("./humann2_out.${pair_idRna}/${pair_idRna}.no_rRNA_genefamilies.tsv") into rGf
     file ("./humann2_out.${pair_idRna}/${pair_idRna}.no_rRNA_pathabundance.tsv")
     file ("./humann2_out.${pair_idRna}/${pair_idRna}.no_rRNA_pathcoverage.tsv")

    """
    bash -x ${baseDir}/mb.run.2.sh
    """
}
/*
 * Step 5. HumanN2 merge matched transcriptome normalization
 */
process mergeNormrna {
    module 'bioBakery/1.7'
    publishDir params.rnaOutdir, mode: 'copy'
    input:
     file(gf) from rGf.collect()
    output:
     file("humann2*DNAnormalizedRNA*tsv")
     """
#Join HUMAnN2 output per sample into one table.
    mkdir -p humann2_4normalize_final_out
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_resultsRna/ --file_name pathabundance --output humann2_4normalize_final_out/humann2_pathabundance.tsv
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_resultsRna/ --file_name genefamilies --output humann2_4normalize_final_out/humann2_genefamilies.tsv
#join matched genome humann2 output per sample into one table
    mkdir -p humann2_final_matched.genome.out
    runbioBakery humann2_join_tables -s --input ${baseDir}/genome_normalization/ --file_name pathabundance --output humann2_final_matched.genome.out/humann2_pathabundance.matchedGenome.tsv
    runbioBakery humann2_join_tables -s --input ${baseDir}/genome_normalization/ --file_name genefamilies --output humann2_final_matched.genome.out/humann2_genefamilies.matchedGenome.tsv
#DEBUG!!copy humann2_genefamilies.tsv, humann2_pathabundance.tsv,*matchedGenome.tsv  to basedir
    #mkdir -p ${baseDir}/humann2_4normalize_final_out/
    #cp humann2_4normalize_final_out/*.tsv ${baseDir}/humann2_4normalize_final_out/
    #mkdir -p ${baseDir}/humann2_final_matched.genome.out
    #cp humann2_final_matched.genome.out/*tsv ${baseDir}/humann2_final_matched.genome.out/
    #mkdir -p ${baseDir}/humann2_genomeNormalized
    #runbioBakery humann2_rna_dna_norm -d ${baseDir}/humann2_final_matched.genome.out/humann2_pathabundance.matchedGenome.tsv -r ${baseDir}/humann2_4normalize_final_out/humann2_pathabundance.tsv -o ${baseDir}/humann2_genomeNormalized/humann2_pathabundance.DNAnormalizedRNA
   #runbioBakery humann2_rna_dna_norm -d ${baseDir}/humann2_final_matched.genome.out/humann2_genefamilies.matchedGenome.tsv -r ${baseDir}/humann2_4normalize_final_out/humann2_genefamilies.tsv -o ${baseDir}/humann2_genomeNormalized/humann2_genefamilies.DNAnormalizedRNA
#Re-normalize RNASeq gene family and pathway abundances by corresponding DNA-level outputs to quantify microbial expression independent of gene copy number 
    runbioBakery humann2_rna_dna_norm -d humann2_final_matched.genome.out/humann2_pathabundance.matchedGenome.tsv -r humann2_4normalize_final_out/humann2_pathabundance.tsv -o humann2_pathabundance.DNAnormalizedRNA
    runbioBakery humann2_rna_dna_norm -d humann2_final_matched.genome.out/humann2_genefamilies.matchedGenome.tsv -r humann2_4normalize_final_out/humann2_genefamilies.tsv -o humann2_genefamilies.DNAnormalizedRNA
   """
}
/*
 * Step 6. HumanN2 merge and functional analysis transcriptome
 */
process mergerna {
    module 'bioBakery/1.7'
    publishDir params.rnaOutdir, mode: 'copy'
    input:
     file(hn2rnaLog) from forMqc6.collect()
    output:
     file("humann2*tsv")
     file("*rf.log") into forMqc7
     """
#Join HUMAnN2 output per sample into one table.
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_resultsRna/ --file_name pathabundance --output humann2_pathabundance.tsv
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_resultsRna/ --file_name pathcoverage --output humann2_pathcoverage.tsv
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_resultsRna/ --file_name genefamilies --output humann2_genefamilies.tsv
#Re-normalize gene family and pathway abundances
    runbioBakery humann2_renorm_table --input humann2_pathabundance.tsv --units relab --output humann2_pathabundance_relab.tsv
    runbioBakery humann2_renorm_table --input humann2_genefamilies.tsv --units relab --output humann2_genefamilies_relab.tsv
    ##reformat the humann2 log for multiqc
    for file in ${hn2rnaLog};do
        changePartFileName.rfLog.pl \$file
        rm \$file
    done


     """
}


/*
/*
 * Step 7. QC summary with multiqc for transcriptome
 */


process multiqcrna {
    publishDir params.rnaOutdir, mode:'copy'

    input:
    file('./*') from forMqc7.mix(forMqc4,forMqc5,forMqc8).collect()

    output:
    file('multiqc_report.html')

    script:
    """
    module load multiqc
    multiqc -v .

    """
}


workflow.onComplete {
        println ( workflow.success ? "Done!Open the following report in your browser --> $params.outdir/multiqc_report.html\n$params.rnaOutdir/multiqc_report.rna.html\n" : "Oops .. something went wrong" )
}


