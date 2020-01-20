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
params.rnaReads = "/scratch/general/lustre/u0762203/nextflowTest/rnaTestReads/*_R{1,2}*.fastq.gz"
params.g38Dir = "/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/Human/GRCh38/bowtie2"
params.sampSheet="/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/SampleSequencingDone_Master_190409missing15411X45.txt"
params.outdir = 'dnaResults'
params.rnaOutdir = 'rnaResults'
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
 * Step 5. HumanN2 merge and functional analysis, prepare for matched transcriptome normalization
 */
process merge {
    module 'bioBakery/1.7'
    publishDir "$baseDir/genome_normalization", mode: 'move' , pattern: '*norm*tsv' 
    publishDir params.outdir, mode: 'copy', pattern: 'humann2*tsv' || 'metaphlan2_merged.txt' 
    input:
     file(bugsListTsv) from tomerge.collect()
    output:
     file("humann2*tsv")  
     file("*norm*tsv") into normTsv 
     file("metaphlan2_merged.txt")  
    """
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_results/ --file_name pathabundance --output humann2_pathabundance.tsv  
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_results/ --file_name pathcoverage --output humann2_pathcoverage.tsv
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_results/ --file_name genefamilies --output humann2_genefamilies.tsv
    runbioBakery humann2_renorm_table --input humann2_pathabundance.tsv --units relab --output humann2_pathabundance_relab.tsv
    runbioBakery humann2_renorm_table --input humann2_genefamilies.tsv --units relab --output humann2_genefamilies_relab.tsv
    runbioBakery merge_metaphlan_tables.py ${bugsListTsv} > metaphlan2_merged.txt
    sed -i 's/_metaphlan_bugs_list//g' metaphlan2_merged.txt
    changePartFileName.nf.g2t4metaphlan.pl ${baseDir}/humann2_results ${params.sampSheet} 
    changePartFileName.nf.g2t4abundGenefam.pl ${baseDir}/humann2_results ${params.sampSheet} 
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

/*
 *For transcriptome part:  Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */
Channel
    .fromFilePairs( params.rnaReads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairsRna } 


/*
 * Step 7. remove optical duplicates from NovaSeq 
 */ 
process dedupRna {
    tag "$pair_idRna"
     
    input:
    set pair_idRna, file(reads) from read_pairsRna
 
    output:
    set pair_idRna, "${pair_idRna}.clump*.fq.gz" into totrimRna
    """
    export APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app
    \$APP/BBmap/v38.34/clumpify.sh in1=${reads[0]} in2=${reads[1]} out1=${pair_idRna}.clump1.fq.gz out2=${pair_idRna}.clump2.fq.gz dupedist=12000 dedupe=t optical=t
    """
}


/*
 * Step 8. Trimming with cutadapt
 */
process trimRna {
    tag "$pair_idRna"
    input:
    set pair_idRna, file("${pair_idRna}.clump*.fq.gz") from totrimRna
     
    output:
    set pair_idRna, file("${pair_idRna}.trimmed_R{1,2}.fq") into tostitch
    set pair_idRna, file("${pair_idRna}.trimmed_R1.fq") into tofqcRna
    file("*cutadapt.txt") into forMqc4
    """
     module load cutadapt
     cutadapt -j 0 -O 6 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${pair_idRna}.trimmed_R1.fq -p ${pair_idRna}.trimmed_R2.fq -j ${task.cpus} -O 6 -m 20 ${pair_idRna}.clump*.fq.gz > ${pair_idRna}.cutadapt.txt 
    """

}
/*
 * Step 9. stitch paired end with pear
 */
process stitch {
   module 'pear/0.9.6'
   publishDir "$baseDir/pearLog", mode: 'copy' , pattern: '*stitch.log.txt' 
   input:
   set pair_idRna, file(trimmedReads) from tostitch
  
   output:
   set pair_idRna, file("${pair_idRna}.assembled.fastq") into tosortmerna 
   """
   pear -f ${trimmedReads[0]} -r ${trimmedReads[1]} -j ${task.cpus} -o ${pair_idRna} >> ${pair_idRna}.stitch.log.txt
   """
}

/*
 * Step 10. QC with fastqc
 */
process fastqcRna {

    input:
     set pair_idRna, file("${pair_idRna}.trimmed_R1.fq") from tofqcRna

    output:
     file("fastqc_${pair_idRna}_logs") into forMqc5

    """
    mkdir fastqc_${pair_idRna}_logs
    module load fastqc
    fastqc -o fastqc_${pair_idRna}_logs -f fastq -q ${pair_idRna}.trimmed_R1.fq
    """
}

/*
 * Step 11. sortmerna to remove rRNA
 */
process sortmerna{
   input:
   set pair_idRna, file("${pair_idRna}.assembled.fastq") from tosortmerna
   output:
   set pair_idRna, file("${pair_idRna}.no_rRNA.fastq") into tohumann2rna
   file ("*log") into forMqc8
   """
   module load sortmerna
   export DB=/uufs/chpc.utah.edu/sys/installdir/sortmerna/sortmerna-2.1b/rRNA_databases
   sortmerna --ref \$DB/silva-bac-16s-id90.fasta,\$DB/silva-bac-16s-id90.idx:\
\$DB/silva-bac-23s-id98.fasta,\$DB/silva-bac-23s-id98.idx:\
\$DB/silva-euk-18s-id95.fasta,\$DB/silva-euk-18s-id95.idx:\
\$DB/silva-euk-28s-id98.fasta,\$DB/silva-euk-28s-id98.idx \
--reads ${pair_idRna}.assembled.fastq --aligned ${pair_idRna}.rRNA --other ${pair_idRna}.no_rRNA --fastx --log -a ${task.cpus} -v
   grep -v '^\$' ${pair_idRna}.no_rRNA.fastq > tmpfile
   mv tmpfile ${pair_idRna}.no_rRNA.fastq
   """
}
/*
 * Step 12. HumanN2 classification for paired transcriptome
 */
process humann2rna {
    module 'bioBakery/1.7'
    publishDir "$baseDir/humann2_resultsRna", mode: 'copy'  
    input:
     set pair_idRna, file("${pair_idRna}.no_rRNA.fastq") from tohumann2rna
     file(nt) from normTsv.collect() //just to make sure transcriptome run start after genome run is finished
    output:
     file("./humann2_out/${pair_idRna}.no_rRNA_humann2_temp/${pair_idRna}.no_rRNA.log") into forMqc6
     file ("./humann2_out/${pair_idRna}.no_rRNA_genefamilies.tsv") into rGf  
     file ("./humann2_out/${pair_idRna}.no_rRNA_pathabundance.tsv")
     file ("./humann2_out/${pair_idRna}.no_rRNA_pathcoverage.tsv")
 
    """
    prName=`echo ${pair_idRna}| cut -d "_" -f 1`
    runbioBakery humann2 --threads ${task.cpus} --input ${pair_idRna}.no_rRNA.fastq --output humann2_out/ --protein-database /uufs/chpc.utah.edu/sys/installdir/bioBakery/databases/uniref/ --nucleotide-database=/uufs/chpc.utah.edu/sys/installdir/bioBakery/databases/chocophlan --taxonomic-profile ${baseDir}/genome_normalization/\${prName}.norm.tsv --prescreen-threshold 0.0001  
    """
}

/*
 * Step 13. HumanN2 merge matched transcriptome normalization
 */
process mergeNormrna {
    module 'bioBakery/1.7'
    publishDir params.rnaOutdir, mode: 'copy' 
    input:
     file(gf) from rGf.collect()
    output:
     file("humann2*DNAnormalizedRNA*")  
    """
#Join HUMAnN2 output per sample into one table.
    mkdir -p humann2_4normalize_final_out
    runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_resultsRna/humann2_out/ --file_name pathabundance --output humann2_4normalize_final_out/humann2_pathabundance.tsv
runbioBakery humann2_join_tables -s --input ${baseDir}/humann2_resultsRna/humann2_out/ --file_name genefamilies --output humann2_4normalize_final_out/humann2_genefamilies.tsv
#join matched genome humann2 output per sample into one table
   mkdir -p humann2_final_matched.genome.out
    runbioBakery humann2_join_tables -s --input ${baseDir}/genome_normalization/ --file_name pathabundance --output humann2_final_matched.genome.out/humann2_pathabundance.matchedGenome.tsv
    runbioBakery humann2_join_tables -s --input ${baseDir}/genome_normalization/ --file_name genefamilies --output humann2_final_matched.genome.out/humann2_genefamilies.matchedGenome.tsv
#Re-normalize RNASeq gene family and pathway abundances by corresponding DNA-level outputs to quantify microbial expression independent of gene copy number 
    runbioBakery humann2_rna_dna_norm -d humann2_final_matched.genome.out/humann2_pathabundance.matchedGenome.tsv -r humann2_4normalize_final_out/humann2_pathabundance.tsv -o humann2_pathabundance.DNAnormalizedRNA
    runbioBakery humann2_rna_dna_norm -d humann2_final_matched.genome.out/humann2_genefamilies.matchedGenome.tsv -r humann2_4normalize_final_out/humann2_genefamilies.tsv -o humann2_genefamilies.DNAnormalizedRNA
   
    """
}
/*
 * Step 14. HumanN2 merge and functional analysis transcriptome 
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
#reformat the humann2 log for multiqc
    for file in ${hn2rnaLog};do
        /uufs/chpc.utah.edu/common/home/u0762203/bin/microBiom/changePartFileName.rfLog.pl \$file
        rm \$file
    done

   
    """
}

/*
/*
 * Step 14. QC summary with multiqc for transcriptome
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

