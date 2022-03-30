#!/usr/bin/env nextflow

nextflow.enable.dsl=2 


/* 
 * Proof of concept Nextflow based humann3 shortgun metagenomic pipeline
 * 
 */ 

 
/*
 * Defines some parameters and specify 
 * read pairs,config file, metadata and outdir in job dir  by using the command line options
 */
//phixDb = "/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/FastQ_Screen_Genomes/PhiX/phi_plus_SNPs"
//hg38Db = "/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/FastQ_Screen_Genomes/Human/hg38.standard.bowtie2"
phixDbDir = "/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/FastQ_Screen_Genomes/PhiX"
hg38DbDir = "/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/FastQ_Screen_Genomes/Human"
proDb = "/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/db/uniref"
ntDb = "/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/db/chocophlan"
bt2Db = "/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/db/bowtie2db/mpa_v30_CHOCOPhlAn_201901"
//params.refStrain = "/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/refStrain/Bifidobacterium_breve/GCF_001281425.1_ASM128142v1_genomic.fna"
//chocFile = "/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/db/bowtie2db/mpa_v30_CHOCOPhlAn_201901/mpa_v30_CHOCOPhlAn_201901.pkl"
fqsDb = "/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/FastQ_Screen_Genomes"
log.info """
         Metagenome  P I P E L I N E:short-gun    
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

process dedup {
    
    input:
    tuple val(pair_id), path(read) 
    output:
    tuple val(pair_id), path("${pair_id}.clump*.fq.gz"), emit: fq 
    tuple val(pair_id), path("${pair_id}.clump1.fq.gz") , emit: fqc
 
    script:

    """
    #path of clumpify.sh
    /bbmap/clumpify.sh in1=${read[0]} in2=${read[1]} out1=${pair_id}.clump1.fq.gz out2=${pair_id}.clump2.fq.gz dupedist=12000 dedupe=t optical=t

    """
}
process trim {
    //stageInMode 'copy'
    //publishDir "${params.outdir}/knReadCounts", mode: 'copy', pattern: '*kneaddata_read_counts.phiX.human.txt'
    input:
    tuple val(pair_id), path("${pair_id}.clump*.fq.gz") 
    path(phixDb)
    path(hg38Db)
    output:
    tuple val(pair_id), path("${pair_id}.fq"), emit: fqh 
    tuple val(pair_id), path("${pair_id}.1.fq"), emit: fqs 
    path("./kneaddata_out/${pair_id}.clump1_kneaddata.log") , emit: report
    //path("${pair_id}.kneaddata_read_counts.phiX.human.txt")
 
    script:

    """
    #bowtie2 and Trimmomatic-0.36 need to be installed and know the path

    kneaddata -i ${pair_id}.clump1.fq.gz -i ${pair_id}.clump2.fq.gz -o kneaddata_out/ -db ${phixDb}/phi_plus_SNPs -db ${hg38Db}/hg38.standard.bowtie2 --bowtie2 /usr/bin/ --trimmomatic /Trimmomatic-0.36/ --bypass-trf -t ${task.cpus} --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" --bowtie2-options "--very-sensitive --dovetail"  
    kneaddata_read_count_table --input kneaddata_out --output ${pair_id}.kneaddata_read_counts.phiX.human.txt
    mkdir kneaddata_out/phix.human.contam_seq
    mv kneaddata_out/*_contam*fastq kneaddata_out/phix.human.contam_seq
    cat kneaddata_out/*_kneaddata_*_1.fastq > ${pair_id}.1.fq
    cat kneaddata_out/*_kneaddata_*fastq  > ${pair_id}.fq

    """
}

process fastqc {
    input:
    tuple val(pair_id), path("${pair_id}.clump1.fq.gz") 
    output:
    path("fastqc_${pair_id}_logs") 
 
    script:

    """
    mkdir fastqc_${pair_id}_logs
    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${pair_id}.clump1.fq.gz
    """
}

process fqscreenQc {
    input:
    tuple val(pair_id), path("${pair_id}.1.fq") 
    path(fqsConfig)
    path(fqsDb)
    output:
    path("${pair_id}.1_screen.txt") 
 
    script:

    """
    /FastQ-Screen-0.15.2/fastq_screen --conf ${fqsConfig} --aligner bowtie2 --threads ${task.cpus} ${pair_id}.1.fq
    """
}

process hmn3 {
    input:
    tuple val(pair_id), path("${pair_id}.fq") 
    path(proDb)
    path(ntDb)
    path(bt2Db)
    output:
     path("./humann3_out/${pair_id}_humann_temp/${pair_id}_metaphlan_bugs_list.tsv"), emit: bList  
     path("./humann3_out/${pair_id}_genefamilies_relab.tsv"), emit: gfr 
     path("./humann3_out/${pair_id}_pathabundance_relab.tsv"), emit: par 
     path("./humann3_out/${pair_id}_pathcoverage.tsv"), emit: pc 
    script:
 
    """
    mkdir -p humann3_out/
    humann --threads ${task.cpus} --input ${pair_id}.fq --output humann3_out/ --protein-database ${proDb} --nucleotide-database ${ntDb} --metaphlan-options "--bowtie2db ${bt2Db}" --prescreen-threshold 0.00001
    humann_renorm_table -i humann3_out/${pair_id}_genefamilies.tsv -o humann3_out/${pair_id}_genefamilies_relab.tsv -u relab
    humann_renorm_table -i humann3_out/${pair_id}_pathabundance.tsv -o humann3_out/${pair_id}_pathabundance_relab.tsv -u relab
    """
}


process strainPh {
    input:
    tuple val(pair_id), path("${pair_id}.fq")
    path(bt2Db)
    output:
    path("strainPh/${pair_id}.sam.bz2")
    script:

    """
    mkdir -p strainPh
    metaphlan ${pair_id}.fq --input_type fastq --bowtie2db ${bt2Db} -s strainPh/${pair_id}.sam.bz2 --bowtie2out strainPh/${pair_id}.bowtie2.bz2 -o strainPh/${pair_id}_profile.tsv


    """
}
process mergeHmn {
    publishDir "${params.outdir}", mode: 'copy' 
    input:
    path(bugsListTsv) 
    path(genefamTsv) 
    path(pathabundTsv) 
    path(pathcovTsv) 
    output:
    path("humann*tsv")
    path("humann*spf")
    path("metaphlan_merged.spf")

    script:

    """
    mkdir hmn3
    cp ${pathabundTsv} hmn3/
    cp ${pathcovTsv} hmn3/
    cp ${genefamTsv} hmn3/
    humann_join_tables -s -i hmn3/ --file_name pathabundance_relab -o humann_pathabundance.tsv
    humann_join_tables -s -i hmn3/ --file_name pathcoverage -o humann_pathcoverage.tsv
    humann_join_tables -s -i hmn3/ --file_name genefamilies_relab -o humann_genefamilies.tsv

    humann_split_stratified_table -i humann_pathabundance.tsv -o .
    humann_split_stratified_table -i humann_genefamilies.tsv -o .
    humann_split_stratified_table -i humann_pathcoverage.tsv -o .

    sed 's/_Abundance-RPKs//g' humann_genefamilies_unstratified.tsv | sed 's/# Gene Family/GeneFamily/' | sed 's/_kneaddata//g' > humann_genefamilies_unstratified.spf
    sed 's/_Abundance//g' humann_pathabundance_unstratified.tsv | sed 's/# Pathway/Pathway/' | sed 's/_kneaddata//g' > humann_pathabundance_unstratified.spf

    merge_metaphlan_tables.py ${bugsListTsv} > metaphlan_merged.txt
    sed -i 's/_metaphlan_bugs_list//g' metaphlan_merged.txt
    #scp metaphlan_to_stamp.pl to zion to bind into container...
    metaphlan_to_stamp.pl metaphlan_merged.txt | sed 's/_kneaddata//g' > metaphlan_merged.spf
    


    """
}



process mergeStrain {
    publishDir "${params.outdir}", mode: 'copy' 
    input:
    path(strainPhSamFiles)
    path(bt2Db) 
    path(refStrain)
    path(metadata)
    
    output:
    path("sams_simplified/*.sam.bz2")
    //path("consensus_markers/*pkl")
    //path("strainPhOutput/RAxML*")
    //path("strainPhOutput/*info")
    //path("strainPhOutput/*polymorphic")
    //path("strainPhOutput/*aln")
    script:

    """
    
    mkdir -p sams_simplified
    cp ${strainPhSamFiles} sams_simplified/
    #mkdir -p consensus_markers
    #sample2markers.py -i sams_simplified/*.sam.bz2 -o consensus_markers -n ${task.cpus}
    #mkdir -p clade_markers
    #mkdir -p strainPhOutput
    #clade can be changed according to request
    #extract_markers.py -c s__Eubacterium_rectale -d ${bt2Db}/mpa_v30_CHOCOPhlAn_201901.pkl -o clade_markers
    #strainphlan -s consensus_markers/*.pkl -m clade_markers/s__Eubacterium_rectale.fna -r ${refStrain} -o strainPhOutput -c s__Eubacterium_rectale --phylophlan_mode fast -n ${task.cpus} -d ${bt2Db}/mpa_v30_CHOCOPhlAn_201901.pkl
    #cd strainPhOutput/
    #add_metadata_tree.py -t strainPhOutput/RAxML_bestTree.s__Eubacterium_rectale.StrainPhlAn3.tre -f ${metadata}
    #plot_tree_graphlan.py -t strainPhOutput/RAxML_bestTree.s__Eubacterium_rectale.StrainPhlAn3.tre.metadata -m indv --leaf_marker_size 60 --legend_marker_size 60


    """
}

process multiqc {
    publishDir "${params.outdir}", mode:'copy'

    input:
    path mqcFiles

    output:
    file('multiqc_report.html') optional true

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    multiqc -v .

    """
}

workflow {
    dedup(read_pairs)
    trim(dedup.out.fq, phixDbDir, hg38DbDir)
    fastqc(dedup.out.fqc)
    fqscreenQc(trim.out.fqs, params.fqsConfig, fqsDb)
    hmn3(trim.out.fqh, proDb, ntDb, bt2Db)
    //strainPh(trim.out.fqh, bt2Db)
    mergeHmn(hmn3.out.bList.collect(), hmn3.out.gfr.collect(), hmn3.out.par.collect(), hmn3.out.pc.collect())
    //mergeStrain(strainPh.out.collect(), bt2Db, params.refStrain, params.metadata)
    multiqc(trim.out.report.mix(fastqc.out).mix(fqscreenQc.out).collect())
}




