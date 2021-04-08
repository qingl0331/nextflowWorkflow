# nextflow workflows
1. The fqScreen.nf is a QC pipeline for paired end reads for exome, genome or transcriptome sequencing. It includes steps of deduplication, adapter and quality trimming, QC with fastqc for reads quality, QC with fastqscreen for sample composition and removal of contaminated reads (the example is for human reads, but you can update the databases in the fastqscreen config file). 
2. The metagenomePipe.nf is a metagenomic pipeline for paired end reads in metagenome sequencing. Currently, it mainly includes QC steps and HUMAnN2 classification and functional analysis.
3. The metatranscriptomePipe.nf is a metatranscriptome pipeline for paired end reads in metatranscriptome sequencing. Please note that it will run the metagenome analysis and use the output to normalize the paired metatranscriptome. Accessory scripts are provided in the accessoryScripts folder.  
4. The loh.nf is a clinical genomic pipeline that detects LOH event from paired normal and tumor samples. The accessory scripts are also provided in accessoryScripts.
5. The 16sPipeline folder contains the QIIME2 based 16S rRNA data processing and analysis pipelines.
6. The bacteriaWGS_alignVar.nf is a GATK4 based  all-inclusive bacteria variant call pipeline, which includes alignment, QC, variant call and variant annotation. 
7. The bacteriaWGS_alignVar.dsl2.nf is a DSL2 incorporated, modulized version for bacteriaWGS_alignVar.nf, compatible with nextflow version 20.10 or latter. You need execute it with -profile slurm option. 
    

Both the loh and microbiome pipelines have a singularity version for slurm systems in the corresponding folders.  To run the microbiome pipeline, you shoud 
* download the scripts in accessoryScripts sub-directory and put them in your path; Download the scripts/sample sheets in microbiomeSingularityPipeline sub-directory and put them in your working directory;
* run: nextflow run metagnom.master.nf -with-trace -resume --reads "path to Fastqs/*_R{1,2}*.fastq.gz" --outdir "output directory" 
* or nextflow run metatrans.master.nf -with-trace -resume --reads "path to Fastqs/*_R{1,2}*.fastq.gz" --rnaReads "path to RNASeq Fastqs/*_R{1,2}*.fastq.gz" --outdir "output directory for genome" --rnaOutdir "output directory for transcriptome" --sampSheet "xx/xx"
