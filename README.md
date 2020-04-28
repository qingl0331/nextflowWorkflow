# nextflow workflows
1. The fqScreen.nf is a QC pipeline for paired end reads for exome, genome or transcriptome sequencing. It includes steps of deduplication, adapter and quality trimming, QC with fastqc for reads quality, QC with fastqscreen for sample composition and removal of contaminated reads (the example is for human reads, but you can update the databases in the fastqscreen config file). 
2. The metagenomePipe.nf is a metagenomic pipeline for paired end reads in metagenome sequencing. Currently, it mainly includes QC steps and HUMAnN2 classification and functional analysis.
3. The metatranscriptomePipe.nf is a metatranscriptome pipeline for paired end reads in metatranscriptome sequencing. Please note that it will run the metagenome analysis and use the output to normalize the paired metatranscriptome. Accessory scripts are provided in the accessoryScripts folder.  
4. The loh.nf is a clinical genomic pipeline that detects LOH event from paired normal and tumor samples. The accessory scripts are also provided in accessoryScripts. 

Both the loh and microbiome pipelines have a singularity version for slurm systems in the corresponding folders.  To run the microbiome pipeline, you shoud 
* download the scripts in accessoryScripts and microbiomeSingularityPipeline sub-directories
* run: nextflow run metagnom.master.nf -with-trace -resume --reads <path to Fastqs> --outdir <output directory> 
* or nextflow run metatrans.master.nf -with-trace -resume --reads "xx/xx" --rnaReads "xx/xx" --outdir "xx/xx" --rnaOutdir "xx/xx" --sampSheet "xx/xx"
