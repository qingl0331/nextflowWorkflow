# Setting up jobs
Download the accessory scripts and put them in your path;
Download the scripts/sample sheets in microbiomeSingularityPipeline and put all the scripts in your working directory;
Type in your terminal (better run it in a screen or with nohup):               
    $nextflow run metagnom.master.nf -with-trace -resume --reads "path to Fastqs/*_R{1,2}*.fastq.gz" --outdir "output director"             
    $nextflow run metatrans.master.nf -with-trace -resume --reads "path to Fastqs/*_R{1,2}*.fastq.gz" --rnaReads "path to RNASeq Fastqs/*_R{1,2}*.fastq.gz" --outdir "output directory for genome" --rnaOutdir "output directory for transcriptome" --sampSheet "xx/xx"
