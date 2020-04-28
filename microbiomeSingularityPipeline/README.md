# Setting up jobs
Download the accessory scripts and put them in your path;
Download the scripts/sample sheets in microbiomeSingularityPipeline and put all the scripts in your working directory;
Type in your terminal (better run it in a screen or with nohup):               
    $nextflow run metagnom.master.nf -with-trace -resume --reads "xx/xx" --outdir "xx/xx"             
    $nextflow run metatrans.master.nf -with-trace -resume --reads "xx/xx" --rnaReads "xx/xx" --outdir "xx/xx" --rnaOutdir "xx/xx" --sampSheet "xx/xx"
