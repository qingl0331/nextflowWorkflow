# metagenome pipeline
1. The metagenome pipeline is humann3 based metagenome short-gun seq data processing pipeline.
2. The pipeline is implemented with nextflow, kneaddata and fastqscreen, etc for quick data processing, data cleaning and QC. 
3. All the software packages needed for the pipeline is installed in a container, so you can directly run it locally or on any HPC platform.

    

* run on chpc:
* For example, in screen:
>
$module load nextflow
>
$nextflow run mg.master.singu.dsl2.hn3.nf --reads "/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/metagenome/raw_data/*_R{1,2}*.fastq.gz" --fqsConfig "/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/metagenome/fastq_screen.conf" --outdir "./results" -resume
>
