#!/bin/bash
#SBATCH --account=hci-rw
#SBATCH --partition=hci-rw
#SBATCH -N 1
#SBATCH -t 24:00:00

set -e; start=$(date +'%s'); rm -f FAILED COMPLETE QUEUED; touch STARTED

# Huntsman Cancer Institute
# This fires a Sentieon CNV workflow generating germline Copy Number Variant calls. See https://support.sentieon.com/appnotes/cnv/ 


# This uses the fast Sentieon apps where appropriate.  Uses modules loaded on redwood.chpc.utah.edu

# 1) Install and load sentieon, bedtools, and snakemake as modules
module use /uufs/chpc.utah.edu/common/PE/proj_UCGD/modulefiles/$UUFSCELL &> /dev/null
module load sentieon  &> /dev/null
module load bedtools  &> /dev/null
module load nextflow/20.10


# 2) Define file paths to the TNRunner data bundle downloaded and uncompressed from https://hci-bio-app.hci.utah.edu/gnomex/?analysisNumber=A5578 . The second is the path to your data.
dataBundle=/uufs/chpc.utah.edu/common/PE/hci-bioinformatics1/TNRunner
myData=/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/qli
fasta=$dataBundle/Indexes/B38IndexForBwa-0.7.17/hs38DH.fa
avatarPon=$myData/cnvBenchmark/avatarPon
exomBed=$myData/cnvBenchmark/hg38NimIdtMergedPad150bp.bed
ucscAno=$myData/db/hg38RefSeq8Aug2018_Merged.ucsc.gz

# 3) Create a folder named as you would like the analysis name to appear, this along with the genome build will be prepended onto all files, no spaces, change into it. 

# 4) Soft link your: normal bam file and its associated index into the job dir.

# 5) Copy over the workflow docs:  xxx.README.sh, and xxx.nf into the job directory.

# 6) Launch the xxx.README.sh via sbatch or run it on your local server.  

# 7) If the run fails, fix the issue and restart.  Nextflow should pick up where it left off.



#### No need to modify anything below ####

echo -e "\n---------- Starting -------- $((($(date +'%s') - $start)/60)) min"

# Read out params
name=${PWD}

unset OMP_NUM_THREADS
allThreads=`nproc`
allRam=$(expr `free -g | grep -oP '\d+' | head -n 1` - 2)


# Print params
echo
echo -n name"         : "; echo $name
echo -n threads"      : "; echo $allThreads
echo -n ram"          : "; echo $allRam
echo -n host"         : "; echo $(hostname)

# Launch nextflow


bam=$(ls $name/*bam)
bamIdx=$(ls $name/*bai)
nextflow run copyAnalysis.nf -with-trace -resume --bam "$bam" --bamIdx "$bamIdx" --fasta "$fasta" --avatarPon "$avatarPon" --exomBed "$exomBed" --ucscAno "$ucscAno" 


echo -e "\n---------- Complete! -------- $((($(date +'%s') - $start)/60)) min total"

# Final cleanup
mkdir -p RunScripts
mv -f copyAnalysis* RunScripts/
mv -f slurm* Logs/ || true
rm -f FAILED STARTED DONE RESTART*
touch COMPLETE 

