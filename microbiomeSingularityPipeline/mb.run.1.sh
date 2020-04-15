#!/bin/bash
set -e; start=$(date +'%s');  touch STARTED




#### Do just once ####

# 1) Install Singularity (https://www.sylabs.io) or load via a module then define the path to the executable
module load singularity/3.5.2
singExec=/uufs/chpc.utah.edu/sys/installdir/singularity3/std/bin/singularity

# 2) Define file paths to "mount" in the container.
g38Dir=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/Human/GRCh38/bowtie2
runDir=/scratch/general/lustre/u0762203/microbiomePipeSingu/run
hn2pdb=/uufs/chpc.utah.edu/sys/installdir/bioBakery/databases/uniref
hn2ndb=/uufs/chpc.utah.edu/sys/installdir/bioBakery/databases/chocophlan
srRnaDb=/uufs/chpc.utah.edu/sys/installdir/sortmerna/sortmerna-2.1b/rRNA_databases
# 3) Modify the workflow xxx.sing file setting the paths to the required resources. These must be within the mounts.

# 4) Build the singularity container, and define the path to the xxx.sif file, do just once after each update.
#$singExec pull docker://qingl0331/qinglhci2019:microbiomePipe
container=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/qinglhci2019_microbiomePipe.sif


#### Do for every run ####

# 1) Create a folder named as you would like the analysis name to appear. This must reside in the mount paths.

# 2) Copy over the workflow docs: xxx.sing, xxx.sh into the job directory.
# 3) Launch the xxx.sh via slurm's sbatch or run it on your local server.  




#### No need to modify anything below ####

echo -e "\n---------- Starting -------- $((($(date +'%s') - $start)/60)) min"

# Read out params
#name=${PWD##*/} # last dir of current dir, which is single patient dir
jobDir=`readlink -f .` #full path of current dir
SINGULARITYENV_jobDir=$jobDir  SINGULARITYENV_g38Dir=$g38Dir SINGULARITYENV_runDir=$runDir SINGULARITYENV_hn2pdb=$hn2pdb SINGULARITYENV_hn2ndb=$hn2ndb SINGULARITYENV_srRnaDb=$srRnaDb SINGULARITYENV_container=$container $singExec exec --containall --bind $g38Dir $container bash $runDir/mb.1.sing

echo -e "\n---------- Complete! -------- $((($(date +'%s') - $start)/60)) min total"

touch COMPLETE 

