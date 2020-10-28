#!/bin/bash

#SBATCH --account=hci-kp
#SBATCH --partition=hci-kp
#SBATCH -J q2pipe
#SBATCH --nodes=1
#SBATCH -e ./stderr.txt


# NOTES: Preprocesses (trim,import, merge, denoise, taxonomy) of 16S sequences using qiime2 from multiple sequencing lanes and merging based on sample if needed.
# Because this is to test and comapre different primers, trimming with different primers is done before importing to qiimme2 in a seperate script.
# Manifest files are made on the fly, but asssume all seq files in raw directories should be included.

#### User-defined variables #################################
# Raw sequences directory containing paired .fastq.gz files. List any number of directories in quotes between parentheses with space separation.
#RAWDIRS=( )
READSDIRS=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/18094primerTest/trimmedReads/
# Directory on scratch file system for intermediate files (will be created if not existing)
SCRATCH=/scratch/general/nfs1/u0762203/microbiome/16s/tmp
# The results directory for key outputs. Generally, your project directory or within it. (will be created if not existing)
ResultDir=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/18094primerTest/q2.result
# The qiime2 container image (update version tag as needed. Year.Month after "core:")
ContainerImage=docker://qiime2/core:2020.2
# In order to group samples with mulitple seq runs, and analyze each seq run, have maps listing by gnomexID and by sampleID and define column needed to group.
MAPBYGNOMEXID=
GROUPMAP=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/15215Rtest/sample-metadata.tsv
GROUPBYCOL=
##############################################################
#### Param (depends on the primer used)###################################################
Nproc=`nproc`
echo $Nproc
Nproc=$((Nproc-2))
echo $Nproc
JoinedTrimLength=TBD
# Default minimum size after merge is set to Rprimer+Fprimer+10. Permissive, but remove most primer dimers
MinMergeLength=200
MaxDiffs=30
# Classifier will likely need to be retrained with differnt qiime versions as scikit classifier different versions are not compatible.
CLASSIFIER=TBD
##############################################################
#### Setup ###################################################
module load singularity
mkdir -p ${SCRATCH}
mkdir -p ${SCRATCH}/tmp_XDG
mkdir -p ${SCRATCH}/tmp_sing
mkdir -p ${ResultDir}/q2_viz; mkdir -p ${ResultDir}/metadata
#############################################################
#### SINGULARITYENV (space & runtime issues) ################
XDG_RUNTIME_DIR=${SCRATCH}/tmp_XDG
SINGTEMPDIR=${SCRATCH}/tmp_sing
export SINGULARITYENV_XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR}
export SINGULARITYENV_TMPDIR=${SINGTEMPDIR}
#############################################################

echo "VERSION: SINGULARITY: `singularity --version`"
echo "VERSION: QIIME2: `singularity exec docker://qiime2/core:2020.2 qiime --version`"

# Part 1: Make manifest file on the fly, place in metadata file in project directory.
# Should allow for slightly changing naming schemes by core, but always assumes "_R1" somewhere in filename indicates read1 though.
echo "sample-id,absolute-filepath,direction" > ${ResultDir}/metadata/manifest_16SRawFiles.txt

for dir in ${READSDIRS[@]}
do
    cd ${dir}
    for f in *_R1.fq
    do
        SAMPLEID=`basename ${f%%_*}`
        echo "${SAMPLEID},${PWD}/${f},forward" >> ${ResultDir}/metadata/manifest_16SRawFiles.txt
        echo "${SAMPLEID},${PWD}/${f/_R1/_R2},reverse" >> ${ResultDir}/metadata/manifest_16SRawFiles.txt
    done
done
MANIFEST=${ResultDir}/metadata/manifest_16SRawFiles.txt

# Part 2: Import sequences (this is slow)
cd ${SCRATCH}

echo "TIME: START import = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ${MANIFEST} \
--output-path ${SCRATCH}/inseqs-demux.qza \
--input-format PairedEndFastqManifestPhred33

singularity exec ${ContainerImage} qiime demux summarize \
  --i-data ${SCRATCH}/inseqs-demux.qza \
  --o-visualization ${SCRATCH}/inseqs-demux-summary.qzv  
echo "TIME: END import = `date +"%Y-%m-%d %T"`"

# Part 3: Clean, denoise, filter chimeras, create table and phylogeny, call taxonomy

echo "TIME: START  merge = `date +"%Y-%m-%d %T"`"

singularity exec ${ContainerImage} qiime vsearch join-pairs \
--i-demultiplexed-seqs ${SCRATCH}/inseqs-demux.qza \
--o-joined-sequences PE-demux_trim_join.qza \
--p-minmergelen ${MinMergeLength} \
--verbose \

singularity exec ${ContainerImage} qiime demux summarize \
  --i-data PE-demux_trim_join.qza \
  --o-visualization PE-demux_trim_join.qzv
  
# At this stage plots should be inspected to infer the JoinedTrimLength for deblur. It should hold whenever using the same primer set though with good quality seq runs.

singularity exec ${ContainerImage} qiime quality-filter q-score-joined \
--i-demux PE-demux_trim_join.qza \
--o-filtered-sequences PE-demux_trim_join_filt.qza \
--o-filter-stats PE-demux_trim_join_filt_stats.qza \
--p-min-quality 10 
 
singularity exec ${ContainerImage} qiime demux summarize \
  --i-data PE-demux_trim_join_filt.qza \
  --o-visualization PE-demux_trim_join_filt.qzv

echo "TIME: END merge = `date +"%Y-%m-%d %T"`"
# Copy visualization results to ResultDir for inspection
cp ${SCRATCH}/*.qzv ${ResultDir}/q2_viz/
