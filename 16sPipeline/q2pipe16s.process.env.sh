#!/bin/bash

#SBATCH --account=hci-kp
#SBATCH --partition=hci-kp
#SBATCH -J q2prPipe
#SBATCH --nodes=1
#SBATCH -e ./stderr.txt


# Because this uses commercial V3V4 primers, trimming with  primers and adapters are done before importing to qiimme2 

#### User-defined variables #################################
workDir=${PWD}
READSDIR=$workDir/results/fqscreenFilt/
# Directory on scratch file system for intermediate files (will be created if not existing)
SCRATCH=$workDir/tmp
# The results directory for key outputs. Generally, your project directory or within it. (will be created if not existing)
ResultDir=$workDir/q2.dada2.result
# The qiime2 container image (update version tag as needed. Year.Month after "core:")
ContainerImage=/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/containers/qimme2.core_2021.11.sif
MAPBYGNOMEXID=$workDir/sampleSheet.txt
##############################################################
Nproc=`nproc`
echo $Nproc
Nproc=$((Nproc-2))
echo $Nproc

#dada2 QC parameters:Position at which forward/reverse read sequences should be truncated due to decrease in quality, need to be adjusted according to inspection of inseqs-demux-summary.qzv 
ptf=275
ptr=267


# Classifier will likely need to be retrained with differnt qiime versions as scikit classifier different versions are not compatible.
CLASSIFIER=/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/feature-classifiers/classifier.silva.2011.11.qza



##############################################################
#### Setup ###################################################
module load singularity
#############################################################
#### SINGULARITYENV (space & runtime issues) ################
XDG_RUNTIME_DIR=${SCRATCH}/tmp_XDG
SINGTEMPDIR=${SCRATCH}/tmp_sing
export SINGULARITYENV_XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR}
export SINGULARITYENV_TMPDIR=${SINGTEMPDIR}
#############################################################
echo "VERSION: SINGULARITY: `singularity --version`"


#  Clean, denoise, create table and phylogeny, call taxonomy
echo "TIME: START  dada2 denoise = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime dada2 denoise-paired \
--i-demultiplexed-seqs ${SCRATCH}/inseqs-demux.qza \
--p-trunc-len-f $ptf \
--p-trunc-len-r $ptr \
--p-n-threads $Nproc \
--o-table ${SCRATCH}/table_bySeqID.qza \
--o-representative-sequences ${SCRATCH}/repseq.qza \
--o-denoising-stats ${SCRATCH}/table_bySeqID_stats.qza

singularity exec ${ContainerImage} qiime feature-table summarize \
--i-table ${SCRATCH}/table_bySeqID.qza \
--o-visualization ${SCRATCH}/table_bySeqID.qzv \
--m-sample-metadata-file ${MAPBYGNOMEXID}


echo "TIME: END dada2 denoise = `date +"%Y-%m-%d %T"`"

singularity exec ${ContainerImage} qiime feature-table tabulate-seqs \
--i-data ${SCRATCH}/repseq.qza \
--o-visualization ${SCRATCH}/repseq.qzv



echo "TIME: START phylogeny = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences ${SCRATCH}/repseq.qza \
--o-alignment ${SCRATCH}/aligned_repseq.qza \
--o-masked-alignment ${SCRATCH}/masked_aligned_repseq.qza \
--o-tree ${SCRATCH}/tree_unroot.qza \
--o-rooted-tree ${SCRATCH}/tree_root.qza

echo "TIME: END phylogeny = `date +"%Y-%m-%d %T"`"

echo "TIME: START taxonomy = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime feature-classifier classify-sklearn \
--i-classifier ${CLASSIFIER} \
--i-reads ${SCRATCH}/repseq.qza \
--o-classification ${SCRATCH}/taxonomy.qza \
--p-n-jobs $Nproc

singularity exec ${ContainerImage} qiime metadata tabulate \
--m-input-file ${SCRATCH}/taxonomy.qza \
--o-visualization ${SCRATCH}/taxonomy.qzv

singularity exec ${ContainerImage} qiime taxa barplot \
--i-table ${SCRATCH}/table_bySeqID.qza \
--i-taxonomy ${SCRATCH}/taxonomy.qza \
--m-metadata-file ${MAPBYGNOMEXID} \
--o-visualization ${SCRATCH}/taxbarplots_AllSampNoFilter.qzv


echo "TIME: END taxonomy = `date +"%Y-%m-%d %T"`"


# Copy visualization results to ResultDir for inspection
cp ${SCRATCH}/*.qzv ${ResultDir}/q2_viz/
cp ${SCRATCH}/*able_bySeqID.qza ${ResultDir}/
cp ${SCRATCH}/*epseq.qza ${ResultDir}/
cp ${SCRATCH}/*axonomy.qza ${ResultDir}/
cp ${SCRATCH}/*ree_root.qza ${ResultDir}/
cp ${SCRATCH}/*ree_unroot.qza ${ResultDir}/

