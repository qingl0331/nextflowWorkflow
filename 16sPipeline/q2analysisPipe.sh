#!/bin/bash

#SBATCH --account=hci-kp
#SBATCH --partition=hci-kp
#SBATCH -J q2pipe
#SBATCH --nodes=1
#SBATCH -e ./stderr.txt


# NOTES: 16S analysis pipeline- Alpha and beta diversity analyses, Differential abundance measurements, xxx.
# Take the output (table_bySeqID.qza, tree_root.qza) from preprocessing pipeline q2pipe16sMultipleLanes.sh as input.

#### User-defined variables #################################
table=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/15215Rtest/q2.bl.fake.result/table_byColoCareID.qza
tax=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/15215Rtest/q2.bl.fake.result/taxonomy.qza
rootedTree=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/15215Rtest/q2.bl.fake.result/tree_root.qza
sampMetadata=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/15215Rtest/fake.map_SampleID_2_ColocareID.txt
SCRATCH=/scratch/general/nfs1/u0762203/microbiome/16s/tmp
## Need inspect table_bySeqID.qzv to decide the  --p-sampling-depth 
pSampDep=11144 
pMaxDep=28000 
# The qiime2 container image (update version tag as needed. Year.Month after "core:")
ContainerImage=docker://qiime2/core:2020.2
##############################################################
##############################################################
#### Setup ###################################################
module load singularity
mkdir -p ${SCRATCH}
mkdir -p ${SCRATCH}/tmp_XDG
mkdir -p ${SCRATCH}/tmp_sing
#############################################################
#### SINGULARITYENV (space & runtime issues) ################
XDG_RUNTIME_DIR=${SCRATCH}/tmp_XDG
SINGTEMPDIR=${SCRATCH}/tmp_sing
export SINGULARITYENV_XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR}
export SINGULARITYENV_TMPDIR=${SINGTEMPDIR}
#############################################################

echo "VERSION: SINGULARITY: `singularity --version`"
echo "VERSION: QIIME2: `singularity exec docker://qiime2/core:2020.2 qiime --version`"


# Part 1: alpha and beta diversity analysis
#computing diversity metrics
: '
singularity exec ${ContainerImage} qiime diversity core-metrics-phylogenetic \
--i-phylogeny $rootedTree \
--i-table $table \
--p-sampling-depth $pSampDep \
--m-metadata-file $sampMetadata \
--output-dir core-metrics-results
# '

#test for associations between categorical metadata columns and alpha diversity data.Here we use Faith Phylogenetic Diversity (a measure of community richness) and evenness metrics. For continuous sample metadata columns, use qiime diversity alpha-correlation.
: '
singularity exec ${ContainerImage} qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
--m-metadata-file $sampMetadata \
--o-visualization core-metrics-results/faith-pd-group-significance.qzv


singularity exec ${ContainerImage} qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/evenness_vector.qza \
--m-metadata-file $sampMetadata \
--o-visualization core-metrics-results/evenness-group-significance.qzv
# '
#test for associations between categorical metadata columns and beta diversity data using PERMANOVA (apply this to our unweighted UniFrac distances). It will test whether distances between samples within a group, such as samples from the same timepoint (e.g.,baseline), are more similar to each other than they are to samples from the other groups. For continous metadata, use qiime diversity mantel and qiime diversity bioenv. 
: '
singularity exec ${ContainerImage} qiime diversity beta-group-significance \
--i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file $sampMetadata \
--m-metadata-column Timepoint \
--o-visualization core-metrics-results/unweighted-unifrac-timepoint-significance.qzv \
--p-pairwise
# '
#This method cannot operate on a grouping vector with only unique values
: '
singularity exec ${ContainerImage} qiime diversity beta-group-significance \
--i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file $sampMetadata \
--m-metadata-column PrimerID \
--o-visualization core-metrics-results/unweighted-unifrac-primer-significance.qzv \
--p-pairwise
# '
#use the Emperor tool to explore principal coordinates (PCoA) plots in the context of sample metadata. With --p-custom-axes, which is very useful for exploring time series data. See MP example.

 
#explore alpha diversity as a function of sampling depth using alpha-rarefaction visualizer
: '
singularity exec ${ContainerImage} qiime diversity alpha-rarefaction --i-table $table --i-phylogeny $rootedTree --p-max-depth $pMaxDep --m-metadata-file $sampMetadata --o-visualization alpha-rarefaction.qzv
# '

#Taxonomic analysis already in pre-process pipeline


#Differential abundance testing with ANCOM.ANCOM assumes that few (less than about 25%) of the features are changing between groups.--> filter (on body site)before run it (see command on evernote).

: '
singularity exec ${ContainerImage} qiime composition add-pseudocount \
--i-table $table \
--o-composition-table core-metrics-results/comp-table.qza

singularity exec ${ContainerImage} qiime composition ancom \
--i-table core-metrics-results/comp-table.qza \
--m-metadata-file $sampMetadata \
--m-metadata-column Timepoint \
--o-visualization core-metrics-results/ancom-timepoint.qzv

# '

#performing a differential abundance test at a specific taxonomic level, for example, collapse our feature table at the family level (i.e. level 5 of the Greengenes taxonomy)

#: '
singularity exec ${ContainerImage} qiime taxa collapse \
--i-table $table \
--i-taxonomy $tax \
--p-level 5 \
--o-collapsed-table core-metrics-results/timepoint-table-l5.qza  

singularity exec ${ContainerImage} qiime composition add-pseudocount \
--i-table core-metrics-results/timepoint-table-l5.qza \
--o-composition-table core-metrics-results/comp-tp-table-l5.qza

singularity exec ${ContainerImage} qiime composition ancom \
--i-table core-metrics-results/comp-tp-table-l5.qza \
--m-metadata-file $sampMetadata \
--m-metadata-column Timepoint \
--o-visualization core-metrics-results/l5-ancom-timepoint.qzv

# '
