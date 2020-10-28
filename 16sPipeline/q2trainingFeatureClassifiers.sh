#!/bin/bash

#SBATCH --account=hci-kp
#SBATCH --partition=hci-kp
#SBATCH -J q2tc
#SBATCH --nodes=1
#SBATCH -e ./stderr.txt


# NOTES: Training feature classifiers for Greengenes (16S rRNA), Silva (16S/18S rRNA), UNITE (fungal ITS).

#### User-defined variables #################################
ggRef=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/analysisPipDev/training-feature-classifiers/markerGeneRefDb/gg_13_8_otus/rep_set/99_otus.fasta
ggTax=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/analysisPipDev/training-feature-classifiers/markerGeneRefDb/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt

silvaRef=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/analysisPipDev/training-feature-classifiers/markerGeneRefDb/SILVA_132_QIIME_release/rep_set/rep_set_all/99/silva132_99.fna
silvaTax=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/analysisPipDev/training-feature-classifiers/markerGeneRefDb/SILVA_132_QIIME_release/taxonomy/taxonomy_all/99/taxonomy_all_levels.txt


unitRef=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/analysisPipDev/training-feature-classifiers/markerGeneRefDb/sh_qiime_release_s_04.02.2020/developer/sh_refs_qiime_ver8_dynamic_s_04.02.2020_dev.uc.fasta
unitTax=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/analysisPipDev/training-feature-classifiers/markerGeneRefDb/sh_qiime_release_s_04.02.2020/developer/sh_taxonomy_qiime_ver8_dynamic_s_04.02.2020_dev.txt

repSeq=/uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/neli/microbiome/16SrRNA/15215Rtest/q2.bl.fake.result/repseq.qza

FPrimerSeq=TGCCTACGGGNBGCASCAG
RPrimerSeq=GCGACTACNVGGGTATCTAATCC
MinMergeLength=189
MaxMergeLength=600

SCRATCH=/scratch/general/nfs1/u0762203/microbiome/16s/tmp


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


#import gg data into QIIME 2 Artifacts
: '
singularity exec ${ContainerImage} qiime tools import \
--type ''FeatureData[Sequence]'' \ #need remove extra quotes
--input-path $ggRef \
--output-path ggOutput/99_otus.qza
# '

: '
singularity exec ${ContainerImage} qiime tools import \
--type ''FeatureData[Taxonomy]'' \ #need remove extra quotes
--input-format HeaderlessTSVTaxonomyFormat \
--input-path $ggTax \
--output-path ggOutput/ref-taxonomy.qza
# '


#Extract reference reads (do it only for bacterial 16S rRNA)
: '
singularity exec ${ContainerImage} qiime feature-classifier extract-reads \
--i-sequences ggOutput/99_otus.qza \
--p-f-primer $FPrimerSeq \
--p-r-primer $RPrimerSeq \
--p-min-length $MinMergeLength \
--p-max-length $MaxMergeLength \
--o-reads ggOutput/ref-seqs.qza
# '
#Train the classifier for primer test - don't extract ref reads but use 99_otus.qza directly

#: '
singularity exec ${ContainerImage} qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ggOutput/99_otus.qza \
--i-reference-taxonomy ggOutput/ref-taxonomy.qza \
--o-classifier ggOutput/pTclassifier.qza
# '

#Train the classifier

: '
singularity exec ${ContainerImage} qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ggOutput/ref-seqs.qza \
--i-reference-taxonomy ggOutput/ref-taxonomy.qza \
--o-classifier ggOutput/classifier.qza
# '

#test the classifier 
: '
singularity exec ${ContainerImage} qiime feature-classifier classify-sklearn \
--i-classifier ggOutput/classifier.qza \
--i-reads $repSeq \
--o-classification ggOutput/taxonomy.qza

singularity exec ${ContainerImage} qiime metadata tabulate \
--m-input-file ggOutput/taxonomy.qza \
--o-visualization ggOutput/taxonomy.qzv

# '

#Fit a classifier for the UNITE Database (fungal ITS)
: '
singularity exec ${ContainerImage} qiime tools import \
--type ''FeatureData[Sequence]'' \  #need remove extra quotes
--input-path $unitRef \
--output-path unitOutput/unite-ver8-dynamic-dev-seqs.qza


singularity exec ${ContainerImage} qiime tools import \
--type ''FeatureData[Taxonomy]'' --input-format HeaderlessTSVTaxonomyFormat \  #need remove extra quotes
--input-path $unitTax \
--output-path unitOutput/unite-ver8-dynamic-dev-tax.qza

singularity exec ${ContainerImage} qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads unitOutput/unite-ver8-dynamic-dev-seqs.qza \
--i-reference-taxonomy unitOutput/unite-ver8-dynamic-dev-tax.qza \
--o-classifier unitOutput/unite-ver8-dynamic-dev-classifier.qza 
# '


