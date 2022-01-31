#!/bin/bash

#SBATCH --account=hci-kp
#SBATCH --partition=hci-kp
#SBATCH -J q2prPipe
#SBATCH --nodes=1
#SBATCH -e ./stderr.txt


# NOTES: Preprocesses (repair,import, check parameter for dada2) of 16S sequences using qiime2 from multiple sequencing lanes and merging based on sample if needed.
# Because this uses commercial V3V4 primers, trimming with  primers and adapters are done before importing to qiimme2 
# Manifest files are made on the fly.

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
fastq_screen=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app/fastq_screen/v0.14.0/fastq_screen
config=./fastq_screen.conf
repairSh=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app/BBmap/v38.90/repair.sh

# Classifier will likely need to be retrained with differnt qiime versions as scikit classifier different versions are not compatible.
CLASSIFIER=/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/feature-classifiers/classifier.silva.2011.11.qza



# Part 1: Pairing filtered FASTQ files
echo "TIME: START repair = `date +"%Y-%m-%d %T"`"
mkdir -p repairFilt
EXTENSION=.tagged_filter.fastq
ls ${READSDIR}/*$EXTENSION | sed -e "s/R[12]$EXTENSION//" -e 's/.*\///g' | uniq | parallel --no-notice -j $Nproc $repairSh \
 in1=${READSDIR}/{}R1$EXTENSION \
 in2=${READSDIR}/{}R2$EXTENSION \
 out1=repairFilt/{}repr_R1.fq out2=repairFilt/{}repr_R2.fq outs=repairFilt/{}repr_orph.fq repair


echo "TIME: END repair = `date +"%Y-%m-%d %T"`"

# Part 2: QC with fastqscreen after repair.
mkdir -p fqscreenQCafFtRp
module load bowtie2

echo "TIME: START  QC = `date +"%Y-%m-%d %T"`"
for f in repairFilt/*repr_R1.fq
do
	$fastq_screen --conf ${config} --threads $Nproc --outdir fqscreenQCafFtRp $f
done


#multiqc
module load multiqc
multiqc .
module unload multiqc

echo "TIME: END QC  = `date +"%Y-%m-%d %T"`"
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

# Part 3: Make manifest file on the fly, place in metadata file in project directory.
# Should allow for slightly changing naming schemes by core, but always assumes "_R1" somewhere in filename indicates read1 though.
echo "sample-id,absolute-filepath,direction" > ${ResultDir}/metadata/manifest_16SRawFiles.txt

cd ./repairFilt
for f in *_R1.fq
do
	SAMPLEID=`basename ${f%%_*}`
	echo "${SAMPLEID},${PWD}/${f},forward" >> ${ResultDir}/metadata/manifest_16SRawFiles.txt
	echo "${SAMPLEID},${PWD}/${f/_R1/_R2},reverse" >> ${ResultDir}/metadata/manifest_16SRawFiles.txt
done
MANIFEST=${ResultDir}/metadata/manifest_16SRawFiles.txt

# Part 4: Import sequences (this is slow)
cd ../

echo "TIME: START import = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ${MANIFEST} \
--output-path ${SCRATCH}/inseqs-demux.qza \
--input-format PairedEndFastqManifestPhred33
# At this stage plots should be inspected to infer the --p-trunc-len-f, --p-trunc-len-r, --p-trim-left-f and --p-trim-left-r for dada2. It should hold whenever using the same primer set though with good quality seq runs.
singularity exec ${ContainerImage} qiime demux summarize \
  --i-data ${SCRATCH}/inseqs-demux.qza \
  --o-visualization ${SCRATCH}/inseqs-demux-summary.qzv  
echo "TIME: END import = `date +"%Y-%m-%d %T"`"

# Copy visualization results to ResultDir for inspection
cp ${SCRATCH}/*.qzv ${ResultDir}/q2_viz/

#dada2 QC parameters:Position at which forward/reverse read sequences should be truncated due to decrease in quality, need to be adjusted according to inspection of inseqs-demux-summary.qza (shoot for middle box 20 - 30 maybe)

