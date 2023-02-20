#!/bin/bash

#SBATCH --account=hci-kp
#SBATCH --partition=hci-kp
#SBATCH -J q2prPipe
#SBATCH --nodes=1
#SBATCH -e ./stderr.1.txt


# NOTES: Preprocesses (repair,import, check parameter for dada2) of 16S sequences using qiime2 from multiple sequencing lanes.

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
config=/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/16s/fastq_screen.conf
repairSh=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app/BBmap/v38.90/repair.sh

#### Param (depends on the primer used)###################################################
FPrimerSeqToTrim=GTGCCAGCMGCCGCGGTAA
RPrimerSeqToTrim=GGACTACHVGGGTWTCTAAT


# Classifier will likely need to be retrained with differnt qiime versions as scikit classifier different versions are not compatible.
CLASSIFIER=/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/feature-classifiers/classifier.2011.11.qza



# Part 1: Pairing filtered FASTQ files
echo "TIME: START repair = `date +"%Y-%m-%d %T"`"
mkdir -p repairFilt
EXTENSION=.tagged_filter.fastq.gz
ls ${READSDIR}/*$EXTENSION | sed -e "s/R[12]$EXTENSION//" -e 's/.*\///g' | uniq | parallel --no-notice -j $Nproc $repairSh \
 in1=${READSDIR}/{}R1$EXTENSION \
 in2=${READSDIR}/{}R2$EXTENSION \
 out1=repairFilt/{}repr_R1.fq.gz out2=repairFilt/{}repr_R2.fq.gz outs=repairFilt/{}repr_orph.fq repair

# may change before proceed
##  may need a customized script to filter barcode reads according to filtered repaired reads

echo "TIME: END repair = `date +"%Y-%m-%d %T"`"

# Part 2: QC with fastqscreen after repair.
mkdir -p fqscreenQCafFtRp
module load bowtie2

echo "TIME: START  QC = `date +"%Y-%m-%d %T"`"

for f in repairFilt/*repr_R1.fq.gz
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
mkdir -p ${ResultDir}/q2_viz
#############################################################
#### SINGULARITYENV (space & runtime issues) ################
XDG_RUNTIME_DIR=${SCRATCH}/tmp_XDG
SINGTEMPDIR=${SCRATCH}/tmp_sing
export SINGULARITYENV_XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR}
export SINGULARITYENV_TMPDIR=${SINGTEMPDIR}
#############################################################
echo "VERSION: SINGULARITY: `singularity --version`"

# Part 3: Import the multiplexed sequences. Demultiplexing and Trimming Adapters from Reads with q2-cutadapt
# Should allow for slightly changing naming schemes by core, but always assumes "_R1" somewhere in filename indicates read1 though.
mkdir  -p importPe
cp ./repairFilt/*_R1.fq.gz forward.fastq.gz
gunzip forward.fastq.gz
cp ./repairFilt/*_R2.fq.gz reverse.fastq.gz
gunzip reverse.fastq.gz
cp /uufs/chpc.utah.edu/common/home/u0762203/u0762203/project/hammer/reads/*I1.fastq.gz barcodes.fastq.gz
gunzip barcodes.fastq.gz
# The order of the records in the fastq.gz files defines the association between a sequence read and its barcode read (i.e. the first barcode read corresponds to the first sequence read, the second barcode to the second read, and so on.). need filter the barcode read according to the filtered and repaired first sequence read

filterBarcodeReadsAccording2fwdReads.pl forward.fastq barcodes.fastq reverse.fastq 
rm forward.fastq
rm barcodes.fastq
rm reverse.fastq
mv forward.rmTag.fastq forward.fastq
mv reverse.rmTag.fastq reverse.fastq
mv barcodes.filtered.fastq barcodes.fastq
gzip barcodes.fastq
gzip forward.fastq
gzip reverse.fastq
mv barcodes.fastq.gz importPe/
mv forward.fastq.gz importPe/
mv reverse.fastq.gz importPe/

echo "TIME: START import = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime tools import --type EMPPairedEndSequences \
	--input-path importPe \
	--output-path ${ResultDir}/multiplexed-seqs.qza
echo "TIME: END import = `date +"%Y-%m-%d %T"`"

echo "TIME: START demultiplex = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime demux emp-paired \
	 --i-seqs ${ResultDir}/multiplexed-seqs.qza \
	 --m-barcodes-file ${MAPBYGNOMEXID} \
	 --m-barcodes-column BarcodeSequence \
	 --o-per-sample-sequences ${ResultDir}/demultiplexed-seqs.qza \
	 --o-error-correction-details ${ResultDir}/errorCorrectionDetails.qza \
	 --p-rev-comp-mapping-barcodes \
         --verbose

echo "TIME: END demultiplex = `date +"%Y-%m-%d %T"`"

#trim adapters from demultiplexed reads 

echo "TIME: START trim = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime cutadapt trim-paired \
	--i-demultiplexed-sequences ${ResultDir}/demultiplexed-seqs.qza \
	--p-front-f ${FPrimerSeqToTrim} \
        --p-front-r ${RPrimerSeqToTrim} \
	--o-trimmed-sequences ${ResultDir}/demux_trim.qza \
	--p-cores $Nproc

# At this stage plots should be inspected to infer the --p-trunc-len-f, --p-trunc-len-r, --p-trim-left-f and --p-trim-left-r for dada2. It should hold whenever using the same primer set though with good quality seq runs.
singularity exec ${ContainerImage} qiime demux summarize \
  --i-data ${ResultDir}/demux_trim.qza \
  --o-visualization ${ResultDir}/q2_viz/demux_trim.qzv  
echo "TIME: END trim = `date +"%Y-%m-%d %T"`"


#dada2 QC parameters:Position at which forward/reverse read sequences should be truncated due to decrease in quality, need to be adjusted according to inspection of Interactive Quality Plot tab in the demux_trim.qzv (shoot for middle box 20 - 30--> ~28? maybe)

