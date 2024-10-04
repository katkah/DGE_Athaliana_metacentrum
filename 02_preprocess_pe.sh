#!/bin/bash
#PBS -l select=1:ncpus=10:mem=45gb:scratch_local=100gb
#PBS -l walltime=30:00:00
#PBS -N 02_preprocess


### Variables
INPUT_DIR="/storage/storage_name/home/user/RNAseq/raw_fastq"
OUTPUT_DIR="/storage/storage_name/home/user/RNAseq/preprocessed_data"
OUTPUT_DIR_FASTQC="/storage/storage_name/home/user/RNAseq/preprocessed_qc"
ADAPTERS="/storage/storage_name/home/user/adapters_merge.fa"

APPENDIX1="R1.fastq.gz"
APPENDIX2="R2.fastq.gz"
THREADS=10
# Number of threads to use


### copy inputs
cd $SCRATCH
mkdir -p $SCRATCH/raw_data
cp $INPUT_DIR/* $SCRATCH/raw_data
cp $ADAPTERS $SCRATCH/ 
ADAPTERS=$SCRATCH/$(basename $ADAPTERS)

mkdir -p $SCRATCH/preprocessed_data
mkdir -p $SCRATCH/preprocessed_qc

### Trimmomatic fastqc: Preprocess and qc
source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load fastqc
module load trimmomatic

for i in $SCRATCH/raw_data/*$APPENDIX1
do
	READ_FOR=$(basename $i)
	READ_REV=${i%$APPENDIX1}$APPENDIX2
	READ_REV=$(basename $READ_REV)

	echo "Now I am processing PE reads $SCRATCH/raw_data/$READ_FOR $SCRATCH/raw_data/$READ_REV - Trimming"

	trimmomatic PE -threads ${THREADS} \
        $SCRATCH/raw_data/$READ_FOR \
        $SCRATCH/raw_data/$READ_REV \
        $SCRATCH/preprocessed_data/${READ_FOR%$APPENDIX1}R1_trim.fastq.gz \
        $SCRATCH/preprocessed_data/${READ_FOR%$APPENDIX1}R1un_trim.fastq.gz \
        $SCRATCH/preprocessed_data/${READ_REV%$APPENDIX2}R2_trim.fastq.gz \
        $SCRATCH/preprocessed_data/${READ_REV%$APPENDIX2}R2un_trim.fastq.gz \
        SLIDINGWINDOW:4:15 \
        HEADCROP:9 \
        ILLUMINACLIP:$ADAPTERS:2:30:10\
		MINLEN:35 &> $SCRATCH/preprocessed_data/${READ_REV%$APPENDIX2}trim.log
	echo "Done processing PE reads $READ_FOR $READ_REV - Trimming"

done

fastqc --outdir $SCRATCH/preprocessed_qc --format fastq --threads $THREADS $SCRATCH/preprocessed_data/*_trim.fastq.gz


### multiqc
module add mambaforge
mamba activate multiqc_v1.12_py3.7
multiqc -o $SCRATCH/preprocessed_qc $SCRATCH/preprocessed_qc


### copy results
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR_FASTQC
cp $SCRATCH/preprocessed_data/* $OUTPUT_DIR
cp $SCRATCH/preprocessed_qc/* $OUTPUT_DIR_FASTQC

clean_scratch
