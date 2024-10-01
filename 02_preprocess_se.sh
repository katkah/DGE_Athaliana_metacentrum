#!/bin/bash 
#PBS -l select=1:ncpus=12:mem=45gb:scratch_local=150gb
#PBS -l walltime=30:00:00
#PBS -N 02_preprocess

# script uses trimmomatic to preprocess SE (single end) sequencing data

### Variables
# set variables to data location e.g.:
INPUT_DIR="/storage/storage_name/home/user/RNAseq/raw_fastq"
OUTPUT_DIR="/storage/storage_name/home/user/RNAseq/preprocessed_data"
OUTPUT_DIR_FASTQC="/storage/storage_name/home/user/RNAseq/preprocessed_qc"
ADAPTERS="/storage/storage_name/home/user/adapters_merge.fa"
# suffix of files that are to be processed, usually fastq.gz, but can be different
APPENDIX=".fastq.gz"
# Number of threads to use
THREADS=12


### copy inputs to scratch
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

for i in $SCRATCH/raw_data/*$APPENDIX
do
	READ_FOR=$(basename $i)

	echo "Now I am processing SE reads $SCRATCH/raw_data/$READ_FOR - Trimming"
    #set trimmomatic parameters based on the trimmomatic manual and your data requirements
	trimmomatic SE -threads $THREADS $SCRATCH/raw_data/$READ_FOR \
		$SCRATCH/preprocessed_data/${READ_FOR%$APPENDIX}_trim.fastq.gz \
		SLIDINGWINDOW:4:20 \
        HEADCROP:12 \
		MINLEN:35 &> $SCRATCH/preprocessed_data/${READ_FOR%$APPENDIX}_trim.log
	echo "Done processing SE reads $READ_FOR - Trimming"

done

fastqc --outdir $SCRATCH/preprocessed_qc --format fastq --threads $THREADS $SCRATCH/preprocessed_data/*_trim.fastq.gz

### multiqc
module load conda-modules
conda activate multiqc_v1.12_py3.7
multiqc -o $SCRATCH/preprocessed_qc $SCRATCH/preprocessed_qc
conda deactivate

### copy results
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR_FASTQC
cp $SCRATCH/preprocessed_data/* $OUTPUT_DIR
cp $SCRATCH/preprocessed_qc/* $OUTPUT_DIR_FASTQC

clean_scratch
