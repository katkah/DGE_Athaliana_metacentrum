#!/bin/bash
#PBS -l select=1:ncpus=6:mem=35gb:scratch_local=150gb
#PBS -l walltime=8:00:00
#PBS -N 01_initialQuality

# script uses fastqc to check initial quality


### Variables
# set variables to data location e.g.:
INPUT_DIR="/storage/storage_name/home/user/RNAseq/raw_fastq"
OUTPUT_DIR="/storage/storage_name/home/user/RNAseq/raw_qc"
# suffix of files that are to be processed, usually fastq.gz, but can be different
APPENDIX="fastq.gz"
# number of threads to be used
THREADS=6


### Copy inputs to scratch
cp $INPUT_DIR/*$APPENDIX $SCRATCH/
cd $SCRATCH/


### FastQC QC
mkdir -p $SCRATCH/fastqc

source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load fastqc
fastqc -t $THREADS -f fastq -o $SCRATCH/fastqc *$APPENDIX

module load conda-modules
conda activate multiqc_v1.12_py3.7
multiqc -o $SCRATCH/fastqc $SCRATCH/fastqc
conda deactivate

cd $SCRATCH/fastqc

### Copy results
mkdir -p $OUTPUT_DIR

cp -r $SCRATCH/fastqc $OUTPUT_DIR/

clean_scratch

