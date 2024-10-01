#!/bin/bash
#PBS -l select=1:ncpus=6:mem=25gb:scratch_local=30gb
#PBS -l walltime=4:00:00
#PBS -N 03_prepare_index_rsem
#


# script creates rsem index, run only once for a given genome 

### Variables
GENOME="/storage/storage_name/home/user/RNAseq/TAIR10_release58.fa"
GTF="/storage/storage_name/home/user/RNAseq/Arabidopsis_thaliana.TAIR10.58.gtf"
INDEX_OUTPUT_DIR="/storage/storage_name/home/user/RNAseq/"
THREADS=6
NAME="Arabidopsis_thaliana.TAIR10.58"

### Copy inputs
cp $GTF $SCRATCH/
cp $GENOME $SCRATCH/
cd $SCRATCH/
GENOME=$(basename $GENOME)
GTF=$(basename $GTF)



### Indexing
source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load rsem
mkdir $SCRATCH/RSEM_index
rsem-prepare-reference --gtf $GTF $GENOME $SCRATCH/RSEM_index/$NAME


### Copy results
mkdir -p $INDEX_OUTPUT_DIR
cp -r $SCRATCH/RSEM_index $INDEX_OUTPUT_DIR/

clean_scratch
