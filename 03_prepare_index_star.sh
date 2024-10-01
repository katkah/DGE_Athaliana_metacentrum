#!/bin/bash
#PBS -l select=1:ncpus=6:mem=25gb:scratch_local=30gb
#PBS -l walltime=4:00:00
#PBS -N 03_prepare_index_star
#

# script creates star index, run only once for a given genome and a read length 

### Variables
# set your variables, e.g.
GENOME="/storage/storage_name/home/user/RNAseq/TAIR10_release58.fa"
GTF="/storage/storage_name/home/user/RNAseq/Arabidopsis_thaliana.TAIR10.58.gtf"
INDEX_OUTPUT_DIR="/storage/storage_name/home/user/RNAseq/STAR_index"
#read length -1 
RD_LENGTH=99 
THREADS=6 # Number of threads to use



### Copy inputs
cp $GTF $SCRATCH/
cp $GENOME $SCRATCH/
cd $SCRATCH/
GENOME=$(basename $GENOME)
GTF=$(basename $GTF)


### STAR create index

source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load star

mkdir -p $SCRATCH/STAR_index

STAR --runMode genomeGenerate --runThreadN $THREADS --genomeDir $SCRATCH/STAR_index --genomeFastaFiles $GENOME \
--sjdbGTFfile $GTF --sjdbOverhang $RD_LENGTH --genomeSAindexNbases 12 
# --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ


### Copy results
mkdir -p $INDEX_OUTPUT_DIR

cp -r $SCRATCH/STAR_index/* $INDEX_OUTPUT_DIR/

clean_scratch
