#!/bin/bash -e
#PBS -l select=1:ncpus=6:mem=35gb:scratch_local=150gb
#PBS -l walltime=38:00:00
#PBS -N 04_alignment_star_se

# script uses star to align SE (single end) sequencing data

### Variables
# set variables to data location e.g.:
# set the path to star_index created by 03_prepare_index_star.sh
GENOME_DIR="/storage/storage_name/home/user/RNAseq/STAR_index"
# set the path to gtf file
GTF="/storage/storage_name/home/user/RNAseq/Arabidopsis_thaliana.TAIR10.58.gtf"
INPUT_DIR="/storage/storage_name/home/user/RNAseq/preprocessed_data"
OUTPUT_DIR="/storage/storage_name/home/user/RNAseq/alignments_star_se"

# the suffix of preprocessed data to be aligned by star, usually it is _trim.fastq.gz but can be different
APPENDIX="_trim.fastq.gz"

THREADS=6

############################################################################################
### Copy inputs
mkdir $SCRATCH/processed_data
cp $INPUT_DIR/* $SCRATCH/processed_data/
cp -r $GENOME_DIR/ $SCRATCH/
cp $GTF $SCRATCH/
GTF=$(basename $GTF)

####################################################################################################
###STAR
source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load star

mkdir $SCRATCH/alignments
for i in $SCRATCH/processed_data/*$APPENDIX
do
    READ_FOR=$(basename $i)
	OUT_PREFIX=${READ_FOR%$APPENDIX}
	echo "Now I am processing SE reads $SCRATCH/processed_data/$READ_FOR"
	STAR --runThreadN $THREADS\
		--genomeDir $SCRATCH/STAR_index\
		--readFilesIn $SCRATCH/processed_data/$READ_FOR\
        --readFilesCommand zcat\
		--outSAMtype BAM SortedByCoordinate\
        --sjdbGTFfile $SCRATCH/$GTF\
		--outFileNamePrefix $SCRATCH/alignments/$OUT_PREFIX\
		--outFilterMultimapNmax 20\
		--outFilterMismatchNoverReadLmax 1.00\
		--outFilterMismatchNmax 999\
        --quantTranscriptomeBan IndelSoftclipSingleend\
        --quantMode TranscriptomeSAM GeneCounts
	echo "Done processing SE reads $SCRATCH/processed_data/$READ_FOR"
done



###samtools
cd $SCRATCH/alignments

source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load samtools

for i in $SCRATCH/alignments/*sortedByCoord.out.bam
do
    samtools index $i
done


### Copy results
mkdir -p $OUTPUT_DIR

cp -r $SCRATCH/alignments/* $OUTPUT_DIR/

clean_scratch
