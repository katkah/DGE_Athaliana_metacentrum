#!/bin/bash -e
#PBS -l select=1:ncpus=6:mem=36gb:scratch_local=200gb
#PBS -l walltime=24:00:00
#PBS -N 05_gene_counts


# Gene/isoform counting - RNA-Seq Single end
#
# Input is from STAR mapping to the genome and the transcriptome
# You should know the strandness of the experiment (either by kit name or check STAR gene counts result)
# How to interpret STAR gene counts to figure out the strandness is discussed here https://www.biostars.org/p/218995/
#
#for paired-end use:
#    --paired-end parameter with RSEM



### Variables
#set variables e.g.
INPUT_DIR="/storage/storage_name/home/user/RNAseq/alignments_star_se/transcriptome"
OUTPUT_DIR="/storage/storage_name/home/user/RNAseq/gene_counts_rsem"
OUTPUT_DIR_QC="/storage/storage_name/home/user/RNAseq/qc_rsem"

STRAND="fwd" # [fwd, rev, none] strandedness of reads 

GENOME="/storage/storage_name/home/user/RNAseq/TAIR10_release58.fa"
GTF="/storage/storage_name/home/user/RNAseq/Arabidopsis_thaliana.TAIR10.58.gtf"
RSEM_GENOME_INDEX="/storage/storage_name/home/user/RNAseq/RSEM_index"
RSEM_GENOME_INDEX_NAME="Arabidopsis_thaliana.TAIR10.58"

THREADS=6
RSEM_RANDOM=123456



### Load tools
source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load rsem
RSEM=rsem-calculate-expression

module load samtools
SAMTOOLS=samtools

module load conda-modules
conda activate multiqc_v1.12_py3.7
MULTIQC=multiqc


### Copy inputs
cp $GTF $SCRATCH/
cp $GENOME $SCRATCH/
cp $RSEM_GENOME_INDEX/* $SCRATCH/

cp $INPUT_DIR/*out.bam $SCRATCH/ 

cd $SCRATCH/

### Genome and annotation preparation
GENOME=$(basename $GENOME)
GTF=$(basename $GTF)
RSEM_GENOME_INDEX=$RSEM_GENOME_INDEX_NAME


### Gene/isoform counting

# Count reads using RSEM
mkdir -p $SCRATCH/gene_counts/rsem

echo "Using $RSEM_RANDOM as --seed"

for i in *out.bam
do
	echo "Started RSEM counting $i"
	$RSEM -p $THREADS \
    --alignments \
    --estimate-rspd \
    --calc-ci \
    --seed $RSEM_RANDOM \
    --no-bam-output \
    --ci-memory 32000 \
    --strandedness $STRAND \
    $i $RSEM_GENOME_INDEX $SCRATCH/gene_counts/rsem/${i%.*}.rsem
	echo "Done RSEM counting $i"
done

mkdir $SCRATCH/RSEM_gc

$MULTIQC -o $SCRATCH/RSEM_gc $SCRATCH/gene_counts/rsem/ &

### Copying results
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR_QC

cp -r $SCRATCH/gene_counts/ $OUTPUT_DIR/
cp -r $SCRATCH/RSEM_gc/ $OUTPUT_DIR_QC/

clean_scratch
