# DGE_Athaliana_metacentrum

### Prerequisites for Analysis

Before starting the analysis, ensure you are familiar with running scripts on MetaCentrum and have a basic understanding of bash scripting. Helpful resources include:

- [MetaCentrum Beginner’s Guide](https://wiki.metacentrum.cz/wiki/Pruvodce_pro_zacatecniky)
- [Introduction to Bash Scripting](https://www.geeksforgeeks.org/bash-scripting-introduction-to-bash-and-bash-scripting/)

### Initial Quality Check

**Script: `01_initialQuality.sh`**

1. Set the paths and variables within the script.
2. Estimate the computational resources required based on your dataset size, and define these at the top of the script.
3. Create a directory in your home space (e.g., `RNAseq/scripts/`).
4. Copy the scripts to MetaCentrum and place them in the directory.
5. Run the `dos2unix` command on the script to convert any Windows-specific characters to Unix format:
   ```
   dos2unix 01_initialQuality.sh
   ```
6. Submit the script to the queue using the `qsub` command:
   ```
   qsub 01_initialQuality.sh
   ```
7. Verify that the script is running with the following command (replace `username` with your actual username):
   ```
   qstat -u username
   ```

### Preprocessing Data

**Single-End Data Script: `02_preprocess_se.sh`**

This script uses Trimmomatic to preprocess single-end data. Follow these steps:

1. Set the appropriate paths and variables.
2. Review the results of `01_initialQuality.sh` and adjust Trimmomatic parameters accordingly.
3. For guidance on interpreting FastQC results, refer to this [FastQC assessment guide](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html).
4. Find the Trimmomatic manual [here](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

#### Important Considerations:
- **Adapter Content:** If the adapter content is high, use a FASTA file to specify adapters for removal. Commonly used file: `adapters_merge.fa`. Note that using many adapters may slow down the process (2-3 days for large datasets). Options include:
  - Waiting for the process to complete.
  - Using tools like TrimGalore or cutadapt.
  - Reducing the number of adapters in the FASTA file.
  - Using built-in adapters provided by Trimmomatic.
- **GC Content:** Check for smooth, unimodal GC content, and look for potential contamination.
- **Sequence Quality:** Assess the sequence quality and consider trimming the reads if needed.

The preprocessed data will be in the `preprocessed_data` and `preprocessed_qc` directories. Check the `preprocessed_qc` folder to confirm the data is ready for mapping.

**Paired-End Data Script: `02_preprocess_pe.sh`**

This script is for paired-end data, also using Trimmomatic. The steps are similar to the single-end data script.

### Mapping Data to a Reference Genome

Before alignment, download the Arabidopsis genome (or another reference genome) from [Araport](https://www.araport.org) or [EnsemblPlants](https://plants.ensembl.org). Below are the steps to prepare the TAIR version 58 genome, though newer versions are available:

#### Genome Preparation Steps:

1. Download the genome files:
   ```bash
   wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.Pt.fa.gz
   wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.Mt.fa.gz
   wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa.gz
   ...
   ```
2. Concatenate the chromosomes into one file:
   ```bash
   cat Arabidopsis_thaliana.TAIR10.dna.chromosome.*.fa.gz > TAIR10_release58.fa.gz
   ```
3. Download the genome annotation in GFF3 format from EnsemblPlants.
4. Convert the GFF3 file to GTF format using `agat`:
   ```bash
   agat_convert_sp_gff2gtf.pl --gtf_version 3 --gff Arabidopsis_thaliana.TAIR10.58.gff3 -o Arabidopsis_thaliana.TAIR10.58.gtf
   ```

### Mapping Reads with STAR

**Script: `03_prepare_index_star.sh`**

- This script creates the STAR index, which is required for alignment. You can reuse the index for datasets with the same genome and read length. Refer to the [STAR manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for detailed instructions.

**Single-End Alignment Script: `04_alignment_star_se.sh`**

- This script aligns single-end reads to both the transcriptome and genome using STAR. If `--quantMode TranscriptomeSAM GeneCounts` is set, it will also count reads, helping to determine strandedness for the next step (RSEM).

**Paired-End Alignment Script: `04_alignment_star_pe.sh`**

- This script is similar to the single-end alignment script but is tailored for paired-end reads.

### Counting Reads with RSEM

**Single-End Read Counting Script: `05_rsem_se.sh`**

- This script counts reads per isoform and gene from single-end alignments. The results will have a suffix `.Aligned.toTranscriptome.out.rsem.genes`, with the TPM value in the 6th column.

**Paired-End Read Counting Script: `05_rsem_pe.sh`**

- This script performs the same function for paired-end data.
