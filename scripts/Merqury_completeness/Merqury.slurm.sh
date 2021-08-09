#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J merqury #jobname
#SBATCH --mem=150g
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -p normal #queue
#SBATCH -A ebp
#SBATCH --mail-user=your_name@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -D /home/groups/earthbiogenome/results/path/to/directory/

# ----------------Load Modules--------------------

module load merqury/1.3-IGB-gcc-8.2.0

# ----------------Commands------------------------

# Count kmers in TellSeq reads:

# adjust k value depending on estimated genome size:
#kmer length" field with 19, 20 or 21 depending on if the size of your working genome is < 600 Mbp, 600 Mbp to ~2.2 Gbp or >= 2.3 Gbp

# Adjust prefix/filenames as needed.

meryl k=XX count output EBXX_TellSeq_read1.meryl /home/groups/earthbiogenome/data/TellSeq_not_interleaved_reads/EBXX_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz

meryl k=XX count output EBXX_TellSeq_read2.meryl /home/groups/earthbiogenome/data/TellSeq_not_interleaved_reads/EBXX_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz

# Output from count is 2 new directories; one for each read file.

# Join the kmer output into one dataset 
# Usage: meryl union-sum output  <joined meryl output directory> <prefix of separate kmer directories>

meryl union-sum output EBXX_TellSeq.meryl EBXX_TellSeq_read*.meryl

# If running Merqury in a new results directory, symlink Sealer scaffold file into new working directory:
# Otherwise, comment out the following line:
ln -s /path/to/Sealer.fa/ /path/to/merqury/output/

# run merqury on joined meryl dataset and Sealer scaffold fasta file
# Usage: merqury.sh <Joined Meryl output directory> <EBXX Sealer scaffold fasta file> <prefix for Merqury output files>

merqury.sh EBXX_TellSeq.meryl EBXX_Sealer_scaffold.fa EBXX_merqury