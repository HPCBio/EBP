#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J bwa_for_blobtools #jobname
#SBATCH --mem=100g
#SBATCH -N 1 #nodes
#SBATCH -n 24 #cpus-per-task
#SBATCH -p hpcbio #queue
#SBATCH --mail-user=your_mail@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -D /home/groups/path/to/working/directory/

# ----------------Load Modules--------------------

module load BWA/0.7.17-IGB-gcc-8.2.0

# ----------------Commands------------------------

# map short reads to genome assembly
# bwa mem will have to be adjusted if short reads are interleaved, use -p
# adapter/barcodes should probably be removed
  

bwa mem -t 24 supernova_fasta_file.fa /path/to/read1.fastq.gz /path/to/read2.fastq.gz  > prefix.sam

# ----------------Load Modules--------------------

module purge

module load SAMtools/1.11-IGB-gcc-8.2.0

# ----------------Commands------------------------

samtools view -@ 12 -bo prefix.bam prefix.sam

samtools sort -@ 12 -o prefix_sorted.bam prefix.bam 

samtools index -@ 12 prefix_sorted.bam