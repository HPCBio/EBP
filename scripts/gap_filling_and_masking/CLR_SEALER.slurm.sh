#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J EB20_SEALER #jobname
#SBATCH --mem=50g
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -p normal #queue
#SBATCH -A ebp
#SBATCH --mail-user=kosterbu@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -D /home/groups/earthbiogenome/results/20210428_EB20_RAILS_SEALER/

# ----------------Load Modules--------------------

module load abyss/2.2.5-IGB-gcc-8.2.0

# ----------------Commands------------------------

#usage: options -b <Bloom filter size> -k <size of kmer> -k <size of kmer> -j <number of threads/cpus> --mask <mask new and changed bases as lower case> -o <output prefix> -S <path to fasta to gap-fill> <paths to short_reads>

abyss-sealer -b20G -k64 -k96 -j 12 --mask -o EB20_Salsa_Sealer \
-S ./EB20_salsa_scaffolds_FINAL_noWrap.fasta \
/home/groups/earthbiogenome/data/HiC_reads/EB20_ACTTGA_L003_R1_001.fastq.gz \
/home/groups/earthbiogenome/data/HiC_reads/EB20_ACTTGA_L003_R2_001.fastq.gz \
/home/groups/earthbiogenome/data/TellSeq_not_interleaved_reads/EB20_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
/home/groups/earthbiogenome/data/TellSeq_not_interleaved_reads/EB20_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz