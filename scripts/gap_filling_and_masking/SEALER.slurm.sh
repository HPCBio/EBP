#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J SEALER #jobname
#SBATCH --mem=50g
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -p normal #queue
#SBATCH -A ebp
#SBATCH --mail-user=youremail@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -D /path/to/working/directory/

# ----------------Load Modules--------------------

module load abyss/2.2.5-IGB-gcc-8.2.0

# ----------------Commands------------------------

#usage: options -b <Bloom filter size> -k <size of kmer> -k <size of kmer> -j <number of threads/cpus> --mask <mask new and changed bases as lower case> -o <output prefix> -S <path to fasta to gap-fill> <paths to short_reads>

abyss-sealer -b20G -k64 -k96 -j 12 --mask -o EBXX_Sealer \
-S ./EBXX_rails.scaffolds.fa \
/home/groups/earthbiogenome/data/HiC_reads/EBXX.R1.fastq.gz \
/home/groups/earthbiogenome/data/HiC_reads/EBXX.R2.fastq.gz \
/home/groups/earthbiogenome/data/TellSeq_not_interleaved_reads/EBXX.R1.fastq.gz \
/home/groups/earthbiogenome/data/TellSeq_not_interleaved_reads/EBXX.R2.fastq.gz