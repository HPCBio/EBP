#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J EB15_Asset_TellSeq_reads #jobname
#SBATCH --mem=100g #memory requested
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -p normal #queue
#SBATCH -A ebp
#SBATCH --mail-user=kosterbu@illinois.edu #email address
#SBATCH --mail-type=ALL #get all job notifications
#SBATCH -D /home/groups/earthbiogenome/results/20210420_EB15_Asset_Merqury/

# ----------------Load Modules--------------------

module load asset/1.0.3-IGB-gcc-8.2.0
module load minimap/2.17-IGB-gcc-8.2.0
module load BWA/0.7.17-IGB-gcc-8.2.0
module load SAMtools/1.11-IGB-gcc-8.2.0

# ----------------Commands------------------------

# TellSeq reads: align TellSeq reads to assembly fasta file with bwa; output is bam file

##bwa index EB15_Sealer_scaffold.fa

#ln -s /home/groups/earthbiogenome/data/TellSeq_like_Supernova_reads/EB15_nobc_4tenx_S1_L001_R1_001.fastq.gz .
#ln -s /home/groups/earthbiogenome/data/TellSeq_like_Supernova_reads/EB15_nobc_4tenx_S1_L001_R2_001.fastq.gz .

# Usage: 10x [options...] <in1.fq/in1.fq.gz> <in2.fq/in2.fq.gz>

# 10x -l 18 -m 19 -c -p EB15_TellSeq_10x_format_trimmed EB15_nobc_4tenx_S1_L001_R1_001.fastq.gz EB15_nobc_4tenx_S1_L001_R2_001.fastq.gz

bwa mem -t 12 EB15_Sealer_scaffold.fa \
 EB15_TellSeq_10x_format_trimmed_1.fq.gz \
 EB15_TellSeq_10x_format_trimmed_2.fq.gz | samtools view -b - > EB15_TellSeq_10x_format_trimmed.bam

# Convert bam to bed format
# Usage: ast_10x [options] <GAP_BED> <BAM_FILEs> ...

ast_10x ./EB15_Sealed_gaps.bed EB15_TellSeq_10x_format_trimmed.bam > ./EB15_Sealer_scaff_TellSeq_10x_format_trimmed.bed 2>EB15_ast_10x_TellSeq_10xFormat.log
