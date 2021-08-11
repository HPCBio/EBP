#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J EB20_pilon_TellSeq_reads #jobname
#SBATCH --mem=250g #memory requested
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -A ebp #account
#SBATCH -p normal #queue
#SBATCH -w compute-7-4
#SBATCH --mail-user=kosterbu@illinois.edu #email address
#SBATCH --mail-type=ALL #get all job notifications
#SBATCH -D /home/groups/earthbiogenome/results/20210407_EB20_ArrowPolished_purged_scaffolding/

# ----------------Load Modules--------------------

module load BWA/0.7.17-IGB-gcc-8.2.0
module load SAMtools/1.11-IGB-gcc-8.2.0 

# ----------------Commands------------------------

SCRATCH=/scratch/kosterbu/EB20/

#mkdir -p $SCRATCH

#cp /home/groups/earthbiogenome/results/20210407_EB20_ArrowPolished_purged_scaffolding/EB20_redbean_g2.7_80X_arrow_polished_purged_purged.fa $SCRATCH

#cp /home/groups/earthbiogenome/data/TellSeq_not_interleaved_reads/EB20_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz $SCRATCH
#cp /home/groups/earthbiogenome/data/TellSeq_not_interleaved_reads/EB20_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz $SCRATCH

cd $SCRATCH

export TMPDIR=/scratch/kosterbu/EB20/

# index the assembly for bwa
#bwa index EB20_redbean_g2.7_80X_arrow_polished_purged_purged.fa

# run bwa, use 20 threads for aligning and sorting and set memory for samtools at 3G per thread
#bwa mem -t 12 EB20_redbean_g2.7_80X_arrow_polished_purged_purged.fa \
# EB20_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
# EB20_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz | samtools sort -@12 -m 3G -o EB20_TellSeq_sort.bam -

#samtools index EB20_TellSeq_sort.bam

#module purge

# ----------------Load Modules--------------------

module load pilon/1.23-Java-1.8.0_152 

# ----------------Commands------------------------

# run pilon


java -Xmx250G -jar /home/apps/software/pilon/1.23-Java-1.8.0_152/pilon-1.23.jar --genome EB20_redbean_g2.7_80X_arrow_polished_purged_purged.fa \
    --bam EB20_TellSeq_sort.bam \
    --changes --diploid --threads 12 \
    --output EB20_Arrow_purged_TellSeq_polished
    
cp EB20* /home/groups/earthbiogenome/results/20210407_EB20_ArrowPolished_purged_scaffolding/
