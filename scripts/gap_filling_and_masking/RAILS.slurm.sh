#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J RAILS #jobname
#SBATCH --mem=50g
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -p normal #queue
#SBATCH -A ebp
#SBATCH --mail-user=your_email@illinois.edu #email address
#SBATCH --mail-type=ALL
#SBATCH -D /home/groups/earthbiogenome/results/path/to/your/working/directory/

# ----------------Load Modules--------------------

module load seqkit/0.12.1

# ----------------Commands------------------------

# modify your email in parameters above
# modify your working directory in parameters above

# change "prefix" to correct file names

# Convert wrapped SALSA scaffold FASTA file to single line FASTA file

seqkit seq -w 0 scaffolds_FINAL.fasta > prefix_salsa_scaffolds_FINAL_noWrap.fasta

module purge

# ----------------Load Modules--------------------

module load RAILS/1.5.1-IGB-gcc-8.2.0
module load BWA/0.7.17-IGB-gcc-8.2.0 
module load picard/2.10.1-Java-1.8.0_152

# ----------------Commands------------------------

gunzip /home/groups/earthbiogenome/data/PacBio_HiFi_reads/prefix.fastq.gz

ln -s /home/groups/earthbiogenome/data/PacBio_HiFi_reads/prefix.fastq .

# Convert fastq to fasta (script outputs a fasta file with ".fasta" appended to end of fastq filename)

perl /home/groups/earthbiogenome/src/Fastq2Fasta.pl prefix.fastq

# Run RAILS
# "Usage: $(basename $0) <FASTA assembly .fa> <FASTA long sequences .fa> <anchoring sequence length eg. 250> <min sequence identity 0.95> <max. softclip eg. 250bp> <min. number of read support eg. 2> <long read type eg.: ont, pacbio, nil> <path to samtools>"

sh /home/groups/earthbiogenome/src/runRAILSminimapSTREAM.sh prefix_salsa_scaffolds_FINAL_noWrap.fasta prefix.fastq.fasta 1000 0.95 250 3 pacbio /home/apps/software/SAMtools/1.11-IGB-gcc-8.2.0/bin/