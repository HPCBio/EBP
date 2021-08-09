#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J diamond_for_blobtools #jobname
#SBATCH --mem=300g
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -p hpcbio #queue
#SBATCH --mail-user=your_email@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -D /home/groups/path/to/working/directory/

# ----------------Load Modules--------------------

module load DIAMOND/2.0.6-IGB-gcc-8.2.0

# ----------------Commands------------------------

# run this on /scratch

SCRATCH=/scratch/path/to/run/diamond/

mkdir -p $SCRATCH

# stage assembly to /scratch

cp supernova_fasta_file.fa $SCRATCH

# stage UniRef100 diamond db to /scratch

cp uniref100.dmnd.gz $SCRATCH 

cd $SCRATCH

# Run diamond in "blastx" mode:

gunzip uniref100.dmnd.gz

diamond blastx -p 12 -t $SCRATCH -d uniref100.dmnd -q supernova_fasta_file.fa --evalue 1e-25 --sensitive --outfmt 100 -o matches.daa 

diamond view -a matches.daa -f 6 -o matches.m8 
#-a is view input file, -f 6 is tab-delimited format, -o is view output file

# Copy matches* files back to working directory

cp matches* /home/groups/path/to/working/directory/

#clean up /scratch as needed