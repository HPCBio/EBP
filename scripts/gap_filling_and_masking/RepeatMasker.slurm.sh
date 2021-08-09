#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=10g #memory requested
#SBATCH -N 1 #nodes
#SBATCH -n 4 #cpus-per-task
#SBATCH -A ebp #account
#SBATCH -p normal #queue
#SBATCH --mail-type=ALL #get all job notifications
#SBATCH -J EBxx_RepeatMasker #jobname
#SBATCH --mail-user=youremail@illinois.edu #email address
#SBATCH -D /path/to/working/directory/

# ----------------Declare variables--------------------
# GENOME should be the scaffolds.fa file after the Sealer step


GENOME=/full/path/to/input_scaffolds.fa
OUTPUTDIR=$(basename $GENOME)
SPECIES=insecta

# ----------------Load Modules--------------------

module load RepeatMasker/4.0.7-IGB-gcc-4.9.4-Perl-5.24.1

# ----------------Commands------------------------

#usage: RepeatMasker -species <taxon> <Sealer_scaffold.fasta file>

cd $OUTPUTDIR

RepeatMasker -species $SPECIES $GENOME