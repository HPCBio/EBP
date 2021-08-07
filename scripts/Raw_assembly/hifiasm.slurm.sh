#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=50g # memory requested
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -A ebp #account
#SBATCH -p normal #queue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@illinois.edu
#SBATCH -J EBxx_hifiasm #jobname
#SBATCH -D /path/to/your/working/directory/

### edit lines 11-13 above and variables below

### next lines are included to help troubleshoot the run
### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -x
set -euo pipefail
IFS=$'\n\t'
echo `hostname`

echo `date`
echo ################################ Declare variables and sanity check
### Declare variables
### These variables need editing
### HiFi is the file of HiFi reads
### OUTPUTDIR is the folder where results will be written to
 
PREFIX=EBxx
OUTPUTDIR=/path/to/results/<date>_${PREFIX}_hifiasm
HiFi=/path/to/HiFi_consensus.fastq

if [ ! -s "$HiFi" ]
then
	die "$HiFi file not found"
fi
if [ ! -d "$OUTPUTDIR" ]
then
	mkdir -p $OUTPUTDIR
fi

echo `date`
echo ################################ Run hifiasm

module load hifiasm/0.13-IGB-gcc-8.2.0

hifiasm -o $PREFIX}_hifiasm -t $SLURM_NPROCS $HiFi

echo `date`
echo ################################ Run gfa2fa 

module purge

module load gfatools/0.4-IGB-gcc-4.9.4

gfatools gfa2fa $PREFIX}_hifiasm.p_ctg.gfa > $PREFIX}_hifiasm.p_ctg.fasta


