#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=100g # memory requested
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -A ebp #account
#SBATCH -p normal #queue
#SBATCH --mail-type=ALL
#SBATCH -J EBxx_jellyfish #jobname
#SBATCH --mail-user=youremail@illinois.edu
#SBATCH -D /path/to/working/directory/


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
### R1 and R2 are the files of TellSeq reads
### OUTPUTDIR is the folder where results will be written to
 
PREFIX=EBxx
OUTPUTDIR=/path/to/results/<date>_${PREFIX}_Jellyfish
R1=/path/to/TellSeq_R1.fastq
R2=/path/to/TellSeq_R2.fastq

if [ ! -s "$R1" ]
then
	die "$R1 file not found"
fi
if [ ! -s "$R2" ]
then
	die "$R2 file not found"
fi
if [ ! -d "$OUTPUTDIR" ]
then
	mkdir -p $OUTPUTDIR
fi

echo `date`
echo ################################ Run Jellyfish

module load jellyfish/2.3.0-IGB-gcc-8.2.0

cd $OUTPUTDIR

ln -s $R1 ./
ln -s $R2 ./

## run step1 - jellyfish count

jellyfish count -C -m 21 -s 5G -t $SLURM_NPROCS  -o ${PREFIX}_bothTrimmedReads.jf *.fastq

## run step2 - jellyfish histo

jellyfish histo -t $SLURM_NPROCS ${PREFIX}_bothTrimmedReads.jf > ${PREFIX}_bothTrimmedReads.histo 

## the file ending in histo is the file that is used to run GenomeScope
