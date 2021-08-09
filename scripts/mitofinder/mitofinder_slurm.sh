#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=50g  #memory requested
#SBATCH -N 1       #nodes
#SBATCH -n 4       #cpus-per-task
#SBATCH -p normal  #queue
#SBATCH -A ebp     #account
#SBATCH --mail-type=ALL #get all job notifications
#SBATCH --mail-user=youremail@illinois.edu #email address
#SBATCH -J EBxx_MitoFinder #jobname
#SBATCH -D /working/dir

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
### GENOME1 is the polised genome after the sealer step with full path
### GENOME2 is the raw genome after the hifiasm step with full path
### PREFIX1 is the prefix of output files for the first mitofinder cmd with GENOME1
### PREFIX2 is the prefix of output files for the second mitofinder cmd with GENOME2
### SPECIES is the genome file in genbank format of the closes species to the genome

PREFIX1=EB31_MitoFinder_polished_genome
PREFIX2=EB31_MitoFinder_raw_genome
GENOME1=/full/path/to/EB31_Sealer_scaffold.fa
GENOME2=/full/path/to/EB31_hifiasm.p_ctg.fasta
SPECIES=/full/path/NC_053918.1_Parazyginella_tiani_mitochondrion.gb
OUTPUTDIR=/home/groups/earthbiogenome/results/date_EBxx_MitoFinder
GENOMENAME1=$(basename $GENOME1)
GENOMENAME2=$(basename $GENOME2)
SPECIESNAME=$(basename $SPECIES)

if [ ! -s "$GENOME1" ]
then
	die "$GENOME1 file not found"
fi
if [ ! -s "$GENOME2" ]
then
	die "$GENOME2 file not found"
fi
if [ ! -s "$SPECIES" ]
then
	die "$SPECIES file not found"
fi
if [ ! -d "$OUTPUTDIR" ]
then
	mkdir -p $OUTPUTDIR
fi

cd $OUTPUTDIR
ln -s $GENOME1 ./$GENOMENAME1
ln -s $GENOME2 ./$GENOMENAME2
ln -s $SPECIES ./$SPECIESNAME


echo `date`
echo ################################ Step 1: run MitoFinder  with $GENOME1

module load singularity

singularity run /home/groups/hpcbio/singularity/MitoFinder_v1.4 \
 -j $PREFIX1 -a $GENOME1 \
 -r $SPECIESNAME \
 -p $SLURM_NPROCS -o 5 --max-contig-size 50000

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo MitoFinder cmd failed with $GENOME1
exit $exitcode
fi

echo `date`
echo ################################ Step 2: run MitoFinder  with $GENOME2

singularity run /home/groups/hpcbio/singularity/MitoFinder_v1.4 \
 -j $PREFIX2 -a $GENOME2 \
 -r $SPECIESNAME \
 -p $SLURM_NPROCS -o 5 --max-contig-size 50000

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo MitoFinder cmd failed with $GENOME1
exit $exitcode
fi
