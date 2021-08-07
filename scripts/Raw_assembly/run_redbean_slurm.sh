#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH --mem=750g
#SBATCH -N 1 #nodes
#SBATCH -n 36 #cpus-per-task
#SBATCH -A ebp #account
#SBATCH -p normal #queue
#SBATCH --mail-type=ALL
#SBATCH -J EBxx_Redbean_80X #jobname
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
### CLR1.. CLRn are the files of CLR reads
### OUTPUTDIR is the folder where results will be written to
### DEPTH is the read depth
### SIZE is the estimated genome size

PREFIX=EBxx
CLR1=/path/to/CLR_set1_subreads.fastq
CLR2=/path/to/CLR_set1_subreads.fastq
DEPTH=80x
SIZE=2.5g
OUTPUTDIR=/path/to/results/<date>_${PREFIX}_Redbean_${SIZE}_${DEPTH}
SCRATCH=/scratch/$SLURM_JOB_NAME

if [ ! -s "$CLR1" ]
then
	die "$CLR1 file not found"
fi
if [ ! -s "$CLR2" ]
then
	die "$CLR2 file not found"
fi
if [ ! -d "$OUTPUTDIR" ]
then
	mkdir -p $OUTPUTDIR
fi
if [ ! -d "$SCRATCH" ]
then
	mkdir -p $SCRATCH
fi

echo `date`
echo ################################ Step 1: stage in files to scratch

cd $SCRATCH

export TMPDIR=$SCRATCH

cp  $CLR1 $SCRATCH
cp  $CLR2 $SCRATCH


echo `date`
echo ################################ Step2: run wtdbg2


# ----------------Load Modules--------------------

module load wtdbg2/2.5-IGB-gcc-4.9.4


wtdbg2 -x sq -g $SIZE -X $DEPTH -t $SLURM_NPROCS \
 -i $CLR1 \
 -i $CLR2 \
 -fo ${PREFIX}_Redbean_${SIZE}_${DEPTH}

echo `date`
echo ################################ Step3: run wtpoa

wtpoa-cns -t $SLURM_NPROCS -i ${PREFIX}_Redbean_${SIZE}_${DEPTH}.ctg.lay.gz -fo ${PREFIX}_Redbean_${SIZE}_${DEPTH}.ctg.fa

echo `date`
echo ################################ Step4: Stage out files to outputdir

cp ${PREFIX}_Redbean_* $OUTPUTDIR