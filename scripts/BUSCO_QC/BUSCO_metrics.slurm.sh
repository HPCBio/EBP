#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=10g # memory requested
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -A ebp #account
#SBATCH -p normal #queue
#SBATCH --mail-type=ALL
#SBATCH -J EBxx_BUSCO_metrics #jobname
#SBATCH --mail-user=your_email@illinois.edu
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
### GENOME the filename of the genome
### OUTPUTDIR is the folder where results will be written to


PREFIX=EBxx
GENOME=/path/to/genome.fasta
GENOMENAME=$(basename $GENOME)
OUTPUTDIR=/path/to/outputdir
BUSCODB=/path/to/BUSCO/data/lineages/insecta_odb10
SCRIPT=src/assemblathon/assemblathon_stats.pl
SCRIPTDIR=$(dirname $SCRIPT)
CONFIGDIR=src/config/
BUSCOCONFIG=src/config.ini
BUSCOPREFIX=BUSCO_$GENOMENAME

if [ ! -s "$GENOME" ]
then
	die "$GENOME file not found"
fi
if [ ! -d "$OUTPUTDIR" ]
then
	mkdir -p $OUTPUTDIR
fi


echo `date`
echo ################################ Step 1: run assemblathon

module load Perl/5.28.1-IGB-gcc-8.2.0

module load BUSCO/4.1.4-IGB-gcc-8.2.0-Python-3.7.2

cd $OUTPUTDIR

ln -s $GENOME ./

perl -I $SCRIPTDIR $SCRIPT $GENOMENAME > ${GENOMENAME}.stats

echo `date`
echo ################################ Step 2: run BUSCO

export AUGUSTUS_CONFIG_PATH=$CONFIGDIR

export BUSCO_CONFIG_FILE=$BUSCOCONFIG

busco --offline -c $SLURM_NPROCS -i $GENOMENAME -o $BUSCOPREFIX -l $BUSCODB -m genome -f 