#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=10g # memory requested
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -A ebp #account
#SBATCH -p normal #queue
#SBATCH --mail-type=ALL
#SBATCH -J EBxx_BUSCO_protein
#SBATCH --mail-user=your_email@illinois.edu
#SBATCH -D /home/groups/earthbiogenome/results/path/to/your/working/directory/

### edit lines 11-13 above and lines 31-32 below

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
### OUTPUTDIR specify the path to the folder with GeneMark results as shown below
### PROT is the input file for the BUSCO cmd. It should be in the cufflinks folders as shown below

OUTPUTDIR=/home/groups/earthbiogenome/results/20210430_EB29_GeneMark/
PROT=$OUTPUTDIR/cufflinks/EB29_proteins.fa
STATS=$OUTPUTDIR/cufflinks/EB29_proteins.fa.stats

if [ ! -d "$OUTPUTDIR" ]
then
	die "$OUTPUTDIR path not found"
fi
if [ ! -s "$PROT" ]
then
	die "$PROT file not found"
fi

# ----------------Load Modules--------------------

module load seqkit/0.12.1
module load BUSCO/4.1.4-IGB-gcc-8.2.0-Python-3.7.2
export AUGUSTUS_CONFIG_PATH=/home/groups/earthbiogenome/src/config/
export BUSCO_CONFIG_FILE=/home/groups/earthbiogenome/src/config.ini

# ----------------Run the command--------------------

cd $OUTPUTDIR

seqkit stats -j $SLURM_NPROCS -T $PROT > $STATS

busco  \
-i $PROT \
-o BUSCO_proteins_vs_insecta \
-l /home/apps/software/BUSCO/4.1.4-IGB-gcc-8.2.0-Python-3.7.2/data/lineages/insecta_odb10 \
-m protein -f --offline -c $SLURM_NPROCS

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo busco failed
exit $exitcode
fi

