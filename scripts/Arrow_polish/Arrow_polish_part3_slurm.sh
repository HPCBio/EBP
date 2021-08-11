#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=1000g
#SBATCH -N 1      #nodes
#SBATCH -n 12     #cpus-per-task
#SBATCH -p hpcbio #queue name
#SBATCH -A ebp    #account
#SBATCH --mail-type=ALL
#SBATCH -D /home/groups/earthbiogenome/results/20210224_EB20_80X_ArrowPolish
#SBATCH -J EB20_Arrow_polish_part3 #jobname
#SBATCH --mail-user=grendon@illinois.edu


### edit lines 11-13 above and lines 34-38 below

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
### GENOME is the genome after the redbean step with full path
### CLRdir is the directory where the subread files are 
### PREFIX is the prefix of output files
### outputdir is the folder where results will be written to
 
PREFIX=EB20
outputdir=/home/groups/earthbiogenome/results/20210224_${PREFIX}_ArrowPolish
GENOME=/home/groups/earthbiogenome/results/20210204_EB20_redbean/EB20_redbean_g2.7_80X.ctg.fa.gz
CLRdir=/home/groups/earthbiogenome/data/PacBio_CLR_reads
GENOMENAME=${PREFIX}_inputgenome.fa
ListContigs=${PREFIX}_ListContigs.txt
BAM=${PREFIX}_aligned_sorted.bam
BAMdir=${outputdir}/tmp_polish
OUTPUTGENOME=${PREFIX}_arrowPolished.fa
SCRATCH=/scratch/$SLURM_JOB_NAME


if [ ! -s "$GENOME" ]
then
	die "$GENOME file not found"
fi
if [ ! -d "BAMdir" ]
then
	die "$BAMdir path not found"
fi
if [ ! -s "$outputdir/$ListContigs" ]
then
	die "$outputdir/$ListContigs file not found"
fi
if [ ! -d "$SCRATCH" ]
then
	mkdir -p $SCRATCH
fi

echo `date`
echo ################################ Step 1: stage in files to scratch

cd $SCRATCH

mkdir tmp_polish
cd tmp_polish

cp $outputdir/$ListContigs ./
cp $BAMdir/*fasta ./
truncate -s 0 cat.log

echo `date`
echo ################################ step 2 concatenate files

while read contig
do
input=$(ls ${contig}*fasta | tr '\n' ' ' )
echo start processing $input >> cat.log
cat $input >> $OUTPUTGENOME
echo done processing $input  >> cat.log
done < $ListContigs

if [ ! -s "$OUTPUTGENOME" ]
then
echo ERROR: $OUTPUTGENOME empty file
exit $exitcode
fi

cp $OUTPUTGENOME $SCRATCH

echo `date`
echo ################################ step 3 QC polished assembly


module purge
module load Perl/5.28.1-IGB-gcc-8.2.0
module load BUSCO/4.1.4-IGB-gcc-8.2.0-Python-3.7.2

cd $SCRATCH

perl -I /home/groups/earthbiogenome/src/assemblathon2-analysis/ \
/home/groups/earthbiogenome/src/assemblathon2-analysis/assemblathon_stats.pl \
$OUTPUTGENOME > ${OUTPUTGENOME}.stats

################## run BUSCO ##########################

export AUGUSTUS_CONFIG_PATH=/home/groups/earthbiogenome/src/config/

export BUSCO_CONFIG_FILE=/home/groups/earthbiogenome/src/config.ini

busco --offline -c $SLURM_NPROCS -i $OUTPUTGENOME \
-o BUSCO_${PREFIX}_arrowPolished \
-l /home/apps/software/BUSCO/4.1.4-IGB-gcc-8.2.0-Python-3.7.2/data/lineages/endopterygota_odb10 \
-m genome -f 

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo busco failed
exit $exitcode
fi

echo `date`
echo ################################ Step 4: stage out  results

buscodir=$(ls ${BUSCO}* | tr '\n' ' ' )
cp ${OUTPUTGENOME}* $outputdir
cp -R $buscodir $outputdir

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo stage out failed
exit $exitcode
fi

echo `date`
echo ################################ done. reset and exit

cd $outputdir
#rm -rf $SCRATCH


