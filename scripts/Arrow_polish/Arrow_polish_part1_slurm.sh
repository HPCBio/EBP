#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=1000g
#SBATCH -N 1      #nodes
#SBATCH -n 36     #cpus-per-task
#SBATCH -p hpcbio #queue name
#SBATCH -A ebp    #account
#SBATCH --mail-type=ALL
#SBATCH -D /home/groups/earthbiogenome/results/20210224_EB20_80X_ArrowPolish
#SBATCH -J EB20_Arrow_polish_part1 #jobname
#SBATCH --mail-user=grendon@illinois.edu

# ----------------Load Modules--------------------

module load  smrtlink/9.0.0.92188


### edit lines 11-13 above and lines 38-41 below

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

if [ ! -d "$CLRdir" ]
then
	die "$CLRdir path not found"
fi
if [ ! -s "$GENOME" ]
then
	die "$GENOME file not found"
fi

if [ ! -d "$SCRATCH" ]
then
	mkdir -p $SCRATCH
fi

if [ ! -d "$outputdir" ]
then
	mkdir -p $outputdir
fi


echo `date`
echo ################################ Step 1: stage in files to scratch

cd $SCRATCH

export TMPDIR=$SCRATCH

cp  $CLRdir/${PREFIX}*.subreads.bam $SCRATCH
zcat -d $GENOME > ${SCRATCH}/$GENOMENAME
grep ">" $GENOMENAME | cut -d  ' ' -f1 | sed 's/>//' > $ListContigs

echo `date`
echo ################################ Step 2: create xml file with bam files of CLR reads

dataset create --type SubreadSet --name ${PREFIX}_CLR ${PREFIX}_CLR.xml  *.subreads.bam

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo dataset create cmd failed
exit $exitcode
fi

    
echo `date`
echo ################################ Step 3: align CLR reads to draft assembly

pbmm2 align \
--log-level INFO --preset SUBREAD \
--sort -m 5G -j $SLURM_NPROCS -J 8 -l 1000 \
$GENOMENAME \
${PREFIX}_CLR.xml \
$BAM \
 > ${PREFIX}_pbmm2.log 2>&1


exitcode=$?
if [ $exitcode -ne 0 ]
then
echo pbmm2 align cmd failed
exit $exitcode
fi

if [ ! -s $BAM ]
then
die "ERROR: $BAM file not created"
fi

echo `date`
echo ################################ Step 4: stage out alignment results

cp ${PREFIX}_pbmm2.log $outputdir
cp $ListContigs $outputdir
cp ${BAM}* $outputdir

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

