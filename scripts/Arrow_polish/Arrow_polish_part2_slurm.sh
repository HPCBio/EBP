#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=1000g
#SBATCH -N 1      #nodes
#SBATCH -n 12     #cpus-per-task
#SBATCH -p hpcbio #queue name
#SBATCH -A ebp    #account
#SBATCH --mail-type=ALL
#SBATCH -D /home/groups/earthbiogenome/results/20210224_EB20_80X_ArrowPolish
#SBATCH -J EB20_Arrow_polish_part2 #jobname
#SBATCH --mail-user=grendon@illinois.edu
#SBATCH --array=1-totcontigs%4

### edit lines 11-14 above and lines 34-38 below
### totcontigs is the number of contigs in the input assembly
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
### ListContigs is a text file; each line is a list of contigs that will go into the same chunk
 
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
if [ ! -d "$outputdir" ]
then
	die "$outputdir path not found"
fi
if [ ! -s "$outputdir/$ListContigs" ]
then
	die "$outputdir/$ListContigs file not found"
fi
if [ ! -s "$outputdir/$BAM" ]
then
	die "$outputdir/$BAM file not found"
fi
if [ ! -d "$BAMdir" ]
then
	mkdir -p $BAMdir
fi
if [ ! -d "$SCRATCH" ]
then
	mkdir -p $SCRATCH
fi

echo `date`
echo ################################ Step 1: stage in files to scratch

cd $SCRATCH

export TMPDIR=$SCRATCH

if [ ! -s "$BAM" ]
then
cp $outputdir/$BAM $SCRATCH/
fi

if [ ! -s "$ListContigs" ]
then
cp $outputdir/$ListContigs $SCRATCH/
fi

contig=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" $ListContigs)
contigBAM=${contig}_input.bam
contigOUT=${contig}_arrow_polished.fasta

echo `date`
echo ################################ step 2 extract bam for this chunck only

## example: samtools view -b -o   chunk_101.bam  EB20_redbean_g2.7_80X_sorted.bam  ctg3  ctg4  ctg5  ctg6

module load SAMtools/1.12-IGB-gcc-8.2.0

samtools view -b -o $contigBAM $BAM $contig


exitcode=$?
if [ $exitcode -ne 0 ]
then
echo samtools view cmd failed for $contig
exit $exitcode
fi

echo `date`
echo ################################ step 3 polish with arrow


module purge
module load smrtlink/10.0.0.108728

gcpp --algorithm=arrow --log-level INFO --log-file ${contig}_gcpp.log -j $SLURM_NPROCS \
$contigBAM -r  $GENOMENAME -o $contigOUT

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo gcpp cmd failed for $contig
exit $exitcode
fi

echo `date`
echo ################################ Step 4: stage out alignment results

cp ${contig}_gcpp.log $BAMdir
cp ${contigBAM}* $BAMdir
cp ${contigOUT}* $BAMdir

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


