#!/bin/bash
#SBATCH -n 24
#SBATCH --mem=50g
#SBATCH -p normal
#SBATCH -A ebp
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=grendon@illinois.edu
#SBATCH -J EB29_prothint_genemark_braker
#SBATCH -D /home/groups/earthbiogenome/results/20210430_EB29_GeneMark/

### edit lines 7-9

### next lines are included to help troubleshoot the run
### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -x
set -euo pipefail
IFS=$'\n\t'
echo `hostname`

### Declare variables
### These variables need editing
### GENOME is the genome after the Sealer and repeatMasker steps
### PREFIX is the prefix of output files

PREFIX=EB29
GENOME=/home/groups/earthbiogenome/results/20210319_EB29_RAILS_SEALER/EB29_Sealer_scaffold.fa.masked
NEWGENOME=${PREFIX}_shortHeaders.fa
OUPUTDIR=/home/groups/earthbiogenome/results/20210430_EB29_GeneMark/

### these variables do not need editing
PROT=/home/groups/earthbiogenome/data/OrthoDB/Arthropoda/clean_proteins.fasta
PROTDB=proteinDB.fa
SCRATCH=/scratch/$SLURM_JOB_NAME
PROTHINTDIR=${SCRATCH}/prothint
GENEMARKDIR=${SCRATCH}/genemark
BRAKERDIR=${SCRATCH}/braker
CUFFLINKSDIR=${SCRATCH}/cufflinks

### these cp cmd have to be executed only the first time that this module is executed

cp -fR /home/apps/software/GeneMark-ES/4.62-IGB-gcc-8.2.0-Perl-5.28.1/gm_key_64 ~/.gm_key_64
cp -fR /home/groups/earthbiogenome/src/.gm_key ~/.gm_key

### sanity check

if [ ! -s $GENOME ]
then
echo ERROR: $GENOME empty file
exit 1
fi

if [ ! -s $PROT ]
then
echo ERROR: $PROT empty file
exit 1
fi


### Step 1: reformat headers of assembly file
echo `date`
cd $OUPUTDIR

awk '/^>/{$0=">scaffold"++i}1' $GENOME > $NEWGENOME


if [ ! -s $NEWGENOME ]
then
echo ERROR: $NEWGENOME empty file
exit 1
fi


### Step 2: Stage files to scratch
echo `date`

if [ ! -d $SCRATCH ]
then
mkdir -p $PROTHINTDIR 
mkdir -p $GENEMARKDIR
mkdir -p $BRAKERDIR
mkdir -p $CUFFLINKSDIR
fi

cd $SCRATCH

if [ ! -s $NEWGENOME ]
then
ln -s $OUPUTDIR/$NEWGENOME ./$NEWGENOME
ln -s $PROT ./$PROTDB
fi


### Step 3: run prohint
echo `date`

module load ProtHint/2.5.0-IGB-gcc-8.2.0

prothint.py --workdir $PROTHINTDIR --threads $SLURM_NPROCS $NEWGENOME $PROTDB

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo prohint.py cmd failed
exit $exitcode
fi
if [ ! -s ${PROTHINTDIR}/prothint.gff ]
then
echo ERROR: ${PROTHINTDIR}/prothint.gff empty file
exit 1
fi

echo stage out prohint results to outputdir

cp -R $PROTHINTDIR $OUPUTDIR

### Step 4: run GeneMark-EP+
echo `date`
module purge
module load BRAKER/2.1.6-IGB-gcc-8.2.0
module load cdbfasta/20181005-IGB-gcc-8.2.0

gmes_petap.pl --EP $PROTHINTDIR/prothint.gff \
--evidence $PROTHINTDIR/evidence.gff \
--seq $NEWGENOME --soft_mask 1000 --verbose --cores $SLURM_NPROCS

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo  gmes_petap.pl cmd failed
exit $exitcode
fi
if [ ! -s genemark.gtf ]
then
echo ERROR: genemark.gtf empty file
exit 1
fi

## move genemark results to a folder for easy transfer

cp genemark.gtf ${GENEMARKDIR}/
cp gmes.log  ${GENEMARKDIR}/
cp gmhmm.mod ${GENEMARKDIR}/
cp run.cfg ${GENEMARKDIR}/
cp -R data ${GENEMARKDIR}/
cp -R run ${GENEMARKDIR}/
cp -R output ${GENEMARKDIR}/
cp -R info ${GENEMARKDIR}/

echo stage out genemark results to outputdir

cp -R $GENEMARKDIR $OUPUTDIR

### Step 5:  run BRAKER2
echo `date`

braker.pl --species=fly --useexisting \
--genome=$NEWGENOME \
--hints=${PROTHINTDIR}/prothint_augustus.gff \
--geneMarkGtf=${GENEMARKDIR}/genemark.gtf \
--softmasking  \
--AUGUSTUS_CONFIG_PATH=/home/groups/earthbiogenome/src/config/ \
--AUGUSTUS_BIN_PATH=/home/apps/software/augustus/3.3.3-IGB-gcc-8.2.0/bin/ \
--AUGUSTUS_SCRIPTS_PATH=/home/apps/software/augustus/3.3.3-IGB-gcc-8.2.0/scripts/ \
--CDBTOOLS_PATH=/home/apps/software/cdbfasta/20181005-IGB-gcc-8.2.0/bin/ \
--workingdir=$BRAKERDIR \
--cores=$SLURM_NPROCS

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo  braker.pl cmd failed
exit $exitcode
fi
if [ ! -s ${BRAKERDIR}/braker.gtf ]
then
echo ERROR: ${BRAKERDIR}/braker.gtf empty file
exit 1
fi

echo stage out braker results to outputdir

cp -R $BRAKERDIR $OUPUTDIR

### Step 5: run cufflinks to reformat gtf files and generate proteins.fa
echo `date`

module purge
module load Cufflinks/2.2.1

cd $CUFFLINKSDIR

ln -s ${BRAKERDIR}/braker.gtf ./
ln -s ${BRAKERDIR}/augustus.hints.gtf ./

gffread -o ${PREFIX}_augustus.hints.gff3 augustus.hints.gtf

gffread -o ${PREFIX}_braker.gff3 braker.gtf

gffread -g ${SCRATCH}/$NEWGENOME -y tmp.faa ${PREFIX}_braker.gff3

sed 's/\.$//g' tmp.faa > ${PREFIX}_proteins.fa
 
rm tmp.faa

if [ ! -s ${PREFIX}_proteins.fa ]
then
echo ERROR: ${PREFIX}_proteins.fa empty file
exit 1
fi


echo stage out braker results, reset and exit
echo `date`

cp -R $CUFFLINKSDIR $OUPUTDIR
cd $OUPUTDIR

#rm -rf $SCRATCH


