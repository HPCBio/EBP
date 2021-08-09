#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=50g #memory requested
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -p normal #queue
#SBATCH -A ebp
#SBATCH --mail-user=kosterbu@illinois.edu #email address
#SBATCH --mail-type=ALL #get all job notifications
#SBATCH -J EBxx_protein_annotation #jobname
#SBATCH -D /some/folder/with/GeneMark/results

PREFIX=EBxx
OUTPUTDIR=/some/folder/with/GeneMark/results
PROTEINS=$OUTPUTDIR/cufflinks/${PREFIX}_proteins.fa
GFF3=$OUTPUTDIR/cufflinks/${PREFIX}_braker.gff3
UNIPROTDIR=/home/groups/earthbiogenome/data/UniprotDB/uniprot
SCRATCH=/scratch/$SLURM_JOB_NAME
INPROTEIN=$(basename $PROTEINS)
INGFF3=$(basename $GFF3)
SCRIPTDIR=/home/groups/earthbiogenome/src

# ----------------Commands------------------------

# stage data to /scratch


mkdir -p $SCRATCH
cd $SCRATCH
cp $PROTEINS $SCRATCH
cp $GFF3 $SCRATCH
cp $UNIPROTDIR/* $SCRATCH



# ----------------InterProScan protein domains----

# ----------------Load Modules--------------------

module load InterProScan/5.47-82.0-IGB-gcc-8.2.0-Java-15.0.1

# ----------------Commands------------------------

interproscan.sh --tempdir $SCRATCH \
--cpu $SLURM_NPROCS -appl pfam -dp -f TSV -goterms -iprlookup -pa -t p \
-i $INPROTEIN -o ${PREFIX}_output.iprscan

cp ${PREFIX}_output.iprscan $OUTPUTDIR

# ----------------SwissProt blastp to add functional information

# ----------------Load Modules--------------------

module purge
module load BLAST+/2.10.1-IGB-gcc-8.2.0

# ----------------Commands------------------------

blastp -num_threads $SLURM_NPROCS -query $INPROTEIN \
-db uniprot_sprot.db -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 \
-out ${PREFIX}_uniprot_sprot_output.blastp

cp ${PREFIX}_output.blastp $OUTPUTDIR

# ----------------Load Modules--------------------

module purge
module load Perl


# Add blastp UniProt info to fasta file

$SCRIPTDIR/maker_functional_fasta \
uniprot_sprot.fasta \
${PREFIX}_uniprot_sprot_output.blastp \
$INPROTEIN > ${PREFIX}_proteins.putative_function.fasta

# Add blastp UniProt info to gff file.

$SCRIPTDIR/maker_functional_gff \
uniprot_sprot.fasta \
${PREFIX}_uniprot_sprot_output.blastp \
$INGFF3 > ${PREFIX}_proteins.putative_function.gff

# Add InterProScan protein domain information to the gff file.

# Change "geneID=" to "Name=" in the gff file so that the next script works properly.

sed 's/geneID=/Name=/g' ${PREFIX}_proteins.putative_function.gff > ${PREFIX}_proteins.putative_function_rename.gff

$SCRIPTDIR/ipr_update_gff \
${PREFIX}_proteins.putative_function_rename.gff \
${PREFIX}_output.iprscan > ${PREFIX}_proteins.putative_function.domain_added.gff

cp ${PREFIX}_proteins.putative_function* $OUTPUTDIR
