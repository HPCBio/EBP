#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=50g # memory requested
#SBATCH -N 1 #nodes
#SBATCH -n 24 #cpus-per-task
#SBATCH -A ebp #account
#SBATCH -p normal #queue
#SBATCH --mail-type=ALL
#SBATCH -J EBxx_ARKS_LINKS #jobname
#SBATCH -D /path/to/working/directory/
#SBATCH --mail-user=youremail@illinois.edu

# ----------------Load Modules--------------------

module load ARCS/1.2.1-IGB-gcc-8.2.0-Perl-5.28.1

### Modify lines 11-13 above 
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
### GENOME is the input genome for scaffolding
### inputARC is the fie of TellSeq reads formatted for ARCS
### PREFIX is the prefix of output files
### outputdir is the folder where results will be written to
### output_prefix string used by ARCS and LINKS

PREFIX=EBxx
output_prefix=somestring
GENOME=/path/to/assembly_to_scaffold.fa
outputdir=/home/groups/earthbiogenome/results/date_${PREFIX}_ARKS_LINKS
inputARC=/home/groups/earthbiogenome/data/TellSeq_reads/${PREFIX}*4ARCS.fastq.gz
ARCSscript=/home/apps/software/ARCS/1.2.1-IGB-gcc-8.2.0-Perl-5.28.1/Examples/makeTSVfile.py
GENOMENAME=$(basename $GENOME)

# ----------------Commands------------------------

# ----------------

# 1. Create an empty file for input to arcs

touch empty.fof

# ----------------

# 2. Softlink fasta file to working directory if needed.

ln -s $GENOME ./$GENOMENAME

# ----------------

# 3. Run arcs read alignment in ARKS mode: ARKS uses kmer alignment, which is much faster
# arcs requires a special interleaved fastq file from the Tell-Seq reads
# usage: arcs [Options] --arks -f <contig sequence file> <list of linked read files>

arcs -v -b output_prefix -k 60 --min_size=1000 -t $SLURM_NPROCS --arks -f $GENOMENAME 

# ----------------

# 4. run makeTSVfile.py to make properly formatted tab-delimited file for LINKS processing
# usage: makeTSVfile.py <graph_file=prefix_original.gv> <prefix_tigpair_checkpoint.tsv> <contig sequence file>

python $ARCSscript \
${output_prefix}_original.gv \
${output_prefix}.tigpair_checkpoint.tsv \
$GENOMENAME

# ----------------

# 5. run LINKS which does the scaffolding
# usage: LINKS -f <contig sequence file> -s <empty.fof> are required 
# -b is prefix of "*_original.gv file"
# -l  minimum number of links (k-mer pairs) to compute scaffold (default -l 5, optional)
# -a  maximum link ratio between two best contig pairs (default -a 0.3, optional) *higher values lead to least accurate scaffolding*
# -z  minimum contig length to consider for scaffolding (default -z 500, optional)

LINKS -f $GENOMENAME -s empty.fof -b $output_prefix -l 5 -a 0.9 -z 1000

# ----------------

# 6. don't forget to shorten the headers before running BUSCO.

awk '/^>/{$0=">scaffold"++i}1' ${output_prefix}.scaffolds.fa > ${output_prefix}.scaffolds_short_headers.fa
