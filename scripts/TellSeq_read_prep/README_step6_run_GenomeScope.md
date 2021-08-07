# Background 

TellSeq libraries: https://www.universalsequencing.com/technology


TellSeq libraries produce short reads.  These short reads were used in this project for two different purposes: 

1. To estimate genome size using kmer distribution histograms of the short reads and then running GenomeScope 2.0 http://qb.cshl.edu/genomescope/genomescope2.0/

2. To scaffold the assembly with short reads.  The scaffolding toolkit we used was ARCS-LINKS. See page: https://github.com/bcgsc/arcs



# Dependencies

We used jellyfish to calculate the kmer distribution histogram.

- jellyfish version 2.3.0 or higher http://www.genome.umd.edu/jellyfish.html#Release

- slurm workload manager https://slurm.schedmd.com/documentation.html

# Edit the script

The script jellyfish_slurm.sh runs jellyfish on a linux cluster. 

Please edit the script with your information.

Notice that jellyfish expects the input to be uncompressed fastq files.


<pre>
#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=100g # memory requested
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -A ebp #account
#SBATCH -p normal #queue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@illinois.edu
#SBATCH -J EBxx_jellyfish #jobname
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
### R1 and R2 are the files of TellSeq reads
### OUTPUTDIR is the folder where results will be written to
 
PREFIX=EBxx
OUTPUTDIR=/path/to/results/<date>_${PREFIX}_Jellyfish
R1=/path/to/TellSeq_R1.fastq
R2=/path/to/TellSeq_R2.fastq

if [ ! -s "$R1" ]
then
	die "$R1 file not found"
fi
if [ ! -s "$R2" ]
then
	die "$R2 file not found"
fi
if [ ! -d "$OUTPUTDIR" ]
then
	mkdir -p $OUTPUTDIR
fi

echo `date`
echo ################################ Run Jellyfish

module load jellyfish/2.3.0-IGB-gcc-8.2.0

cd $OUTPUTDIR

ln -s $R1 ./
ln -s $R2 ./

## run step1 - jellyfish count

jellyfish count -C -m 21 -s 5G -t $SLURM_NPROCS  -o ${PREFIX}_bothTrimmedReads.jf *.fastq

## run step2 - jellyfish histo

jellyfish histo -t $SLURM_NPROCS ${PREFIX}_bothTrimmedReads.jf > ${PREFIX}_bothTrimmedReads.histo 

## the file ending in histo is the file that is used to run GenomeScope

</pre>

# Run the script

<pre>
 sbatch jellyfish_slurm.sh
</pre>

# Output files

<pre>
${PREFIX}_bothTrimmedReads.jf

${PREFIX}_bothTrimmedReads.histo 
</pre>

The file ending in histo is the file needed to run GenomeScope

Transfer that file to a local folder.

# Run GenomeScope

Open this link in a browser: http://qb.cshl.edu/genomescope/genomescope2.0/

Drop the histo file in the space provided and then click on the Submit button



