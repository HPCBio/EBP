# 1. Background

We performed two types of scaffolding: Scaffolding with TellSeq reads using ARCS-LINKS followed by Scaffolding with Hi-C reads using SALSA

## 1.1  Scaffolding with TellSeq reads using ARCS-LINKS

<pre>
Transposase Enzyme Linked Long-read Sequencing (TELL-Seq™) is a simple and scalable NGS library technology 
that generates barcode linked-reads for genome scale sequencing applications.
</pre>

Source: 

https://www.universalsequencing.com/technology

<pre>
ARCS maps the linked reads to the assembled contigs and constructs intercontig links 
by identifying pairs of contigs whose ends share sequences from the same read pool. 
Scaffolding is then performed with the LINKS scaffolder, 
a tool originally developed for scaffolding assemblies with the help of long-read data. 
A similar approach is used by ARKS, a tool that relies on k-mer matches 
instead of sequence alignment to infer the assignment of the linked reads to assembled contigs
</pre>

Sources:

Modern technologies and algorithms for scaffolding assembled genomes
Jay Ghurye, Mihai Pop
PLOS Published: June 5, 2019
https://doi.org/10.1371/journal.pcbi.1006994 

ARCS-LINKS: https://github.com/bcgsc/arcs


## 1.2 Scaffolding with Hi-C reads using SALSA

<pre>
HiC is a special type of paired-read data that are used to study the three-dimensional structure of chromosomes inside a cell...
In the Hi-C protocol, DNA in the cell nucleus is cross-linked and cut with a restriction enzyme. 
This process generates fragments of DNA that are distally located but physically associated with each other. 
The sticky ends of these fragments are biotinylated and then ligated to form a chimeric circle. 
The resulting circles are sheared and processed into sequencing libraries 
in which individual templates are chimeras of the physically associated DNA molecules.
</pre>


Sources:

Modern technologies and algorithms for scaffolding assembled genomes
Jay Ghurye, Mihai Pop
PLOS Published: June 5, 2019
https://doi.org/10.1371/journal.pcbi.1006994 

SALSA  https://github.com/marbl/SALSA

<pre>
SALSA2 begins with a draft assembly generated from long reads such as Pacific Biosciences [23] or Oxford Nanopore [24]. 
SALSA2 requires the unitig sequences and, optionally, a GFA-formatted assembly graph [25] representing the ambiguous reconstructions. 
Hi-C reads are aligned to the unitig sequences, and unitigs are optionally split in regions lacking Hi-C coverage. 
A hybrid scaffold graph is constructed using both ambiguous edges from the GFA and edges from the Hi-C reads, scoring edges according to a “best buddy” scheme. 
Scaffolds are iteratively constructed from this graph using a greedy weighted maximum matching. 
A mis-join detection step is performed after each iteration to check if any of the joins made during this round are incorrect. 
Incorrect joins are broken and the edges blacklisted during subsequent iterations. 
This process continues until the majority of joins made in the prior iteration are incorrect. 
This provides a natural stopping condition, when accurate Hi-C links have been exhausted. 
</pre>

Source:

https://www.biorxiv.org/content/10.1101/261149v2.full

# 2. Dependencies

- This step should be run after the purging duplicates step of a HiFi assembly
OR after the pilon polishing step of a CLR assembly

- ARCS 1.2.1 https://github.com/bcgsc/arcs

- LINKS  1.8

- makeTSVfile.py  a script that is provided in this repo

- BEDTools 2.28.0

- SAMtools 1.11

- SALSA 20191001 https://github.com/marbl/SALSA

- BWA 0.7.17

- picard 2.10.1

- slurm workload manager https://slurm.schedmd.com/documentation.html

# 3. Scaffolding with TellSeq reads

## 3.1 Setup


- Copy makeTSVfile.py to your path

- The file of TellSeq reads has to be interleaved and each header has to have the correct tag
as described here: https://github.com/bcgsc/arcs section *Using stLFR linked reads*

## 3.2 The steps

- Create an empty file empty.fof

- Run arcs read alignment in ARKS mode

- Run makeTSVfile.py to generate checkpoints

- Run LINKS to apply the scaffolding information to the genome 

<pre>

$ touch empty.fof

$ arcs -v \
  -b output_prefix -k 60 --min_size=1000 \
  -t $treads \
  --arks -f $GENOMENAME $TellSeq4ARCS

$ python makeTSVfile.py \
  ${output_prefix}_original.gv \
  ${output_prefix}.tigpair_checkpoint.tsv \
  $GENOMENAME

$ LINKS -f $GENOMENAME \
  -s empty.fof -b $output_prefix \
  -l 5 -a 0.9 -z 1000

</pre>

## 3.3 Run the script


The script that runs the analysis is called ARKS_LINKS.slurm.sh

Edit the script with the information of your genome and corresponding TellSeq reads.


<pre>
sbatch ARKS_LINKS.slurm.sh
</pre>

The script can take anywhere between 2 hours - 4 hours.

## 3.4 Output files

The output files of the script look like this: https://github.com/bcgsc/arcs/tree/master/Examples/arks_test-demo/output

## 3.5 QC the ARCS-scaffolded assembly

We provide a script BUSCO_metrics_slurm.sh that runs assemblathon.pl to calculate assembly contiguity.

It also runs BUSCO (database: insecta_odb10)  to calculate genome completeness.

# 4. Scaffolding with Omni-C using SALSA

This analysis is performed in two steps:

## 4.1 The mapping step using the arima pipeline. 

- Generate genome.fai

- bwa-index genome
 
- bwa-align TellSeq reads to assembly > aln.bam

- Sort aln.bam  > sorted_aln.bam

- Mark and remove PRC duplicates with picard tools > dedup.bam

- Generate index > dedup.bam.bai


## 4.2 The scaffolding step using SALSA

- Generate bed file from dedup.bam file with BEDTools  > dedup.bed

- Sort bed file by read name -> sorted_dedup.bed

- Run SALSA.

<pre>

$ run_pipeline.py \
 -a $GENOMENAME \
 -l ${GENOMENAME}.fai \
 -b $sorted_dedup.bed \
 -e DNASE -o scaffolds \
 -m yes

</pre>

- If you also have the genome in graph format, then you can add it to the SALSA command with the -g option


## 4.3 QC the SALSA-scaffolded assembly

We provide a script BUSCO_metrics_slurm.sh that runs assemblathon.pl to calculate assembly contiguity.

It also runs BUSCO (database: insecta_odb10)  to calculate genome completeness.












