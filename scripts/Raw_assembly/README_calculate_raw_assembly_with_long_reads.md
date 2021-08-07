# Background

We used PacBio long reads to generate the raw assembly.

There are two kinds of PacBio long reads:

HiFi reads. see https://www.pacb.com/smrt-science/smrt-sequencing/hifi-reads-for-highly-accurate-long-read-sequencing/

CLR (continuous long reads) See https://www.pacb.com/tag/continuous-long-read-clr/

There are many assembly tools for these types of long reads.

We used hifiasm to generate raw assemblies with HiFi reads. https://github.com/chhylp123/hifiasm

We used Redbean (aka wtdbg2) to generate raw assemblies with CLR reads. 
https://github.com/ruanjue/wtdbg2
https://github.com/ruanjue/wtdbg2/blob/master/README-ori.md

# Dependencies

- hifiasm version 0.13 or higher

- gfatools version 0.4  or higher

- wtdbg2 version 2.5 or higher

- seqkit version 0.11 or higher

- GenomeScope version 2 http://qb.cshl.edu/genomescope/genomescope2.0/

- slurm workload manager https://slurm.schedmd.com/documentation.html

# 1. Generate raw assembly using hifiasm and HiFi reads

We used hifiasm to generate raw assemblies with HiFi reads. https://github.com/chhylp123/hifiasm


## 1.1 The hifiasm script

The script hifiasm_slurm.sh runs hifiasm on a linux cluster. 

Please edit the script with your information.

The script consists of two steps:

- Step 1 runs hifiasm and generates the raw assembly. The output is NOT in fasta format.

- Step 2 runs gfatools gfa2fa to reformat the output and write it in fasta format.

## 1.2 Run the script

<pre>
 sbatch hifiasm_slurm.sh
</pre>

The tool is very fast. It should not take more than 6 hours to run.


## 1.3 hifiasm output files

Hifiasm generates different types of assemblies based on the input data. 

It also writes error corrected reads to the prefix.ec.bin binary file and writes overlaps to prefix.ovlp.source.bin and prefix.ovlp.reverse.bin. 

For more details, please see the complete: https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output



# 2. Generate raw assembly using redbean and CLR reads

This tool requires some preprocessing steps.

- Estimate genome size

- Filter CLR reads by length min=5000 max=library-length

- Estimate read coverage


## 2.1 Estimate genome size

See README_step6_run_GenomeScope.md  for full details

Example:

<pre>

Results

GenomeScope version 2.0
input file = user_uploads/UFE4kh5TB9DozL0NkKCq
output directory = user_data/UFE4kh5TB9DozL0NkKCq
p = 2
k = 21

property                      min               max               
Homozygous (aa)               95.9674%          100%              
Heterozygous (ab)             0%                4.03258%          
Genome Haploid Length         1,709,133,088 bp  3,291,686,136 bp  
Genome Repeat Length          723,339,182 bp    1,393,107,168 bp  
Genome Unique Length          985,793,907 bp    1,898,578,968 bp  
Model Fit                     62.3549%          97.7983%          
Read Error Rate               0.303731%         0.303731%

</pre>

The estimated genome size is 3.3Gbp

## 2.2 Calculate stats of CLR reads

There are many tools that can do this job. We used seqkit stats

<pre>

module load seqkit

export SEQKIT_THREADS=2

seqkit stats -j 24 -T *subreads.fastq > read_fate.tsv

cat read_fate.tv

file 	num_seqs 	sum_len 	min_len 	avg_len 	max_len
EB20_m64108_210117_053528.subreads.fastq 	27,570,320 	232,559,004,086 	50 	8,435 	521,818
EB20_m64108_210118_121306.subreads.fastq 	28,535,669 	240,060,443,237 	50 	8,413 	436,859

</pre>

## 2.3 Filter CLR reads by read length 

Denovo assembly tools select the longest reads as *anchors".

We need to make sure that the longest reads are of high quality. 

One way to achieve that goal is to filter out reads that are longer than the length of the library used for sequencing. In our case that would be 60kb.

Assembly tools like redbean filter reads shorter than 5kb anyway; but we may need to switch to a different tool that does not filter reads by length at all. 

Therefore, it is better to do it now.

This is the cmd to filter reads shorter than 5kb and longer than 60kb

<pre>

$ module load seqkit

$ seqkit seq -w 0 -m 5000 -M 60000 EB20_m64108_210117_053528.subreads.fastq  > EB20_m64108_210117_053528.subreads_min5kb_max50kb.fastq &

$ seqkit seq -w 0 -m 5000 -M 60000 EB20_m64108_210118_121306.subreads.fastq  > EB20_m64108_210118_121306.subreads_min5kb_max60kb.fastq &

</pre>

## 2.4 Check read lengths again

<pre>

$ module load seqkit

$ seqkit stats -j 24 *.fastq -T > CLR_stats.tsv

$ cat CLR_stats.tsv

file    format    type    num_seqs    sum_len    min_len    avg_len    max_len
EB20_m64108_210117_053528.subreads_min5kb_max60kb.fastq    FASTQ    DNA    12505818    1.85924E+11    5000    14867    60000
EB20_m64108_210118_121306.subreads_min5kb_max60kb.fastq    FASTQ    DNA    13251586    1.95534E+11    5000    14755.5    60000

</pre>

## 2.5 Calculate read depth

The estimated genome size is 3.3Gbp

For 90x coverage we would need = 90 * (3 300 000 000)
For 80x coverage we would need = 80 * (3 300 000 000)
And so on

Do we have enough read coverage for each one of these levels? The answer is yes we do

## 2.6 Edit the script

The script run_redbean_slurm.sh runs the assembly tool on a cluster.

Please edit the script with your information.

The script consists of two steps:

- Step 1 runs wtdbg2 and generates the raw assembly. The output is NOT in fasta format.

- Step 2 runs wtpoa-cns to reformat the output and write it in fasta format.

## 2.7 Run the script

<pre>
 sbatch run_redbean_slurm.sh
</pre>

This tool takes a very long time to execute.  It could take over a week or more to complete.

## 2.8 Output files

This link explains the output files of wtdbg2 : https://github.com/ruanjue/wtdbg2/blob/master/README-ori.md

The raw assembly should be in the file ending in ctg.fa




