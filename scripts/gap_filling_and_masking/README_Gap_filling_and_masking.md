# 1. Background

Unknown sequences, or gaps, are present in many published genomes across public databases. 
Gap filling is an important finishing step in de novo genome assembly, especially in large genomes. 
The gap filling problem is nontrivial and while there are many computational tools partially solving the problem, 
several have shortcomings as to the reliability and correctness of the output, i.e. the gap filled draft genome.

We performed two types of gap filling: 

- Gap-filling  using RAILS-COBBLER when HiFi reads are available

- Gap-filling  using SEALER when  Hi-C reads are available

We perform one type of genome masking:

- Masking with RepeatMasker. 


## 1.1 Gap-filling  using RAILS-COBBLER when HiFi reads are available.

<pre>
RAILS and Cobbler are genomics application for scaffolding and automated finishing of genome assemblies with long DNA sequences. 
They can be used to scaffold & finish high-quality draft genome assemblies with any long, preferably high-quality, sequences 
such as scaftigs/contigs from another genome draft.

They both rely on accurate, long DNA sequences to patch gaps in existing genome assembly drafts.

Cobbler is a utility to automatically patch gaps (ambiguous regions in a draft assembly, represented by N's) 
It does so by first aligning the long sequences to the assembly, 
tallying the alignments and replacing N's with the sequences from these long DNA sequences.

RAILS is an all-in-one scaffolder and gap-filler. Its process is similar to that of Cobbler. 
It scaffolds your genome draft with the help of long DNA sequences 
(contig sequences are ordered/oriented using alignment information). 
The newly created gaps are automatically filled with the DNA sequence of the provided long DNA sequence.
</pre>

Source:

https://github.com/bcgsc/RAILS


## 1.2 Gap-filling  using SEALER when  Hi-C reads are available

<pre>
Sealer is a scalable tool designed to close gaps within assembly scaffolds. 
Sealer is a stand-alone application of Konnector, a fast and low-memory Bloom filter-based de bruijn graph assembler. 
It identifies problematic regions (those with Ns) in a user-input scaffold file, 
targets the gaps for assembly iteratively using a range of k values and 
consolidates assembled fragments into a new draft genome assembly.
</pre>

Source:

https://www.bcgsc.ca/resources/software/sealer

## 1.3 Order of gap-filling

If both types of reads are availble; then perform the RAILS-COBBLER step first followed by the SEALER step.

## 1.4 Genome Masking

Genome masking is the step that identifies areas of the genome that should be masked after identifying them as 

- interspersed repeats, 

- low complexity regions,

- transposable elements and/or 

- endogenous retroviruses,

Hard masking replaces stretches of the genome with Ns 

Soft masking replaces Upper-case nucleotides with their lower-case counterparts.

Masking a genome is a pre-requisite to gene predictions and genome functional annotation.

See: https://www.repeatmasker.org/


# 2. Dependencies

- This step should be run after purging duplicates and scaffolding.

- Fastq2Fasta.pl a script that is provided in this repo

- seqkit 0.12.1 

- RAILS 1.5.1

- BWA 0.7.17

- picard 2.10.1

- abyss 2.2.5

- RepeatMasker https://www.repeatmasker.org/


- slurm workload manager https://slurm.schedmd.com/documentation.html

## 3.1 Gap-filling  using RAILS-COBBLER and HiFi reads

The script is called RAILS.slurm.sh

The execution could take several hours up to 2 days.

## 3.2 Gap-filling  using SEALER and  Hi-C reads

The script is called SEALER.slurm.sh

The execution could take approx a couple of hours.

## 3.3 Masking

The script is called RepeatMasker.slurm.sh

The execution could take up to 3 or 4 days.

A summary report looks like this:

<pre>
==================================================
file name: EB17_Sealer_OmniC_TellSeq_scaffold.fa
sequences:           146
total length: 1471763525 bp  (1471715325 bp excl N/X-runs)
GC level:         32.77 %
bases masked:  208386386 bp ( 14.16 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
Retroelements       243353     72653975 bp    4.94 %
   SINEs:            32339      2277935 bp    0.15 %
   Penelope           3813       325195 bp    0.02 %
   LINEs:           184446     63041511 bp    4.28 %
    CRE/SLACS            2           87 bp    0.00 %
     L2/CR1/Rex      15095      3164888 bp    0.22 %
     R1/LOA/Jockey   16300      2933565 bp    0.20 %
     R2/R4/NeSL        253        52170 bp    0.00 %
     RTE/Bov-B      100781     39156646 bp    2.66 %
     L1/CIN4            73         4150 bp    0.00 %
   LTR elements:     26568      7334529 bp    0.50 %
     BEL/Pao          3526      2100677 bp    0.14 %
     Ty1/Copia        4050       534118 bp    0.04 %
     Gypsy/DIRS1     15237      4326307 bp    0.29 %
       Retroviral        0            0 bp    0.00 %

DNA transposons     382289     77214231 bp    5.25 %
   hobo-Activator    16314      1202614 bp    0.08 %
   Tc1-IS630-Pogo   258740     65143621 bp    4.43 %
   En-Spm                0            0 bp    0.00 %
   MuDR-IS905            0            0 bp    0.00 %
   PiggyBac           2799       623450 bp    0.04 %
   Tourist/Harbinger   839       148989 bp    0.01 %
   Other (Mirage,      465        45742 bp    0.00 %
    P-element, Transib)

Rolling-circles          0            0 bp    0.00 %

Unclassified:       158272     25114896 bp    1.71 %

Total interspersed repeats:   174983102 bp   11.89 %

Small RNA:           10506       796882 bp    0.05 %

Satellites:           4162       211511 bp    0.01 %
Simple repeats:     580093     26953014 bp    1.83 %
Low complexity:     114453      6289474 bp    0.43 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element

The query species was assumed to be insecta       
RepeatMasker Combined Database: Dfam_Consensus-20170127, RepBase-20181026

run with rmblastn version 2.6.0+

</pre>

## 3.4 QC the assembly

We provide a script BUSCO_metrics_slurm.sh that runs assemblathon.pl to calculate assembly contiguity.

It also runs BUSCO (database: insecta_odb10)  to calculate genome completeness.

The gaps are reported in the assemblathon-generated file in these categories:

- Percentage of assembly in unscaffolded contigs

- contig %non-ACGTN

- Number of contig non-ACGTN nt











