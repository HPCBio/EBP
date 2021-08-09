# 1. Genome assessment


<pre>
Contiguity is related to the size and number of contigs. 
If a goal of assembly is to reflect the contiguity of the genome in vivo, 
then the assembly process seeks to maximize the size of the contigs 
and to minimize the number of contigs to reflect the true number and size of the chromosomes in the organism. 
Contiguity errors may arise due to assembler parameters that allow unrelated contigs to be joined 
or that prevent related contigs from being joined.

Completeness is determined by the content of contigs, especially with regard to gene content. 
A contiguous genome that contains no gene content is not useful for downstream analysis. 
Completeness errors can arise in sequencing (important genes may not be sequenced) 
or they may arise in the assembly process (genes may end up in discarded contigs).

Correctness is concerned with the ordering and location of contigs. 
A correct genome assembly has the same order as the true genome. 
If contigs are incorrect, they may have inversions, relocations, or translocations with respect to the true genome.
</pre>

Source: 

Thrash, A., Hoffmann, F. & Perkins, A. 
Toward a more holistic method of genome assembly assessment. 
BMC Bioinformatics 21, 249 (2020). 
https://doi.org/10.1186/s12859-020-3382-4


## 1.1 Assembly contiguity 

We measure assembly contiguity with asemblathon. Genome Res. 2011 Dec; 21(12): 2224–2241.  doi: 10.1101/gr.126599.111

The assemblathon script calculates several genoeme statistics which  include:

- overall assembly size 
- assembly contiguity (N50, NG50, NA50, or NGA50) 
- the number of contigs, contig length and contig mean length

In computational biology, N50 is a widely used metric for assessing an assembly’s contiguity, 
which is defined by the length of the shortest contig for which longer and equal-length contigs cover at least 50% of the assembly. 
NG50 resembles N50 except for the metric, which relates to the genome size rather than the assembly size. 
NA50 and NGA50 are analogous to N50 and NG50 where the contigs are replaced by blocks aligned to the reference


## 1.2 Assembly completeness

We measure assembly completeness with BUSCO. See https://busco.ezlab.org/busco_userguide.html

<pre>

which searches for benchmarking universal single-copy orthologs (BUSCOs) in an assembly. 
BUSCO measures these orthologs by counting complete single-copy BUSCOs, fragmented BUSCOs, missing BUSCOs, and duplicate BUSCOs. 
The authors of BUSCO note that duplicate BUSCOs may represent misassemblies 
where a heterozygous allele failed to collapse into the assembly properly and was retained as a contig 

</pre>

Source: 

Thrash, A., Hoffmann, F. & Perkins, A. 
Toward a more holistic method of genome assembly assessment. 
BMC Bioinformatics 21, 249 (2020). 
https://doi.org/10.1186/s12859-020-3382-4

We provide the script BUSCO_metrics_slurm.sh to assess assembly contiguity and completeness


# 2. Dependencies

- BUSCO 4.1.4 or higher  https://busco.ezlab.org/busco_userguide.html

- BUSCO insecta-odb10 https://busco.ezlab.org/busco_v4_data.html

- augustus 3.3.3  or higher

- Perl 5.28.1 or higher

- slurm workload manager https://slurm.schedmd.com/documentation.html

# 3. The script

We provide a script BUSCO_metrics_slurm.sh that runs assemblathon.pl and then BUSCO.

## 3.1 Setup

Before running the script, you need to do the following steps

- If the BUSCO toolkit was installed in a path that is write-protected; then you need to copy the config folder to you home folder

- If the Augustus tool was installed in a path that is write-protected; then you need to copy the config.ini file to you home folder

- If the BUSCO database is not found in the default path (BUSCO/data/lineages/insecta_odb10); then you need to specify the full path in the script

- Copy the assemblathon folder (provided with this repo) to your home folder

## 3.2 Run the script

<pre>
sbatch BUSCO_metrics_slurm.sh
</pre>

The script can take anywhere between 2 hours - two days

## 3.3 The outputs.

The assemblathon output will be a single file called yourgenomename.fasta.stats
The top portion looks like this:


<pre>


                                         Number of scaffolds        789
                                     Total size of scaffolds 1647789986
                                            Longest scaffold  102939053
                                           Shortest scaffold      15367
                                 Number of scaffolds > 1K nt        789 100.0%
                                Number of scaffolds > 10K nt        789 100.0%
                               Number of scaffolds > 100K nt        295  37.4%
                                 Number of scaffolds > 1M nt        149  18.9%
                                Number of scaffolds > 10M nt         34   4.3%
                                          Mean scaffold size    2088454
                                        Median scaffold size      48120
                                         N50 scaffold length   39672384
                                          L50 scaffold count         12

...

</pre>

The BUSCO results are written to a folder BUSCO_yourgenome...

The summary file will look like this:

</pre>
# BUSCO version is: 4.1.4 
# The lineage dataset is:  insecta_odb10 (Creation date: 2020-09-10, number of species: 75, number of BUSCOs: 1367)
# Summarized benchmarking in BUSCO notation for file EB17_hifiasm_scaffolded_shortNames.fa
# BUSCO was run in mode: genome

    ***** Results: *****

    C:97.0%[S:91.1%,D:5.9%],F:1.1%,M:1.9%,n:1367       
    1327    Complete BUSCOs (C)               
    1246    Complete and single-copy BUSCOs (S)       
    81    Complete and duplicated BUSCOs (D)       
    15    Fragmented BUSCOs (F)               
    25    Missing BUSCOs (M)               
    1367    Total BUSCO groups searched   
    
</pre>




