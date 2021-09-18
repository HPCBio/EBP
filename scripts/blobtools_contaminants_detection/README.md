# README

<pre>

Modular command-line solution for visualisation, quality control and taxonomic partitioning of genome datasets 

</pre>

blobtools page:  https://github.com/DRL/blobtools

The tutorial is here: https://blobtools.readme.io/docs/my-first-blobplot

We will use blotools to visualize and QC the assembly.

## Steps

- Run Diamond to generate the hits file ( blast.out )

- Run bwa with the TellSeq files to generate the coverage file ( mapping.bam )

- Create a blobtools database with genome blast.out mapping.bam

- Generate blobtools plots with the plot module

- Generate blotools taxonomic tables with the vew  module

The script blobtools.slurm.sh performs all these steps automatically


# Outputs


- Example of output file for taxonomic classification at the phylum level generated with blobtools view

<pre>
## 1.1.1
## assembly	: /home/groups/earthbiogenome/results/20210607_EB31_blobtools/EB31_Sealer_scaffold.fa
## coverage	 bam0 - /home/groups/earthbiogenome/results/20210607_EB31_blobtools/EB31_TellSeq_10x_format_trimmed_sorted.bam
## taxonomy	 tax0 - /home/groups/earthbiogenome/results/20210607_EB31_blobtools/EB31_matches.daa.tagc
## nodesDB	: None
## taxrule	: bestsum
## min_score	: 0.0
## min_diff	: 0.0
## tax_collision_random	: False
##
# name	length	GC	N	bam0	phylum.t.6%s	phylum.s.7%s	phylum.c.8
scaffold1	231675529	0.357	14100	67.6016	Arthropoda	50095.0	3
scaffold2	72918906	0.3524	5700	64.9537	Arthropoda	43550.0	2
scaffold3	67397341	0.3511	11000	62.4446	Arthropoda	60970.0	1
scaffold4	66490383	0.3527	10700	60.9654	Arthropoda	38928.0	2
scaffold5	61989192	0.3531	5700	63.8268	Arthropoda	40262.0	3
scaffold6	61707333	0.3589	5400	62.8597	Arthropoda	78798.0	0
scaffold7	60319155	0.3629	7300	39.6533	Arthropoda	38052.0	2
scaffold8	51257439	0.3497	1700	63.6249	Arthropoda	51718.0	2
scaffold9	50414819	0.3481	7200	62.8704	Arthropoda	44286.0	2
</pre>


- Example of the plot produced for taxonomic assignment at the phylum level too with blobtools plot

<p>
<img align="left" src="/docs/EB31_Sealer_blobplot.blobDB.json.bestsum.phylum.p20.span.100.blobplot.bam0.png" />


</br></br></br>
</p>

<p>
</br></br></br>
</p>

.
.
- Example of the plot produced for read coverage with blobtools plot


<p>
<img align="left" src="/docs/EB31_Sealer_blobplot.blobDB.json.bestsum.superkingdom.p20.span.100.blobplot.read_cov.bam0.png" />
</br></br></br>
</p>

<p>
</br></br></br>
</p>

.
.


# Using blobtools to remove contaminants


The taxonomic classification tables can be used to identify which contigs should be removed from the assembly.

In this example we are looking at the taxonomic classification results of the de-novo genome assembly of an insect;
therefore we should remove contigs assigned to Bacteria and Viruses.

<pre>
scaffold434	46611	0.3227	0	40.5496	Eukaryota	16229.0	0
scaffold435	46521	0.3997	0	39.2916	Eukaryota	3789.0	0
scaffold436	45743	0.362	0	60.5585	no-hit	0.0	0
scaffold437	45547	0.3634	0	47.0578	no-hit	0.0	0
scaffold438	45086	0.336	0	47.1976	Eukaryota	18958.0	1
scaffold439	45076	0.3578	0	67.8022	Eukaryota	7270.0	0
scaffold440	44902	0.3785	0	82.2102	Viruses	3161.0	1
scaffold441	44789	0.3427	0	40.7488	Eukaryota	5523.0	0
scaffold442	44593	0.3663	0	43.4363	no-hit	0.0	0
scaffold443	44380	0.3546	0	61.6374	Eukaryota	18805.0	1
scaffold444	44101	0.3064	0	138.553	Eukaryota	1747.0	0
scaffold445	44043	0.3614	0	49.0197	no-hit	0.0	0
scaffold446	42860	0.3878	0	66.4049	Eukaryota	19479.0	0
scaffold447	42463	0.3708	0	51.3088	Eukaryota	21547.0	0
scaffold448	41918	0.366	0	46.4973	no-hit	0.0	0
scaffold449	41890	0.3319	0	55.1293	Eukaryota	2491.0	0
scaffold450	41464	0.3685	0	76.2572	Viruses	4869.0	0
scaffold451	41195	0.3348	0	40.1883	Bacteria	2986.0	1
</pre>
