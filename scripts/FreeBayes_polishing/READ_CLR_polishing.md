# 1. Background

<pre>
Long-read, single-molecule technologies, such as those produced by Pacific Biosciences8 (PacBio) and Oxford Nanopore9 (ONT), 
have the potential to sequence DNA molecules with lengths in the tens of thousands or hundreds of thousands of bases, 
enabling researchers to assemble large and complex repeats. 
However, both of these technologies have high per-read error rates (on the order of 5–15%), 
which has resulted in the development of ‘correction’ algorithms. 
These algorithms attempt to use consensus base-calls, raw signal data and/or shorter, 
more accurate reads to correct long-read assemblies. 
Examples include Quiver and Arrow for PacBio, Nanopolish for ONT, and Pilon.
</pre>

Source:

Errors in long-read assemblies can critically affect protein prediction
Mick Watson & Amanda Warr 
Nature Biotechnology volume 37, pages 124–126 (2019)
DOI  https://doi.org/10.1038/s41587-018-0004-z


Arrow page:
https://github.com/PacificBiosciences/GenomicConsensus/blob/develop/doc/HowTo.rst

Pilon page:
https://github.com/broadinstitute/pilon/wiki

FreeBayes:
https://github.com/VGP/vgp-assembly/tree/master/pipeline/freebayes-polish
https://github.com/freebayes/freebayes

# 2. Dependencies.

- smrtlink version 9.0.0.92188 or higher

- samtools version 1.12 or higher

- BWA version 0.7.17

- Java 1.8.0_152 

- pilon version 1.23

- Python version 3.7.2

- FreeBayes version 1.3.4

- VCFtools version 0.1.16

- BCFtools version 1.12

- Merfin version 20210507

- meryl version 1.3

- Jellyfish version 2.3.0

- slurm workload manager https://slurm.schedmd.com/documentation.html



# 3. Arrow polishing

Run arrow polishing after the raw assembly is generated with CLR reads.

The arrow polishing tool uses the CLR reads to correct base-errors in the assembly.

The first part aligns the long reads to the raw assembly and generates a sorted.bam file.
The execution of this script could take a few days.

The second part splits the sorted.bam file into chunks; one chunk per contig
and then launching the arrow polishing command on each chunk.
The execution of this script could take a long time as in weeks;
depending on the available compute resources and how fragmented the raw assembly is.

The third part collects all results, concatenates them and generates the arrow_polished genome.
It also runs assemblathon and BUSCO to calculate completeness and contiguity of the new assembly.
The execution of this script could take one to two days.


# 4. Pilon polishing

Run Pilon polishing after the purg_dups step.

Pilon uses  short reads, such as TellSeq, that are highly accurate and almost error-free to improve the assembly.

<pre>

Pilon requires as input a FASTA file of the genome along with one or more BAM files of reads aligned to the input FASTA file. 
Pilon uses read alignment analysis to identify inconsistencies between the input genome and the evidence in the reads. 
It then attempts to make improvements to the input genome, including:

- Single base differences
- Small indels
- Larger indel or block substitution events
- Gap filling
- Identification of local misassemblies, including optional opening of new gaps

Pilon then outputs a FASTA file containing an improved representation of the genome from the read data 
and an optional VCF file detailing variation seen between the read data and the input genome.

</pre>

See: https://github.com/broadinstitute/pilon/wiki


# 5. FreeBayes polishing

Run FreeBayes after gap-closing and masking the CLR assembly

FreeBayes also uses  short reads, such as TellSeq, that are highly accurate and almost error-free to improve the assembly.


See    https://github.com/freebayes/freebayes

There are a couple different ways to actually call the consensus from the called variants.

- Method described in VGP pipeline.

        https://github.com/VGP/vgp-assembly/tree/master/pipeline/freebayes-polish

-  Method that incorporate merfin and jellyfish, suggested by Andrew Severin at Iowa State.

        https://github.com/arangrhie/merfin

We tried both consensus methods as I've seen differences with SCN testing.







