# Background 

<pre>

Many genomics analyses require first establishing a reference genome. 
However, de novo genome assembly is a complicated and computationally intensive process with many assumptions hidden to the user. 
A popular assessment prior to genome assembly is genome profiling, 
where the k-mer frequencies within the sequencing reads are analyzed to efficiently estimate major genome characteristics 
such as genome size, heterozygosity, and repetitiveness. 

</pre>

Source: https://github.com/tbenavi1/genomescope2.0


Short reads from genomic DNA libraries are used to estimate genome size. In this project we used TellSeq libraries: https://www.universalsequencing.com/technology


# Dependencies


- jellyfish version 2.3.0 or higher http://www.genome.umd.edu/jellyfish.html#Release

- slurm workload manager https://slurm.schedmd.com/documentation.html

- GenomeScope2 server:  http://qb.cshl.edu/genomescope/genomescope2.0/


# Steps

- Identify and count kmers

- Calcualte histogram of kmer distribution

- Upload histogram file on the GenomeScope server to estimate genome size


# Edit the script

The script jellyfish_slurm.sh runs jellyfish on a linux cluster. 

Please edit the script with your information.

Notice that jellyfish expects the input to be uncompressed fastq files.


<pre>

jellyfish count -C -m 21 -s 5G -t $SLURM_NPROCS  -o ${PREFIX}_bothTrimmedReads.jf *.fastq

jellyfish histo -t $SLURM_NPROCS ${PREFIX}_bothTrimmedReads.jf > ${PREFIX}_bothTrimmedReads.histo 

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

