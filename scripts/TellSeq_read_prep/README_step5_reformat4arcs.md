# Background

ARCS is a scaffolding tool. The ARCS page: https://github.com/bcgsc/arcs

# ARCS reformatting

ARCS expects one file of interleaved paired reads. 

It also expects the barcode to appear as a BX:Z: tag at the end of the header line of each read like this

<pre>
@A00835:176:HTCNTDRXX:2:2101:1136:1016 BX:Z:GTCAGTGCGGTTAGGATA
TNGTGGATTTAACAAATTTAATTATGCAATATTTCAGTCGGCGACCGCACATTTTTATAAAATAAAGAGATTGCAATTGAGAAGATGGCAAAGGGGAAAGGGAAGAAGTAAATGTGTCGAAAGAATGGGCGGAATTTAAATAAAATGTAA
+
F#FFFFF:FF:FFF:,FFF:FFFFFFFFF:FF:FFFFFF:F::FFFFFFFFFFF,FFFFFF:,:FF::F:,FF,FFFFFF,:FFFF:F:::,,F:FF,,:,FF,,:::F,,,:,F,,,F,,FF,::F,:,,,,,,,F:,,F,,:,F,::,
@A00835:176:HTCNTDRXX:2:2101:1136:1016 BX:Z:GTCAGTGCGGTTAGGATA
GCCGTGTTGATGTGGAAAAACGATAGGAGCGATTAAGATTTGAGAATTATTACTCTATATTATTTTTTATTTAAATTCGTTGGATTATTTCGTGTTTTTTTTTTGATAACTTTGAAATTTGGAATAGTATAAAATGAACTCTATTATTTG
+
,,,:F,:F,:F,FF:F::::,F,FF:,F,,:F,,,,,FFFF,F::F,F,F,:,:,,:,,FF:,,FF:,,:,F:,,F,,,:,,,,FF,,F,,:F,,,,F,F,,FF,,:,,,:,,,,,,F,,,,,,::,F,F,F:,:F,,,F,F,,,::,,,
</pre>

# In-house script interleave_us10x.pl

We wrote a perl script called interleave_us10x.pl that appends the tag to the header of each read and generates an interleaved file.

Run the script like this:

<pre>

mkdir /home/a-m/grendon/EB14_nobc/ARCS_inputs

cd /home/a-m/grendon/EB14_nobc/ARCS_inputs

/home/a-m/grendon/TELL-Seq/conversion_tool/interleave_us10x.pl \
-r1 ../Full/EB14_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq \
-r2 ../Full/EB14_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq \
-index ../Full/EB14_nobc_I1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq \
-out EB14_nobc_T500_4ARCS.fastq

pigz -c -p12 EB14_nobc_T500_4ARCS.fastq > EB14_nobc_T500_4ARCS.fastq.gz

</pre>
