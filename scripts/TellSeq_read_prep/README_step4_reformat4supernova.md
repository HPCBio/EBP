# Background

Supernova expects two input files only, R1 and R2. 

It also expects to find the barcodes at the start of the sequence of each read in R1. 

The usx vendors provided a script that reformats the TellSeq files as 10x-like files.

The usx vendors also provide a database of barcodes.

However, the lengths of these barcodes is different than 10x barcodes.

Moreover, there are more barcodes in this database than in the 10x database.

Therefore; the regular supernova toolkit has to be modified with this custome database in order to work properly.

# Install the conversion toolkit

Copy the toolkit provided by the vendor and uncompress it.  No further installation is needed

<pre>

$ mkdir /home/a-m/grendon/TELL-Seq/conversion_tool/
$ cd /home/a-m/grendon/TELL-Seq/conversion_tool/
$ unzip conversion.zip
$ ls 
ust10x
4M-with-alts-february-2016.txt

</pre>

4M-with-alts-february-2016.txt  is the database of barcodes

ust10x is the reformatting script

# run the conversion script

<pre>

## input files cannot be compressed

$ cd /home/a-m/grendon/EB14_nobc/Full/

$ zcat -d EB14_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz > EB14_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq &

$ zcat -d EB14_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz > EB14_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq &

$ zcat -d EB14_nobc_I1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz > EB14_nobc_I1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq &

mkdir /home/a-m/grendon/EB14_nobc/supernova_inputs

cd /home/a-m/grendon/EB14_nobc/supernova_inputs

/home/a-m/grendon/TELL-Seq/conversion_tool/ust10x -sz 16000000 \
-i1 ../Full/EB14_nobc_I1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq \
-r1 ../Full/EB14_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq \
-r2 ../Full/EB14_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq \
-wl /home/a-m/grendon/TELL-Seq/conversion_tool/4M-with-alts-february-2016.txt 

</pre>

Where:

- -sz  is the expected size of the database
- -i1  is the file with barcodes
- -r1  is the file with forward reads
- -r2  is the file with reverse reads
- -w1 is the filename of the database

# Outputs

Two files that have fixed filenames as shown below:

<pre>

$ ls
R1_sl.fastq.gz.4tenx.fastq
R2_sl.fastq.gz.4tenx.fastq

</pre>


# Rename and compress files

Supernova expects filenames to adhere to a certain pattern. 
It also expects the files to be compressed

<pre>
        <sampleID>_S1_L001_R1_001.fastq.gz
</pre>

Unfortunately the script does not perform these two tasks. You have to do them manually.

<pre>

pigz -c -p12  R1_sl.fastq.gz.4tenx.fastq > EB14_nobc_4tenx_S1_L001_R1_001.fastq.gz &

pigz -c -p12  R2_sl.fastq.gz.4tenx.fastq > EB14_nobc_4tenx_S1_L001_R2_001.fastq.gz &

</pre>