# EBP
EDIT ME!

Text goes here

# The Workflow

<p>
<img align="left" src="./docs/EBP_Workflow_1.png" />

</br></br></br>
</p>



# Denovo genome assembly using HiFi reads

These are the steps:


1. generate raw assembly with hifiasm  (see ./scripts/Raw_assembly/ )

2. purge duplicate contigs (see ./scripts/scaffolding )

3. scaffolding using TellSeq reads (see ./scripts/scaffolding )

4. scaffolding using Omni-C reads (see ./scripts/scaffolding )

5. gap filling (see ./scripts/gap_filling_and_masking/ )

6. masking repeats and low complexity regions (see ./scripts/gap_filling_and_masking/ )

7. functional annotation (see ./scripts/Annotation/ )

8. identify and annotate mitochondrial DNA (see ./scripts/mitofinder/ )

9. assess genome completeness w Merqury (see ./scripts/Merqury_completeness/ )

10. identify contaminants and artifacts in genome (see ./scripts/blobtools_contaminants_detection/ )



# Denovo genome assembly using CLR reads

1. generate raw assembly with Redbean  (see ./scripts/Raw_assembly/ )

2. base-correct assembly with Arrow (see ./scripts/Arrow_polish )

3. purge duplicate contigs (see ./scripts/scaffolding )

4. Pilon polishing (see ./scripts/pilon_polishing )

5. scaffolding using TellSeq reads (see ./scripts/scaffolding )

6. scaffolding using Omni-C reads (see ./scripts/scaffolding )

7. gap filling (see ./scripts/gap_filling_and_masking/ )

8. masking repeats and low complexity regions (see ./scripts/gap_filling_and_masking/ )

9. functional annotation (see ./scripts/Annotation/ )

10. FreeBayes polishing (see ./scripts/FreeBayes_polishing )

11. identify and annotate mitochondrial DNA (see ./scripts/mitofinder/ )

12. assess genome completeness w Merqury (see ./scripts/Merqury_completeness/ )

13. identify contaminants and artifacts in genome (see ./scripts/blobtools_contaminants_detection/ )



# Denovo genome assembly using short reads only

Not this repo