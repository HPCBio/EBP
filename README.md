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

2. purge duplicate contigs (see ./scripts/Raw_assembly/scaffolding )

3. scaffolding using TellSeq reads (see ./scripts/Raw_assembly/scaffolding )

4. scaffolding using Omni-C reads (see ./scripts/Raw_assembly/scaffolding )

5. gap filling (see ./scripts/gap_filling_and_masking/ )

6. masking repeats and low complexity regions (see ./scripts/gap_filling_and_masking/ )

7. functional annotation (see ./scripts/Annotation/ )

8. identify and annotate mitochondrial DNA (see ./scripts/mitofinder/ )

9. assess genome completeness w Merqury (see ./scripts/Merqury_completeness/ )

10. identify contaminants and artifacts in genome (see ./scripts/blobtools_contaminants_detection/ )



# Denovo genome assembly using CLR reads

Picture goes here

Text goes here

# Denovo genome assembly using short reads only

Not this repo