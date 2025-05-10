# Master_thesis_analysis
Showing my custom code as part of my master thesis in bioinformatics, analysing mitochondrial data from a mtDNA-mutator mouse cell culture. The raw data is not yet published and as such will not be part of this repository, rendering the code unfortunately untestable.

The data used in the analysis uploaded here comes from 3 sequencing experiments: Short read DNA sequencing with Illumina, long read DNA sequencing with PacBio and short read RNA sequencing with Illumina. Both DNA datasets were first amplified by PCR into 3 amplicons, both to enrich mtDNA from a low yield as well as to avoid contamination from nuclear NUMTs. The short read DNA data was processed and variant called by the time I joined the project.

My pipeline for the RNA data used Trimmomatic, BWA MEM and Samtools mpileup before I applied my custom script create_RNA_nucleotide_table.sh included in this repository. The long read DNA data had been processed with the ccs pipeline by company that did the sequencing, so for this I only used the python package pysam, as can be seen in the jupyter notebook Single_nucleotide_analysis.ipynb.
