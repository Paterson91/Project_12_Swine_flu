# Config file for running DESeq analysis

#Inputs
gene_counts_in="Counts_stranded_gtf_multi.gene_tidied"
metadata_in="metadata_omit.csv"
#design_in="~phase2_cat + replicate"
design_in="~replicate + phase2_cat"






reference_in="P0"

~group + phase2_cat
~phase2_cat + group


#Outputs


