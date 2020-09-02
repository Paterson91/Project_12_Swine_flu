setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #12 Swine Flu - Amy Thomas/Counts_no_multi/GO_analysis")
library(dplyr)
library(tidyr)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("org.Ss.eg.db")
#BiocManager::install("gtools")
library(gtools)
library(clusterProfiler)
library("org.Ss.eg.db")

#Data input

input = read.csv("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_P0_vs_P1.csv", header = TRUE)
head(input)
nrow(input)
Simple_in = input[,c("X","log2FoldChange","padj")]
colnames(Simple_in) = c("Gene","L2FC","padj")
head(Simple_in)
#Convert to fold change
Simple_in$FC =(2^Simple_in$L2FC)
head(Simple_in)

require("biomaRt")

#Further Annotation
listMarts()
mart <-useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(mart)
mart <-useDataset("sscrofa_gene_ensembl", mart)
filters = listFilters(mart)
filters

annotLookup <- getBM(
  mart=mart,
  attributes=c("hgnc_symbol","hgnc_id","ensembl_gene_id","entrezgene_id","external_gene_name"),
  filter="external_gene_name",
  values=Simple_in$Gene,
  uniqueRows=TRUE)
head(annotLookup)

all_annots = left_join(Simple_in,annotLookup, by = c("Gene"="external_gene_name"))
head(all_annots)

### GO Classification
#Set levels
#Set choice of; CC - Cellular Component, MF - Molecular Function, BP - Biological Process

#filter P Value
Psig_simple = filter(all_annots, padj<0.05)
head(Psig_simple)

ont_selection = "BP"
level_selection = 2

ggo1 <- groupGO(gene     = as.character(Psig_simple$entrezgene_id),
                keyType = "ENTREZID",
                OrgDb    = org.Ss.eg.db,
                ont      = ont_selection,
                level    = level_selection,
                readable = TRUE)

barplot(ggo1, drop=TRUE, showCategory=20, border = c(10,10))
dev.copy(png, width = 1000,height = 1000, paste0(ont_selection, "_lvl_", level_selection, "_Go.png"))
dev.off()
head(ggo1)


#Isolate Genes belonging to Immune Process GO:0006955 

GO_of_int = "GO:0002376"

immune_total = ggo1[GO_of_int,]
immune_sep =immune_total$geneID
immune_subset = as.data.frame(strsplit(immune_sep, split = "/"))
colnames(immune_subset) = GO_of_int
head(immune_subset)
nrow(immune_subset)



