#!/usr/bin/env Rscript
library("tidyverse")
library(dplyr)
library(readr)
setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #12 Swine Flu - Amy Thomas/Counts_no_multi/DESeq2_output/replicate + phase2_cat/Output_Figures/")
getwd()
args = commandArgs(trailingOnly = TRUE)


####################################################################
########################## Input formatting ########################
####################################################################

df_tibble <- list.files(path = "../DESeq2_LRT_OutputTables/", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_cols()

df = as.data.frame(df_tibble,row.names = "X1")

df_sig = filter(df, (padj<0.05))
head(df_sig)

df_simple = df_sig %>%
              dplyr::select(X1, log2FoldChange, log2FoldChange1,log2FoldChange2,log2FoldChange3,log2FoldChange4)
head(df_simple)

rownames(df_simple) = df_simple$X1
df_simple$X1 = NULL

colnames(df_simple) = c("P1","P2","P3","P4","P5")
head(df_simple)

mat_data <- data.matrix(df_simple[,1:ncol(df_simple)])
head(mat_data)

#Immune genes
immmune_sig_lookup = read.csv("../GO_immune_sig.csv", header =T)
head(immmune_sig_lookup)
immune_subset_list = as.character(immmune_sig_lookup$GO.0002376)
ordered_input_immune = mat_data[match(immune_subset_list, rownames(mat_data)),]
head(ordered_input_immune)
nrow(ordered_input_immune)
rownames(ordered_input_immune) = as.character(rownames(ordered_input_immune))
ordered_input_immune
#Note row 15 shows NA for some reason... FIX THIS
ordered_input_immune = na.omit(ordered_input_immune)

####################################################################
############################ Heatmap ###############################
####################################################################
library("RColorBrewer")
library(viridis)
library("gplots")
library(rafalib)

# 100 top expressed genes with heatmap.2
#select <- order(rowMeans(counts(dds_lrt,normalized=T)),decreasing=T)[1:n]
#cols <- palette(magma(256))[as.fumeric(as.character(ColourStrategy))]

heatmap.2(ordered_input_immune, 
          col=magma(256),
          Colv=F,
          dendrogram="row",
          scale="row",
          key=T,
          keysize=1,
          symkey=T,
          lhei = c(1.5,5),
          density.info="none",
          trace="none",
          margins=c(6,10),
          cexCol=1.1,
          labRow=rownames(ordered_input_immune),
          main="\nSignificantly regulated Immune Genes \n- Log2 Fold Change")
dev.copy(png, width = 1000,height = 1000, "Significantly_reg_immune_l2fc-Heatmap_DESeq2.png")
dev.off()

# TAKE 2 - Trying with normalised counts...
####################################################################
########################## Input formatting ########################
####################################################################

normalised_reads = read.csv("../All_normalized_counts_DESeq2.csv", header = T, row.names = 1)
head(normalised_reads)

metadata = read.csv("../../../metadata_omit.csv", header = T)
head(metadata)

phase_order = metadata[order(metadata$phase2_cat),]
order1_all = as.character(phase_order$Sample)

ordered_input = normalised_reads[, match(order1_all, colnames(normalised_reads))]
head(ordered_input)

ColourStrategy = phase_order$phase2_cat

#n = 100


ordered_input_immune_reads = ordered_input[match(immune_subset_list, rownames(ordered_input)),]
ordered_input_immune_reads
nrow(ordered_input_immune_reads)

cols <- palette(magma(256))[as.fumeric(as.character(ColourStrategy))]
mat_data2 <- data.matrix(ordered_input_immune_reads[,1:ncol(ordered_input_immune_reads)])

heatmap.2(mat_data2, 
          col=magma(256),
          Colv=F,
          dendrogram="row",
          scale="row",
          key=T,
          keysize=1,
          symkey=T,
          lhei = c(1.5,5),
          density.info="none",
          labCol = ColourStrategy,
          ColSideColors=cols,
          trace="none",
          margins=c(6,7),
          cexCol=1.1,
          labRow=rownames(mat_data2),
          main="\n\nSignificantly regulated Immune Genes \n- Normalised Reads")
dev.copy(png, width = 1000,height = 1000, "Significantly_reg_immune_reads-Heatmap_DESeq2.png")
dev.off()