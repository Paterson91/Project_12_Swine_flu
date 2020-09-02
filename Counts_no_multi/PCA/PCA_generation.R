#!/usr/bin/env Rscript
library(factoextra) #ggplot2 based PCA graphics

setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #12 Swine Flu - Amy Thomas/Counts_no_multi/PCA")
args = commandArgs(trailingOnly = TRUE)

####################################################################
########################## Input formatting ########################
####################################################################

gene_counts_in = "All_normalized_counts_DESeq2.csv"
metadata_in = "metadata_omit.csv"
sampleFiles<- gene_counts_in

countData <- as.matrix(read.csv(gene_counts_in, header=TRUE,
                                  row.names = 1,
                                  as.is=TRUE))
head(countData)

metadata = read.csv(metadata_in, header = TRUE)

cat("\n\nMetadata input: \n\n")
metadata


#Filter on significant
genes = read.csv("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_P0_vs_P1.csv", header = T)
head(genes)
nrow(genes)

sig_genes = genes[which(genes$padj<0.05),]
head(sig_genes)
nrow(sig_genes)
sig_gene_names = as.character(sig_genes$X)

order1_all = as.character(metadata$Sample)
iso_table_all = countData[, match(order1_all, colnames(countData))]
head(iso_table_all)
colnames(iso_table_all)
metadata$Sample
head(iso_table_all)
nrow(iso_table_all)

iso_table_sig = iso_table_all[sig_gene_names,]
head(iso_table_sig)
nrow(iso_table_sig)

#Transpose to make sample names rownames
df_alex_t = t(iso_table_sig)
ncol(df_alex_t)
nrow(df_alex_t)

# # Filter on output
# keep_lrt <- colSums(df_alex_t) >= 1
# keep = df_alex_t[,keep_lrt]

res.pca_norm_out <- prcomp(df_alex_t, scale. = T)

fviz_eig(res.pca_norm_out)
# project.pca.proportionvariances <- ((res.pca_norm_out$sdev^2) / (sum(res.pca_norm_out$sdev^2)))*100
# barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(res.pca_norm_out$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))


#DESeq generated counts - VSD transformed
# counts = counts(dds_lrt)
# vsd = varianceStabilizingTransformation(counts)
# vsd2 = varianceStabilizingTransformation(dds_lrt)
# res.pca_vsd <- prcomp(vsd, scale. = F)
# fviz_eig(res.pca_vsd)
# plotPCA(vsd2, ntop = 100000, returnData = F, intgroup = c("phase2_cat","group"))


#visualise eigenvalues (scree plot)
scree = fviz_eig(res.pca_norm_out)
ggsave("PCA_scree.tiff", plot = scree)
dev.off()


#PCA for Phase of infection
groups_phase <- as.factor(metadata$phase2_cat)
mycolors_phase <- c("blue", "yellow", "tan1", "tomato3", "red4", "darkgray")

#Inclusion of ellipses and annotations
pca_ind_phase = fviz_pca_ind(res.pca_norm_out,
             col.ind = groups_phase, # Color by the quality of representation
             repel = T, # Avoid text overlapping
             addEllipses = T, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Infection Phase",
             palette = mycolors_phase
)
pca_ind_phase
ggsave("PCA_ind_phase.tiff", plot = pca_ind_phase)
dev.off()

#PCA for Status of infection
groups_inf <- as.factor(metadata$infection_status)
myColors <- c("blue", "red", "darkgray")

#Inclusion of ellipses and annotations
pca_ind_inf = fviz_pca_ind(res.pca_norm_out,
             col.ind = groups_inf, # Color by the quality of representation
             repel = T, # Avoid text overlapping
             addEllipses = T, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Infection Status",
             palette = myColors
)
pca_ind_inf
ggsave("PCA_ind_inf.tiff", plot = pca_ind_inf)
dev.off()


pca_var = fviz_pca_var(res.pca_norm_out,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = F) # Avoid text overlapping
pca_var
ggsave("PCA_var.tiff", plot = pca_var)

#res.pca_norm_out <- prcomp(df_alex_t, scale. = F)

#eigenvalues
eig.val <- get_eigenvalue(res.pca_norm_out)

# Results for Variables
res.var <- get_pca_var(res.pca_norm_out)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

