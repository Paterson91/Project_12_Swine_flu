setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #12 Swine Flu - Amy Thomas/Counts_no_multi/KEGG_pathways")
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

input_files = list.files(pattern="*.csv", path = "../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/", full.names = T)
input_files
all_input = lapply(input_files, read.csv,header = T)
names(all_input) <- gsub(".csv","",
                         list.files(pattern="*.csv", path = "../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/", full.names = F),
                       fixed = TRUE)


#dfList = list(P0_vs_P1_Simple_in,P0_vs_P2_Simple_in,P0_vs_P3_Simple_in,P0_vs_P4_Simple_in,P0_vs_P5_Simple_in)

simplified = lapply(all_input, function(x) {
  x[,c("X","log2FoldChange","padj")]
})
new_col_names = c("external_gene_name","L2FC","padj")
simplified = lapply(simplified, setNames, nm = new_col_names)
#Convert to fold change

simplified = lapply(simplified, function(x) {
  x$FC =(2^x$L2FC);return(x)
})

require("biomaRt")

#Further Annotation
listMarts()
mart <-useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(mart)
mart <-useDataset("sscrofa_gene_ensembl", mart)
filters = listFilters(mart)
filters
annots = lapply(simplified, function(x) {
                  y = getBM(
                  mart=mart,
                  attributes=c("external_gene_name","hgnc_symbol","hgnc_id","ensembl_gene_id","entrezgene_id"),
                  filter="external_gene_name",
                  values=x$external_gene_name,
                  uniqueRows=TRUE)
                })

#Make new names for annot return
# newnames = lapply(names(all_input), function(x){
#   paste0(x, "_annots")
# })
# names(annots) = newnames

#Merge datasets with annotations

merged = Map(merge ,simplified, annots, by.x = "external_gene_name", by.y = "external_gene_name", simplify = FALSE)

head(merged[[3]])

# Plot specific KEGG pathways (with fold change) 
## pathway.id : KEGG pathway identifier

pathways = read.csv("KEGG_immune_ssc.csv", colClasses = c("character","character"), header = T)
pathways

# ## feature 1: numeric vector
# #geneList = d[,2]
# geneList = lapply(merged, function(x){
#   x$L2FC
# })

geneList = lapply(merged, function(x){
  x[c("entrezgene_id","L2FC")]
})

merged_no_na = lapply(geneList, function(x) na.omit(x))

## feature 3: decreasing orde
geneList_sorted = lapply(merged_no_na, function(x){
  x[order(x$L2FC, decreasing = T),]
})
head(geneList_sorted)

# names(geneList) = as.character(merged_no_na$DESeq2_P0_vs_P1)
# 
# lapply(merged_no_na,function(x){
#   names(geneList) = as.character(x$entrezgene_id)
# })
# 
# names(geneList) = lapply(merged_no_na,function(x){
#   as.character(x$entrezgene_id)
# })
# 
# ## feature 2: named vector
# names(geneList) = as.character(d$entrezgene_id)
# 
# ## feature 3: decreasing orde
# geneList_sorted = sort(geneList, decreasing = TRUE)
# #de <- names(geneList)[abs(geneList) > 2]
# head(geneList)
# barplot(sort(geneList, decreasing = T))


library(pathview)
getwd()

P0_vs_P1 = merged_no_na$DESeq2_P0_vs_P1$L2FC
names(P0_vs_P1) = as.character(merged_no_na$DESeq2_P0_vs_P1$entrezgene_id)
head(P0_vs_P1)
dir.create("P0_vs_P1")
setwd("P0_vs_P1/")
for (i in pathways$Pathway){
  print(as.character(i))
  pathview(gene.data = P0_vs_P1, 
           pathway.id = as.character(i), 
           species = "ssc")}

P0_vs_P2 = merged_no_na$DESeq2_P0_vs_P2$L2FC
names(P0_vs_P2) = as.character(merged_no_na$DESeq2_P0_vs_P2$entrezgene_id)
head(P0_vs_P2)
setwd("../")
dir.create("P0_vs_P2")
setwd("P0_vs_P2/")
for (i in pathways$Pathway){
  print(as.character(i))
  pathview(gene.data = P0_vs_P2, 
           pathway.id = as.character(i), 
           species = "ssc")}

P0_vs_P3 = merged_no_na$DESeq2_P0_vs_P3$L2FC
names(P0_vs_P3) = as.character(merged_no_na$DESeq2_P0_vs_P3$entrezgene_id)
head(P0_vs_P3)
setwd("../")
dir.create("P0_vs_P3")
setwd("P0_vs_P3/")
for (i in pathways$Pathway){
  print(as.character(i))
  pathview(gene.data = P0_vs_P3, 
           pathway.id = as.character(i), 
           species = "ssc")}

P0_vs_P4 = merged_no_na$DESeq2_P0_vs_P4$L2FC
names(P0_vs_P4) = as.character(merged_no_na$DESeq2_P0_vs_P4$entrezgene_id)
head(P0_vs_P4)
setwd("../")
dir.create("P0_vs_P4")
setwd("P0_vs_P4/")
for (i in pathways$Pathway){
  print(as.character(i))
  pathview(gene.data = P0_vs_P4, 
           pathway.id = as.character(i), 
           species = "ssc")}

P0_vs_P5 = merged_no_na$DESeq2_P0_vs_P5$L2FC
names(P0_vs_P5) = as.character(merged_no_na$DESeq2_P0_vs_P5$entrezgene_id)
head(P0_vs_P5)
setwd("../")
dir.create("P0_vs_P5")
setwd("P0_vs_P5/")
for (i in pathways$Pathway){
  print(as.character(i))
  pathview(gene.data = P0_vs_P5, 
           pathway.id = as.character(i), 
           species = "ssc")}

# names(merged_no_na)
# 
# for (z in names(merged_no_na)){
#   print(z)
#   nam = paste(z,"_input",sep="")
#   assign(nam, merged_no_na[[z]]$L2FC)
#   names(nam) = as.character(merged_no_na[[z]]$entrezgene_id)
# }
# 
# for (i in pathways$Pathway){
#   print(as.character(i))
#   pathview(gene.data = P0_vs_P1, 
#            pathway.id = as.character(i), 
#            species = "ssc")}
# 


# OLD FASHIONED WAY

P0_vs_P1_in = read.csv("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_P0_vs_P1.csv", header = TRUE)
P0_vs_P1_Simple_in = P0_vs_P1_in[,c("X","log2FoldChange","padj")]
colnames(P0_vs_P1_Simple_in) = c("Gene","L2FC","padj")
head(P0_vs_P1_Simple_in)
#Convert to fold change
P0_vs_P1_Simple_in$FC =(2^P0_vs_P1_Simple_in$L2FC)
head(P0_vs_P1_Simple_in)

P0_vs_P2_in = read.csv("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_P0_vs_P2.csv", header = TRUE)
P0_vs_P2_Simple_in = P0_vs_P2_in[,c("X","log2FoldChange","padj")]
colnames(P0_vs_P2_Simple_in) = c("Gene","L2FC","padj")
head(P0_vs_P2_Simple_in)
#Convert to fold change
P0_vs_P2_Simple_in$FC =(2^P0_vs_P2_Simple_in$L2FC)
head(P0_vs_P2_Simple_in)

P0_vs_P3_in = read.csv("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_P0_vs_P3.csv", header = TRUE)
P0_vs_P3_Simple_in = P0_vs_P3_in[,c("X","log2FoldChange","padj")]
colnames(P0_vs_P3_Simple_in) = c("Gene","L2FC","padj")
head(P0_vs_P3_Simple_in)
#Convert to fold change
P0_vs_P3_Simple_in$FC =(2^P0_vs_P3_Simple_in$L2FC)
head(P0_vs_P3_Simple_in)

P0_vs_P4_in = read.csv("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_P0_vs_P4.csv", header = TRUE)
P0_vs_P4_Simple_in = P0_vs_P4_in[,c("X","log2FoldChange","padj")]
colnames(P0_vs_P4_Simple_in) = c("Gene","L2FC","padj")
head(P0_vs_P4_Simple_in)
#Convert to fold change
P0_vs_P4_Simple_in$FC =(2^P0_vs_P4_Simple_in$L2FC)
head(P0_vs_P4_Simple_in)

P0_vs_P5_in = read.csv("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_P0_vs_P5.csv", header = TRUE)
P0_vs_P5_Simple_in = P0_vs_P5_in[,c("X","log2FoldChange","padj")]
colnames(P0_vs_P5_Simple_in) = c("Gene","L2FC","padj")
head(P0_vs_P5_Simple_in)
#Convert to fold change
P0_vs_P5_Simple_in$FC =(2^P0_vs_P5_Simple_in$L2FC)
head(P0_vs_P5_Simple_in)

library(pathview)

annots =getBM(
    mart=mart,
    attributes=c("external_gene_name","hgnc_symbol","hgnc_id","ensembl_gene_id","entrezgene_id"),
    filter="external_gene_name",
    values=P0_vs_P1_in$X,
    uniqueRows=TRUE)

merged = merge(P0_vs_P1_in, annots, by.x = "X", by.y = "external_gene_name", simplify = FALSE)
merge(paste0(i,"_in"), annots, by.x = "X", by.y = "external_gene_name", simplify = FALSE)

list = c("P0_vs_P1","P0_vs_P2","P0_vs_P3","P0_vs_P4","P0_vs_P5")

for (i in list){
  print(i)
  print(paste0(i,"_in"))
  assign(paste0(i, "_annot"), merge(paste0(i,"_in"), annots, by.x = "X", by.y = "external_gene_name", simplify = FALSE))
  #assign(paste0(i, "_na_omit"), na.omit(i[,c("entrezgene_id","L2FC")]))
}






P0_vs_P1_in[,c("entrezgene_id","L2FC")]



KEGG = all_annots
head(KEGG)
d=na.omit(KEGG[,c("entrezgene_id","L2FC")])
d
## feature 1: numeric vector
#geneList = d[,2]
geneList = d$L2FC
## feature 2: named vector
names(geneList) = as.character(d$entrezgene_id)
## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)
#de <- names(geneList)[abs(geneList) > 2]
head(geneList)
barplot(sort(geneList, decreasing = T))

for (i in pathways$Immune.system){
  print(as.character(i))
  pathview(gene.data = geneList, 
           pathway.id = as.character(i), 
           species = "ssc")
}

# SCRIPT DUMP


### GSEA
#Set levels
#Set choice of; CC - Cellular Component, MF - Molecular Function, BP - Biological Process
# # Or all_annots
# GSEA = all_annots
# head(GSEA)
# d=na.omit(GSEA[,c("entrezgene_id","L2FC")])
# d
# ## feature 1: numeric vector
# #geneList = d[,2]
# geneList = d$L2FC
# ## feature 2: named vector
# names(geneList) = as.character(d$entrezgene_id)
# ## feature 3: decreasing orde
# geneList = sort(geneList, decreasing = TRUE)
# #de <- names(geneList)[abs(geneList) > 2]
# head(geneList)
# barplot(sort(geneList, decreasing = T))
# 
# ego_bp <- gseGO(geneList     = geneList,
#               OrgDb        = org.Ss.eg.db,
#               ont          = "BP",
#               keyType      = "ENTREZID",
#               nPerm        = 1000,
#               minGSSize    = 1,
#               #maxGSSize    = 500,
#               pvalueCutoff = 0.05,
#               pAdjustMethod = "BH",
#               by           = "fgsea",
#               verbose      = TRUE)
# dotplot(ego_bp, showCategory=10)
# dev.copy(png, width = 1000,height = 1000, "BP_GSEA_GO.png")
# dev.off()
# 
# ego_mf <- gseGO(geneList     = geneList,
#               OrgDb        = org.Ss.eg.db,
#               ont          = "MF",
#               keyType      = "ENTREZID",
#               nPerm        = 1000,
#               minGSSize    = 1,
#               #maxGSSize    = 500,
#               pvalueCutoff = 0.05,
#               pAdjustMethod = "BH",
#               by           = "fgsea",
#               verbose      = TRUE)
# dotplot(ego_mf, showCategory=10)
# dev.copy(png, width = 1000,height = 1000, "MF_GSEA_GO.png")
# dev.off()
# 
# ego_cc <- gseGO(geneList     = geneList,
#               OrgDb        = org.Ss.eg.db,
#               ont          = "CC",
#               keyType      = "ENTREZID",
#               nPerm        = 1000,
#               minGSSize    = 1,
#               #maxGSSize    = 500,
#               pvalueCutoff = 0.05,
#               pAdjustMethod = "BH",
#               by           = "fgsea",
#               verbose      = TRUE)
# dotplot(ego_cc, showCategory=10)
# dev.copy(png, width = 1000,height = 1000, "CC_GSEA_GO.png")
# dev.off()



# Load pathview
#BiocManager::install("pathview")

# NEEDS TO BE CONDUCTED PER COMPARISON
# 
# library(pathview)
# 
# kegg_enrich <- enrichKEGG(gene = names(geneList),
#                           organism = 'ssc',
#                           pvalueCutoff = 0.05, 
#                           qvalueCutoff = 0.10)
# # Plot results
# barplot(kegg_enrich, 
#         drop = TRUE, 
#         showCategory = 10, 
#         title = "KEGG Enrichment Pathways",
#         font.size = 8)
# 
# go_enrich <- enrichGO(gene = names(geneList),
#                       OrgDb = 'org.Ss.eg.db', 
#                       readable = T,
#                       ont = "BP",
#                       pvalueCutoff = 0.05, 
#                       qvalueCutoff = 0.10)
# barplot(go_enrich, 
#         drop = TRUE, 
#         showCategory = 10, 
#         title = "GO Biological Pathways",
#         font.size = 8)


# P0_vs_P1_in = read.csv("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_P0_vs_P1.csv", header = TRUE)
# P0_vs_P1_Simple_in = P0_vs_P1_in[,c("X","log2FoldChange","padj")]
# colnames(P0_vs_P1_Simple_in) = c("Gene","L2FC","padj")
# head(P0_vs_P1_Simple_in)
# #Convert to fold change
# P0_vs_P1_Simple_in$FC =(2^P0_vs_P1_Simple_in$L2FC)
# head(P0_vs_P1_Simple_in)
# 
# P0_vs_P2_in = read.csv("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_P0_vs_P2.csv", header = TRUE)
# P0_vs_P2_Simple_in = P0_vs_P2_in[,c("X","log2FoldChange","padj")]
# colnames(P0_vs_P2_Simple_in) = c("Gene","L2FC","padj")
# head(P0_vs_P2_Simple_in)
# #Convert to fold change
# P0_vs_P2_Simple_in$FC =(2^P0_vs_P2_Simple_in$L2FC)
# head(P0_vs_P2_Simple_in)
# 
# P0_vs_P3_in = read.csv("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_P0_vs_P3.csv", header = TRUE)
# P0_vs_P3_Simple_in = P0_vs_P3_in[,c("X","log2FoldChange","padj")]
# colnames(P0_vs_P3_Simple_in) = c("Gene","L2FC","padj")
# head(P0_vs_P3_Simple_in)
# #Convert to fold change
# P0_vs_P3_Simple_in$FC =(2^P0_vs_P3_Simple_in$L2FC)
# head(P0_vs_P3_Simple_in)
# 
# P0_vs_P4_in = read.csv("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_P0_vs_P4.csv", header = TRUE)
# P0_vs_P4_Simple_in = P0_vs_P4_in[,c("X","log2FoldChange","padj")]
# colnames(P0_vs_P4_Simple_in) = c("Gene","L2FC","padj")
# head(P0_vs_P4_Simple_in)
# #Convert to fold change
# P0_vs_P4_Simple_in$FC =(2^P0_vs_P4_Simple_in$L2FC)
# head(P0_vs_P4_Simple_in)
# 
# P0_vs_P5_in = read.csv("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_P0_vs_P5.csv", header = TRUE)
# P0_vs_P5_Simple_in = P0_vs_P5_in[,c("X","log2FoldChange","padj")]
# colnames(P0_vs_P5_Simple_in) = c("Gene","L2FC","padj")
# head(P0_vs_P5_Simple_in)
# #Convert to fold change
# P0_vs_P5_Simple_in$FC =(2^P0_vs_P5_Simple_in$L2FC)
# head(P0_vs_P5_Simple_in)

# names(all_input[1])
# unlist(all_input$DESeq2_P0_vs_P1[["X"]])
# Can be accessed by;
# comparison = "P0_vs_P1"
# paste(comparison,"_input", sep="") = all_input[[comparison]]
# 
# lapply(all_input, assign(paste(p,"_simple_in", sep=""), lapply(p, `[`, c("X","log2FoldChange","padj"))))
# all_input$DESeq2_P0_vs_P1[["X"]]

# for (p in all_input) {
#     all_input$p[["X"]]
#   #print(paste("Comparison of ", p[1]," vs ", p[2]," begun"))
#   #assign(paste("DESeq2_Results_",p[1],"_vs_",p[2], sep=""), results(ddskeep,contrast=c("phase2_cat",p[1],p[2]), alpha = 0.05))
#   #results(ddskeep,contrast=c("phase2_cat",p[1],p[2]), alpha = 0.05)
#   #order results by padj (p-adjusted) value (most significant to least)
#   #assign(paste("DESeq2_Results_Ordered_",p[1],"_vs_",p[2], sep=""), name[order(name["padj"]),])
#   # Merge with counts for full output table
#   #assign(paste("DESeq2_Results_Ordered_Full_",p[1],"_vs_",p[2], sep=""),merge(as.data.frame(name2), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE))
#   #write.csv(name[order(name["padj"]),], file = paste("DESeq2_OutputTables/","DESeq2_Results_Ordered_",p[1],"_vs_",p[2], ".csv",sep=""))
#   #Remove all under 0.05
#   #res1= subset(res1, padj<0.05)  #No need to remove all under 0.05
# }
# 
# as.data.frame(counts(dds,normalized =TRUE))
# x =counts(dds,normalized =TRUE)
# counts(DESeq2_Results_P0_vs_P1)
# getwd()
# 
# 
# comparison = "P0_vs_P1"
# assign(paste0(comparison,"_input", sep=""), read.csv(paste0("../DESeq2_output/replicate + phase2_cat/DESeq2_LRT_OutputTables/DESeq2_", comparison,".csv", sep=""), header = TRUE))
# assign(paste0(comparison,"_simple", sep=""), (paste0(comparison,"_input", sep="")[c("X","log2FoldChange","padj"),]))
# 
# head(P0_vs_P1_input[,c("X","log2FoldChange","padj")])
# paste0(comparison,"_input", sep="")
# head(paste0(comparison,"_input", sep="")[,c("X","log2FoldChange","padj")])
# 
# P0_vs_P1_input[c("X","log2FoldChange","padj"),]
