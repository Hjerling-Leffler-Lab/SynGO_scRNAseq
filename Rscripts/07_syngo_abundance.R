##########################
#description: Testing whether abundance of transcripts in each SynGO term is correlated with MN performance. Generates data needed for Supplementary Figure 2 D,E
#author: Amparo Roig Adam
#last updated: 2022-10-26
#input: 
# - Seurat object containing the count data of all SynGO genes
# - SynGO annotations and all SynGO genes
#output:
# - matrix with he average gene expression for each SynGO term to be used in the colour-coding tool from SynGO
##########################
 
library(SummarizedExperiment)
library(MetaNeighbor)
library(dplyr)
library(ggplot2)
library(rjson)
library(reshape2)
library(Seurat)

theme_set(theme_classic())
theme_update(plot.background = element_rect(fill = "transparent", color = NA, plit))

#Data - Syngo genes per SynGO term - all included in each term (children will overlap)
load("TASIC0/results/02_clustering/04_clustering_Syn.RData")
allsynaptic <- UpdateSeuratObject(Sobject)

BP_annotations <- fromJSON(file = "data/SynGO_genesets_2018-6-21/BP_annotations.json")
CC_annotations <- fromJSON(file = "data/SynGO_genesets_2018-6-21/CC_annotations.json")

SynGOsets <- c(BP_annotations$aggregate$hgnc_symbol, CC_annotations$aggregate$hgnc_symbol)
SynGOsets <- lapply(SynGOsets, function(x) x[x%in%rownames(allsynaptic)]) #eliminate genes not in dataset
SynGOsets <- Filter(function(x) length(x) >= 2, SynGOsets) #eliminate empty GO terms or with less than 2 genes

#Calculate abundance per SynGO set as avrg expression and total sum
allsynaptic$dummyident <- 'all'
avgSynGOs <- as.matrix(lapply(SynGOsets, function(x){
  mean(AverageExpression(allsynaptic, features = unlist(x), group.by = 'dummyident', verbose = F)$RNA)
}))
avgSynGOs = avgSynGOs[is.na(avgSynGOs)==FALSE,]

#save plot to colourcode Sunburst
write.csv(as.matrix(avgSynGOs), file = 'TASIC0/results/avgSynGOs.txt',quote = F, row.names = T, col.names = F )
