##########################
#description: Cell type Discriminatory power of each gene set annotated in SynGO
#author: Amparo Roig Adam
#last updated: 2021-02-25
#input: 
# - NeuronalData; NeuronalMetadata
# - SynGO annotations and all SynGO genes
#output:
# - MetaNeighbor AUROC scores of each gene set annotated in SynGO
##########################

library("rjson")
library('MetaNeighbor')
library('SummarizedExperiment')


#Load data generated in 00_read+filter_NeuronalData.R
load("TASIC0/data/00_Neuronal_Data.RData")

Data = NeuronalData
Metadata = NeuronalMetadata

rm(NeuronalData, NeuronalMetadata)

#load SynGO annotations
BP_annotations <- fromJSON(file = "data/SynGO_genesets_2018-6-21/BP_annotations.json")
CC_annotations <- fromJSON(file = "data/SynGO_genesets_2018-6-21/CC_annotations.json")

SynGenes_BP<-BP_annotations$aggregate$hgnc_symbol$'SYNGO:synprocess'
SynGenes_CC<-CC_annotations$aggregate$hgnc_symbol$'GO:0045202'

#Run Metaneighbor
#to do the analysis we need 'two' datasets because it measures how well it can replicate the same cell type annotation in both 
#we are not using two datasets, so we split the one at hand

set.seed(123)#to make the sampling (shuffling of study id) rerpoducible
study_id = sample(matrix(data = c(rep("train", nrow(Metadata)/3),rep("test", (nrow(Metadata)/3)*2+1)), ncol = 1))
SE_colData <- cbind(Metadata[,c("sample_name", "class","subclass","cluster")], study_id)

m <- data.frame(as.character(Metadata[,"sample_name"]), as.character(Metadata[,"subclass"]))
cell_labels <- table(m)

SynData <- Data[row.names(Data) %in% c(SynGenes_BP, SynGenes_CC),]

SE_SynData <- SummarizedExperiment(assays=list(counts=as.matrix(SynData)), colData = SE_colData ,metadata = Metadata)

#Gene set containing all SynGO terms to test
GOsyn <- c(BP_annotations$aggregate$hgnc_symbol, CC_annotations$aggregate$hgnc_symbol)
GOsyn <- Filter(function(x) length(x) >= 2, GOsyn) #eliminate empty GO terms or with less than 2 genes

rm(BP_annotations, CC_annotations, SynGenes_BP, SynGenes_CC)


AUROC_scores_subclas = MetaNeighbor(dat = SE_SynData,
                                    experiment_labels = as.numeric(factor(SE_SynData$study_id)),
                                    celltype_labels = cell_labels, genesets = GOsyn,
                                    bplot = TRUE, fast_version = F)

#keep size of each gene set
nGenes <- lapply(GOsyn, length)
AUROC_scores_subclass <- cbind(t(data.frame(nGenes, check.names = F)), AUROC_scores_subclass)

#save results
write.csv(AUROC_scores_subclass, file = 'TASIC/results/01_metaneighbor/metaneighbor_subclass_synGOauroc.csv', quote = F)

