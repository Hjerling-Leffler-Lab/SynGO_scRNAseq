##########################
#description: In this script we generate random gene sets (of all expressed or only synaptic genes) of the same sizes as the GO terms in SynGO and run Metaneighbour analysis.
#author: Amparo Roig Adam
#last updated: 2021-02-25
#input: 
# - NeuronalData; NeuronalMetadata
# - SynGO annotations and all SynGO genes
#output:
# - MetaNeighbor AUROC scores of each random gene set and for each cell class
##########################

library("rjson")
library('MetaNeighbor')
library('SummarizedExperiment')

### load data
load("TASIC/data/03_Neuronal_Data.RData")
Data = NeuronalData
Metadata = NeuronalMetadata
rm(NeuronalData, NeuronalMetadata)

### Retreive all SynGO genes

BP_annotations <- fromJSON(file = "data/SynGO_genesets_2018-6-21/BP_annotations.json")
CC_annotations <- fromJSON(file = "data/SynGO_genesets_2018-6-21/CC_annotations.json")

SynGenes_BP<-BP_annotations$aggregate$hgnc_symbol$'SYNGO:synprocess'
SynGenes_CC<-CC_annotations$aggregate$hgnc_symbol$'GO:0045202'

SynGOGenes <- unique(c(SynGenes_BP, SynGenes_CC))

  #keep only those present in the dataset - this will avoid problems with metaneighbour since it relies on correlations of gene expression and can't have NAs in the count matrix
SynGOGenes <- SynGOGenes[SynGOGenes %in% row.names(Data) == TRUE]

#Uncomment to perform the analysis with genes present in the dataset instead of syngo - random sets from all expressed genes
#SynGOGenes <- row.names(Data)

  #Filter genes to have expression in more than 10% of the cells
SynGOGenes <- SynGOGenes[rowSums(Data != 0) > ncol(Data)/10]

###Generate Randomised sets
## Gene set sizes == sizes of SynGO terms.  This is done to have random sets of the same size as the GO terms in SynGO
GOsyn <- c(BP_annotations$aggregate$hgnc_symbol, CC_annotations$aggregate$hgnc_symbol)
GOsyn <- Filter(function(x) length(x) >= 3, GOsyn) #filter too small SynGO terms
sizes <- unique(unlist(lapply(GOsyn, length)))

#Clean up unncessary objects
rm(BP_annotations, CC_annotations, SynGenes_BP, SynGenes_CC, GOsyn)

  #seed for reproducible sampling
set.seed(123)
##Random sets
#number of sets of each size to generate
N = 10000
randomSets <- replicate(n = N, expr = lapply(sizes, function(y) {sample(SynGOGenes, size = y)}), simplify = FALSE)
randomSets <- unlist(randomSets, recursive = F)

#give a name to the gene sets - necessary for MetaNeighbor
names(randomSets) <- paste('ranset', 1:length(randomSets), sep = '')

###METANEIGHBOR of all cell subclasses
set.seed(123)
study_id = sample(matrix(data = c(rep("train", nrow(Metadata)/3),rep("test", (nrow(Metadata)/3)*2+1)), ncol = 1))
SE_colData <- cbind(Metadata[,c("sample_name", "class","subclass","cluster")], study_id)
m <- data.frame(as.character(Metadata[,"sample_name"]), as.character(Metadata[,"subclass"]))
cell_labels <- table(m)
  #only keep expression data of SynGO genes
SynData <- Data[row.names(Data) %in% SynGOGenes,]
SE_SynData <- SummarizedExperiment(assays=list(counts=as.matrix(SynData)), colData = SE_colData ,metadata = Metadata)

#Genesets are in randomSets

#Run metaneighbour
AUROC_scores_subclass = MetaNeighbor(dat = SE_SynData,
                                    experiment_labels = as.numeric(factor(SE_SynData$study_id)),
                                    celltype_labels = cell_labels, genesets = randomSets,
                                    bplot = F, fast_version = T)

#keep size of each gene set
nGenes <- lapply(randomSets, length)
AUROC_scores_subclass <- cbind(t(data.frame(nGenes, check.names = F)), AUROC_scores_subclass)

#save random gene sets
save(randomSets, file = 'TASIC/results/01_metaneighbor/DATE_FILENAME.RData')

#Save table file with the obtained scores
write.csv(AUROC_scores_subclass, file = 'TASIC/results/01_metaneighbor/DATE_FILENAME.csv', quote = F)

