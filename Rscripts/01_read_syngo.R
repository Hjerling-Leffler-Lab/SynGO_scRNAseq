##########################
#author: Amparo Roig Adam
#last updated: 2022-11-09
#input: 
# - SynGO_genesets_2018-6-21 from SynGO 1.0 - https://syngoportal.org/data/SynGO_bulk_download_release_20180731.zip
#output:
# - all synaptic, pre and post synaptic gene sets covered in SynGO
##########################

setwd("SynGO/TASIC0")
library("rjson")

# Read SynGO files
BP_annotations <- fromJSON(file = "data/SynGO_genesets_2018-6-21/BP_annotations.json")
BP_ontology <- fromJSON(file = "data/SynGO_genesets_2018-6-21/BP_ontology.json")
CC_annotations <- fromJSON(file = "data/SynGO_genesets_2018-6-21/CC_annotations.json")
CC_ontology <- fromJSON(file = "data/SynGO_genesets_2018-6-21/CC_ontology.json")

#Biological Proccess synaptic genes
SynGenes_BP<-BP_annotations$aggregate$hgnc_symbol$'SYNGO:synprocess'
PreSynGenes_BP <- BP_annotations$aggregate$hgnc_symbol$'SYNGO:presynprocess'
PostSynGenes_BP <- BP_annotations$aggregate$hgnc_symbol$'SYNGO:postsynprocess'

#Cellular Component synaptic genes
SynGenes_CC<-CC_annotations$aggregate$hgnc_symbol$'GO:0045202'
PreSynGenes_CC <- CC_annotations$aggregate$hgnc_symbol$'GO:0098793'
PostSynGenes_CC <- CC_annotations$aggregate$hgnc_symbol$'GO:0098794'

#All SynGO genes
SynGenes <- unique(c(SynGenes_BP, SynGenes_CC))

#All post and pre synaptic genes
PreSynGenes <- c(unique(x = c(PreSynGenes_BP, PreSynGenes_CC)))
PostSynGenes <- c(unique(x = c(PostSynGenes_BP, PostSynGenes_CC)))

#Save
save(SynGenes, PreSynGenes, PostSynGenes, file = "TASIC/data/02_SynGenes.RData")
