##########################
#author: Amparo Roig Adam
#last updated: 2021-03-06
#input: 
# - Raw count data expression matrices and metadata from Tasic et. al (2018) - https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq
#output:
# - filtered all neuronal count data and metadata 
# - filtered synaptic genes neuronal count data
##########################
#Read Data TASIC 2018 ##########################

options(stringsAsFactors = FALSE)

ALMdata<-read.csv(file="./mouse_ALM_gene_expression_matrices_2018-06-14/mouse_ALM_2018-06-14_exon-matrix.csv")
ALMdata01<-ALMdata[,-1]
ALMsamples<-read.csv(file="./mouse_ALM_gene_expression_matrices_2018-06-14/mouse_ALM_2018-06-14_samples-columns.csv")
ALMgenes<-read.csv(file="./mouse_ALM_gene_expression_matrices_2018-06-14/mouse_ALM_2018-06-14_genes-rows.csv")

VISpdata<-read.csv(file="./mouse_VISp_gene_expression_matrices_2018-06-14/mouse_VISp_2018-06-14_exon-matrix.csv")
VISpdata01<-VISpdata[,-1]
VISpsamples<-read.csv(file="./mouse_VISp_gene_expression_matrices_2018-06-14/mouse_VISp_2018-06-14_samples-columns.csv")
VISpgenes<-read.csv(file="./mouse_VISp_gene_expression_matrices_2018-06-14/mouse_VISp_2018-06-14_genes-rows.csv")

data<-cbind(ALMdata01,VISpdata01)
row.names(data)<-VISpgenes[,"gene_symbol"]

metadata<-read.csv(file="./41586_2018_654_MOESM3_ESM/Supplementary_Table_10_Full_Metadata.csv",sep=";")

#order data
DataOrdered<-data[,metadata[,1]]
#Checkpoint
numero<-5666
metadata[numero,1]
colnames(DataOrdered)[numero]
#Save data
save(DataOrdered,metadata,file = "ALM_VISp_data.RData")

# Filter Data. Filter out the non-neuronal cells from the data##########################

table(metadata[,"class"])
NeuronalData<-DataOrdered[,metadata[,"class"]=="GABAergic" | metadata[,"class"]=="Glutamatergic"]
NeuronalMetadata<-metadata[metadata[,"class"]=="GABAergic" | metadata[,"class"]=="Glutamatergic",]

df <- toupper(row.names(NeuronalData))   
row.names(NeuronalData)<-df

AllSynData <- NeuronalData[row.names(NeuronalData) %in% SynGenes,]

#Save filtered and SynGO data
save(NeuronalData, NeuronalMetadata,file = "TASIC/data/03_Neuronal_Data.RData")
save(AllSynData,file = "TASIC/data/03_Synaptic_Data.RData")
