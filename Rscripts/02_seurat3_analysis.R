##########################
#author: Amparo Roig Adam
#last updated: 2022-11-09
#input: 
  # - geneset = List of gene sets to use containing the genes in each set
  # - NeuronalData and NeuronalMetadata = Raw count data previously filtered to all neurons in the dataset
#output:
  # - tSNE plots for each gene set in input
  # - seurat object containing normalised data and dimensionality reduction to obtain the tSNE
  # - mn.genesets = list of HVGs from each gene set to be used in 04_MN_celltype_discrim.R
  
##########################
library(Seurat)
library(cowplot)
library(ggalt)
library(ggplot2)
library(ggforce)
library(biomaRt)

theme_set(theme_classic())
theme_update(plot.background = element_rect(fill = "transparent", color = NA))

#cluster colors
cols <- read.csv('TASIC0/data/sample_heatmap_plot_data.csv', row.names = 1) #file obtained from Tasic et. al. (2018)
cluster_colors <- as.character(cols$cluster_color)
names(cluster_colors) <- cols$cluster_label

##gene set used
#geneset = c('all', 'allsynaptic', 'presynaptic', 'postsynaptic') these names should be present in the input list of gene sets
load("TASIC0/data/02_gensets.RData") #SynGO gene sets retrieved at 01_read_syngo.R
geneset = c(all = row.names(NeuronalData), allsynaptic = SynGenes[SynGenes%in%row.names(Data)] #keep only those in the data set
            , presynaptic = PreSynGenes, postsynaptic = PostSynGenes) 

#load raw data as filtered and formatted in script 00_read_filter_Tasic2018
load("TASIC0/data/01_Neuronal_Data.RData")

met<-NeuronalMetadata[,-1] #metadata rownames are in first column
row.names(met)<-NeuronalMetadata[,1]

#add cluster colours to metadata
met <- cbind(cols[row.names(met),]$cluster_color, met) 
colnames(met)[1] <- "cluster_color"

## run seurat for each geneset and save high variable genes for quantification in 03_MN_celltype_discrim.R
mn.genesets <- list()
for (gs in names(geneset)){
  n = 2000 #number of variable genes used
  print(paste0(gs,'-', n))
    #data
  Data <- NeuronalData[unlist(geneset[gs])[unlist(geneset[gs])%in%row.names(NeuronalData)],]
  print(table(row.names(Data)=='NA.1')) #check there is no NA rows added
    
  fulldataset <- CreateSeuratObject(counts = Data, meta.data = met, project = paste0("SynGO - ",gs,Sys.Date()))
  AddMetaData(fulldataset, metadata = cols[row.names(met),]$cluster_color, col.name = 'cluster_color')
    
    #Check whether there are NA counts
  print(length(names(which(rowSums(is.na(fulldataset@assays$RNA@counts))>0))))
    
    #preprocess and normalise
  fulldataset <- subset(fulldataset, subset = nFeature_RNA > 10 & nCount_RNA > 1000)
  fulldataset <- NormalizeData(fulldataset, normalization.method = "LogNormalize", scale.factor = 10000)
    
    #calculate HVGs
  fulldataset <- FindVariableFeatures(fulldataset,mean.cutoff = c(0.001,  8), dispersion.cutoff = c(-3, Inf),nfeatures = 2000)
    #VariableFeaturePlot(fulldataset)
    
  fulldataset <- ScaleData(fulldataset, features = row.names(fulldataset))
    
    ###Linear Dimensional reduction
  fulldataset <- RunPCA(fulldataset, features = VariableFeatures(fulldataset),npcs = 60) 
  #PCAPlot(object, dim.1 = 1, dim.2 = 2, pt.size = 0.5, no.legend = T)
    
    ##Determine statistical significant principal components
  #ElbowPlot(fulldataset, ndims = 60)
    
  #fulldataset <- JackStraw(fulldataset, dims = 60, num.replicate = 100)
  #JackStrawPlot(fulldataset, dims = 45:60)
    
      #Dims chosen: Syn --> 42; PreSyn and PostSyn --> 21; All genes --> 60   #this should be checked for other data sets used
  if(gs %in% c('all')){
    pcs = 60
  }
  if(gs %in% c("presynaptic", "postsynaptic")){
    pcs = 21
  }
  if(gs == 'allsynaptic'){
    pcs = 42
  }
  if(gs%in%c('all',"presynaptic", "postsynaptic", 'allsynaptic')==F){
    print('Undefined gene set for PC selection')
    break}
    #TSNE to visualize the clustering
  fulldataset <- RunTSNE(fulldataset,dims = 1:pcs, do.fast = T, perplexity = 50)
  DimPlot(fulldataset, reduction = 'tsne', group.by = 'cluster', shuffle = T, seed = 123, label = F, cols = cluster_colors)+ggtitle(paste0(gs, length(fulldataset@assays$RNA@var.features)))+NoAxes()+NoLegend()
  ggsave(filename = paste0('TASIC0/results/03_clustering/',gsub('-', '', Sys.Date()), '_', gs,'_2000.svg'),bg = 'transparent', width = 4,height = 4.25,dpi = 'retina')
  
  save(fulldataset, geneset,file = paste0('TASIC0/results/03_clustering/',gsub('-', '', Sys.Date()), '_', gs,'_2000.RData'), compress = T )
  
    #save HVGs used for next analysis
  hvg.genes <- list(fulldataset@assays$RNA@var.features)
  names(hvg.genes)<- paste0(gs, '_', n)
  mn.genesets <- c(mn.genesets, hvg.genes)
  }

save(mn.genesets, geneset,file = paste0('TASIC0/results/03_clustering/',gsub('-', '', Sys.Date()), '_allgenesets4mn.RData'), compress = T )

