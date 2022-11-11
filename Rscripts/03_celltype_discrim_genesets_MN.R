##########################
#author: Amparo Roig Adam
#last updated: 2022-11-09
#input: 
# - fulldataset = Seurat object with the original data and metadata
# - mn.genesets = gene sets of HVGs used to generate the tSNE plots in 02_seurat3_analysis.R
# - mitoCarta = all genes in MitoCarta 3.0 - http://www.broadinstitute.org/files/shared/metabolism/mitocarta/mouse.mitocarta3.0.html
#output:
# The following files are saved for each geneset provided
# - MetaNeighbor AUROC scores of each gene set and for each cell class or type
# - Statistical test results comparing all gene sets tested
# - Boxplot of all results
##########################

library(SummarizedExperiment)
library(MetaNeighbor)
library(dplyr)
library(ggplot2)
library(rjson)
library(reshape2)
library(sdamr)
library(ggpubr)

theme_set(theme_minimal())
theme_update(plot.background = element_rect(fill = "transparent", color = NA))

#Seurat all data
load("TASIC0/results/02_clustering/fulldataset.RData") #this should be a seurat object containing all the data and metadata
fulldataset <- UpdateSeuratObject(fulldataset)

#cluster colors
cols <- read.csv('TASIC0/data/sample_heatmap_plot_data.csv', row.names = 1)
cluster_colors <- unique(as.character(cols$cluster_color[cols$class_label %in% c('GABAergic', 'Glutamatergic')]))
names(cluster_colors) <- unique(cols$cluster_label[cols$class_label %in% c('GABAergic', 'Glutamatergic')])
#names(cluster_colors) <- gsub(pattern = ' ', x = names(cluster_colors), replacement = '.')
#names(cluster_colors) <- gsub(pattern = '/', x = names(cluster_colors), replacement = '.')

#class colours - class colours are chosen as the cluster colour for the most abundant cell type withinthe class
most.abundant.celltype.per.class <- data.frame(melt(prop.table(table(fulldataset$cluster, fulldataset$subclass), 2))) %>% group_by(Var2) %>% slice(which.max(value))
class_colors <- cluster_colors[gsub('/','/',gsub(most.abundant.celltype.per.class$Var1, pattern = ' ', replacement = ' '))]
names(class_colors) <- gsub(' ', ' ',most.abundant.celltype.per.class$Var2)
#names(class_colors) <- gsub('/', '.',names(class_colors))

#Dataset#######
load("TASIC0/data/00_Neuronal_Data.RData")
Data = NeuronalData
Metadata = NeuronalMetadata

#Gene sets######
load('TASIC0/results/02_clustering/allgenesets4mn.RData') #gene sets of HVGs used to generate the tSNE plots in 02_seurat3_analysis.R

#filter to exclusively pre/post synaptic genes AND present in the dataset
filt.PreSynGenes = mn.genesets$presynaptic_2000[mn.genesets$presynaptic_2000%in%mn.genesets$postsynaptic_2000==F& mn.genesets$presynaptic_2000%in%mn.genesets$allsynaptic_2000]
filt.PostSynGenes = mn.genesets$postsynaptic_2000[mn.genesets$postsynaptic_2000%in%mn.genesets$presynaptic_2000 == F&mn.genesets$postsynaptic_2000%in%mn.genesets$allsynaptic_2000]

#Extra control - mitochondrial genes
mitoCarta <- readxl::read_xls('data/mouseMitoCarta/Mouse.MitoCarta3.0.xls', sheet = 2)
  #get gene list with symbol matching scRNA-seq data
mitoCarta.filt<- mitoCarta[toupper(mitoCarta$Symbol)%in%row.names(Data),]
mitogenes <- row.names(Data)[row.names(Data)%in%toupper(mitoCarta$Symbol)]

#merge gene sets to quantify
  # for quantification of the cell class/type identity encoded in the genes underlying the tSNE plots in Fig1 + Non synaptic and Mitochondrial gene as alternative gene sets
GeneSets1 <- list(all_2000 = mn.genesets$all_2000, non_syn = mn.genesets$all_2000[mn.genesets$all_2000%in%mn.genesets$allsynaptic_2000==F], allsyn = mn.genesets$allsynaptic_2000, presyn = mn.genesets$presynaptic_2000, postsyn = mn.genesets$postsynaptic_2000, mito = mitogenes)
  # for quantification of cell type identity encoded at the pre and postsynapse independently in excitatpry and inhibitory cells
GeneSets2 <- list( allsyn = mn.genesets$allsynaptic_2000, presynaptic = mn.genesets$presynaptic_2000, postsynaptic = mn.genesets$postsynaptic_2000)
GeneSets3 <- list( allsyn = mn.genesets$allsynaptic_2000, presynaptic = filt.PreSynGenes, postsynaptic = filt.PostSynGenes)

#Build Metaneighbour functions to run, plot and save in desired combination of parameters
MN_calc = function(Data, Metadata, annot.level, GeneSets){
  
  print(paste0('Metaneighbor analysis on ', nrow(Metadata), ' cells'))
  print(paste0('Annotation level ', annot.level))
  print(paste0('GeneSets used: ', paste0(names(GeneSets), collapse= ', ')))
  
  set.seed(123)#to make the sampling (shuffling od study id) rerpoducible
  study_id = sample(matrix(data = c(rep("train", (nrow(Metadata)/3)*2+1),rep("test", nrow(Metadata)/3)), ncol = 1))
  SE_colData <- cbind(Metadata[,c("sample_name", "class","subclass","cluster")], study_id)
  
  m <- data.frame(as.character(Metadata[,"sample_name"]), as.character(Metadata[,annot.level])) ##change this for each of 3 levels analysis
  cell_labels <- table(m)
  
  SynData <- Data[row.names(Data) %in% unique(unlist(GeneSets)),Metadata$sample_name] ####GeneSets$AllGenes
  
  
  SE_SynData <- SummarizedExperiment(assays=list(counts=as.matrix(SynData)), colData = SE_colData ,metadata = Metadata)
  
  
  MN_result = MetaNeighbor(dat = SE_SynData,
                           experiment_labels = as.numeric(factor(SE_SynData$study_id)),
                           celltype_labels = cell_labels, genesets = GeneSets,
                           bplot = F, fast_version = T)
  return(MN_result)
}

MN_plot = function(AUROCs, annot.level, stats){
  if(annot.level == 'subclass'){
    group_cols = class_colors[as.character(melt(as.matrix(AUROCs))$Var2)]
  }
  if(annot.level == 'cluster'){
    group_cols = cluster_colors[as.character(melt(as.matrix(AUROCs))$Var2)]
  }
  p.top <-melt(as.matrix(AUROCs)) %>%
    ggplot(aes(x = Var1, y = value))+
    geom_boxplot(outlier.colour = NA, show.legend = F, notch = T)+
    geom_jitter(position=position_jitter(width = 0.25),alpha=0.6,size = 1,show.legend = F, aes(color = Var2)) + 
    scale_color_manual(values = group_cols)+
    ylab('AUROC scores')+xlab('') + 
    scale_y_continuous(limits = c(min(AUROCs)-0.01,NA))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1 , size = 10), axis.line = element_line(size = 0.25))+
    stat_pvalue_manual(stats[stats$p.adj<0.05,], label = 'p.signif', 
                       y.position = c(1.001), tip.length = 0.005, step.increase = 0.05,size = 3)
  p.down <- melt(as.matrix(AUROCs)) %>%
    ggplot(aes(x = Var1, y = value))+
    geom_boxplot(outlier.colour = NA, show.legend = F, notch = T)+
    geom_jitter(position=position_jitter(width = 0.25),alpha=0.6, size = 1,show.legend = F, aes(color = Var2)) +
    scale_y_continuous(limits = c(0.50,0.50), breaks = 0.5, labels = scales::number_format(accuracy = 0.01))+
    ylab('')+xlab('')+theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 0.25))
  
  return(ggarrange(p.top+theme(axis.text.x = element_blank(), axis.line.x = element_line(size =0.25)),p.down, ncol = 1, heights = c(4,1)))
}

MN_save = function(AUROCs, GeneSets, stats, result_plot, geneset.name, annot.level, subset = '',outdir = 'TASIC0/results/03_cellgroup_discriminability/'){
  write.csv(AUROCs, file = paste0(outdir,gsub('-', '', Sys.Date()),'_',annot.level, '_', geneset.name,subset, '_metaneighbor.csv'), quote = F)
  ggsave(result_plot, filename = paste0(outdir,gsub('-', '', Sys.Date()),'_',annot.level, '_', geneset.name,subset, '_metaneighbor.svg'),device= 'svg', bg = 'transparent', width = 4, height = 5)
  write.csv(stats, file = paste0(outdir,gsub('-', '', Sys.Date()),'_',annot.level, '_', geneset.name,subset, '_metaneighbor_stats.csv'), quote = F)
  save(GeneSets, AUROCs, stats, file = paste0(outdir,gsub('-', '', Sys.Date()),'_',annot.level, '_', geneset.name,subset, '_metaneighbor_results.RData'))
}

#Run metaneighbor analysis on the collected gene sets
for(gs in c('GeneSets1', 'GeneSets2', 'GeneSets3')){
  GeneSets = get(gs)
  if(gs == 'GeneSets1'){ #run analysis in all data for subclass and cluster level of annotation
    for(annot.level in c('subclass', 'cluster')){
      AUROC_scores = MN_calc(Data = NeuronalData, Metadata = NeuronalMetadata, annot.level = annot.level, GeneSets = GeneSets)
      
      stats <- compare_means(value ~ Var1,  data = melt(as.matrix(AUROC_scores)), method = 'wilcox.test', paired = F, p.adjust.method = 'fdr')
      print(stats)
      result_plot = MN_plot(AUROCs = AUROC_scores, annot.level = annot.level , stats = stats)
      print(result_plot)
      MN_save(AUROCs = AUROC_scores, GeneSets = GeneSets, stats = stats, result_plot = result_plot, geneset.name = gs, annot.level = annot.level)
    }
  }
  if(gs%in%c('GeneSets2', 'Genesets3')){ #run analysis in GABA and glutamatergic cells independently for cluster level annotation only
    datasets = c('GABAergic', 'Glutamatergic')
    annot.level = c('subclass', 'cluster')
    for(subset in datasets){
      print(subset)
      annot.level = 'cluster'
      AUROC_scores = MN_calc(Data = NeuronalData, Metadata = NeuronalMetadata[NeuronalMetadata$class == subset,], annot.level = annot.level, GeneSets = GeneSets)
        
      stats <- compare_means(value ~ Var1,  data = melt(as.matrix(AUROC_scores)), method = 'wilcox.test', paired = F, p.adjust.method = 'fdr')
      print(stats)
      result_plot = MN_plot(AUROCs = AUROC_scores, annot.level = annot.level , stats = stats)
      print(result_plot)
      MN_save(AUROCs = AUROC_scores, GeneSets = GeneSets, stats = stats, result_plot = result_plot, geneset.name = gs, annot.level = annot.level, subset = subset)
    }}
  rm(AUROC_scores, stats, result_plot)
  gc()
  }
