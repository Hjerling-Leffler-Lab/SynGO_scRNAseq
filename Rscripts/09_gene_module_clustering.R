##########################
#description: Classification of WGCNA according to their variance across cell types
#author: Amparo Roig Adam
#last updated: 2022-07-13
#input: 
# - Seurat object containing the count data of all SynGO genes
# - WGCNA modules and hclust dendrogram from  04_wgcna
#output:
# - k-means clustering of WGCNA modules according to their variance across neuronal clusters
# - PC1 density plot with clustered gene modules
# - dendrogram with previously generated gene modules and corresponding classification 
# - scatter plot showing two anticorrelated gene modules
##########################
library(WGCNA)
library(cowplot)
#opts
theme_set(theme_minimal())
theme_update(plot.background = element_rect(fill = "transparent", color = NA))
options(stringsAsFactors = FALSE)

#load data
load("TASIC0/results/02_clusterings/04_clustering_Syn.RData")
object <- UpdateSeuratObject(Sobject)

#calculate gene module average expression per cell
MEs0 = moduleEigengenes(t(as.matrix(object@assays$RNA@counts)), new.final.colors, excludeGrey = T)
avg.exp.mods <- MEs0$averageExpr

#calculate gene module variance across clusters
var.clust.mods <- matrix(nrow = length(unique(object$cluster)), ncol = length(unique(new.final.colors[new.final.colors!='grey'])))
colnames(var.clust.mods) <- unique(new.final.colors[new.final.colors!='grey'])
row.names(var.clust.mods)<- unique(object$cluster)

for (genemodule in unique(new.final.colors[new.final.colors!='grey'])){
  for(cluster in unique(object$cluster)){
    aegenemodule <- paste('AE', genemodule, sep ='')
    indat <- avg.exp.mods[row.names(object@meta.data[object$cluster == cluster,]),aegenemodule]
    var.clust.mods[cluster,genemodule] <- var(indat)
  }
}

#remove cell type identity and sort in a decreasing fashion
dat <- apply(t(var.clust.mods), 1, sort, decreasing = T)
row.names(dat)<- paste('top', 1:117, sep ='')

#PCA on the variance of each module in each cluster to find clusters of gene modules
library(ggbiplot)
#run pca on ordered matrix
pca.modules <- prcomp(t(dat), scale = T)#prcomp(apply(moduleTraitCor, 1, sort, decreasing = T), scale = T)
summary(pca.modules)
#loadings
pca.modules[["rotation"]][,'PC1'][c(order(pca.modules[["rotation"]][,'PC1'])[1:5],order(decreasing = T,pca.modules[["rotation"]][,'PC1'])[1:5])]

#pca.modules$rotation<- pca.modules$rotation[1:18,]

pca.plot<- ggbiplot(pca.modules, choices = c(1,2),groups = row.names(pca.modules$x), var.axes = F) +
  scale_color_manual(values = gsub(pattern = 'AE', '',  row.names(pca.modules$x)[order(row.names(pca.modules$x))]))
pca.plot+ggtitle('Merged modules')

#Cluster modules
#DelayedArray::seed(100)
genmods.kmeans <- kmeans(pca.modules$x[,1], centers = 3, iter.max = 100)
genmods.kmeans$cluster
kmean.plot<- ggplot(as.data.frame(cbind(pca.modules$x[,1:2],cluster = genmods.kmeans$cluster)), aes(x = PC1, y = PC2, color = as.character(cluster))) +geom_point()

plot_grid(pca.plot+ggtitle('Merged Modules')+NoLegend(), kmean.plot+NoLegend(), ncol = 1)

#density plot of pc1
library(ggplot2)
library(dplyr)
set.seed(123)

d <- ggplot(as.data.frame(cbind(pca.modules$x[,1],cluster = genmods.kmeans$cluster, modname = names(genmods.kmeans$cluster))),
       aes(x = as.numeric(V1)))+ylab('density')+xlab('PC1 (82.5% explained var.)')+
  geom_density( fill="lightgrey", color="#e9ecef", alpha=0.8, bw = 2)+
  geom_jitter(aes( y = 0.01, color = cluster), height = 0.01, width = 0)+NoLegend()+
  theme(axis.line = element_line(size = 0.5, colour = "black"),axis.ticks = element_line(size = 0.5, color="black") )
#+scale_color_manual(values = gsub(pattern = 'AE', '',  row.names(pca.modules$x)[order(row.names(pca.modules$x))]))
ggsave('TASIC0/results/wgcna/finalmergeed/modclassif_eigenvector.svg', width = 60, height = 60, units = 'mm')


### Compare apparent anticorrelated gene modules - yellow and turquoise
##Module Average Expression
ModuleAvgExpression <- function(object = NA, module = NA, nmod = 0){
  module = as.character(module)
  # Select genes
  gene.set <- get(paste('genes_in_modules', nmod, sep = ''))[,module]
  gene.set <- gene.set[is.na(gene.set) == F]
  gene.set = gene.set[gene.set %in% row.names(object@assays[["RNA"]]@data) == TRUE]
  gene.set = gene.set[gene.set %in% row.names(object@assays[["RNA"]]@data) == TRUE]
  #Normalize
  norm.exp <- object@assays[["RNA"]]@data[gene.set, ]/max(object@assays[["RNA"]]@data[gene.set, ])
  max(norm.exp)
  # Get mean expression of genes of interest per cell
  mean.exp <- colMeans(x = as.matrix(norm.exp[gene.set, ]), na.rm = TRUE)
  # Add mean expression values in 'object@meta.data$gene.set.score'   DO THIS FOR EACH MODULE SI THAT IT CAN BE RETRIEVED LATER (assign)
  modname = module#paste(module, mods, sep = ".")
  object@meta.data$gene.set.score <- mean.exp
  #FeaturePlot(object = object, features = "gene.set.score",cols = c("light grey","red")) + ggtitle(module)
  return(mean.exp)
}

object@meta.data$yellow <- ModuleAvgExpression(object, module = 'yellow', nmod = '.merged')
object@meta.data$turquoise <- ModuleAvgExpression(object, module = 'turquoise', nmod = '.merged')

c <- ggplot(object@meta.data , aes(x = yellow, y = turquoise))+
  geom_point(alpha = 0.4)+geom_smooth(method=lm , color="red", se=FALSE)+ 
  stat_cor(method = "pearson",label.x = 0.45, label.y.npc = "top")+
  theme(text = element_text(size = 5),axis.line = element_line(size = 0.25, colour = "black"),axis.ticks = element_line(size = 0.25, color="black") )

cowplot::plot_grid(d+theme(text = element_text(size = 5),axis.line = element_line(size = 0.25, colour = "black"),axis.ticks = element_line(size = 0.25, color="black") )
                   , c, nrow = 1)
#save figure
ggsave('TASIC0/results/wgcna/finalmergeed/Fig3BC_modcorr.svg', width = 100, height = 50, units = 'mm')

#save data
save(pca.modules, genmods.kmeans, var.clust.mods, kmean.plot, pca.plot, file = 'TASIC0/results/wgcna/finalmergeed/module_classification.RData')

#dendrogram with gene module classification
classif.colors <- merged.colors
classif.colors[classif.colors%in%names(genmods.kmeans$cluster[genmods.kmeans$cluster == 1])] <- '#F9766D'
classif.colors[classif.colors%in%names(genmods.kmeans$cluster[genmods.kmeans$cluster == 2])] <- '#00BA38'
classif.colors[classif.colors%in%names(genmods.kmeans$cluster[genmods.kmeans$cluster == 3])] <- '#619CFF'

svg('TASIC0/results/wgcna/finalmergeed/dendro.svg', width = 8, height = 5)
plotDendroAndColors(dendro = geneTree , cbind(merged.colors, classif.colors),
                    c("Merged modules", "Classification"),
                    dendroLabels = F, cex.dendroLabels = 0.25,hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05,
                    main = '')
dev.off()


