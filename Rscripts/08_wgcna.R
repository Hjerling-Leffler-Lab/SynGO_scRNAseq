##########################
#description: Weighted gene correlation analysis to obtain gene modules and curation of those modules to merge similar ones
#author: Amparo Roig Adam
#last updated: 2021-11-15
#input: 
# - NeuronalData and NeuronalMetadata = Raw count data previously filtered in 00_read+filter_NeuronalData
# - SynGO genes
# - Seurat object containing the count data of all SynGO genes
#output:
# - table of gene modules obtained with the genes within each module
# - geneTree = hclust object containing the WGCNA results
# - tSNE plots showing the average expression of all genes in each gene module
##########################
library(Seurat)
library(ggplot2)
library(cowplot)
library(WGCNA)
options(stringsAsFactors = FALSE)

#Data needed
load("TASIC/data/01_Pre&Post_SynGenes.RData")
load("TASIC/data/00_Neuronal_Data.RData")

#Weighted gene correlation network analysis -  this analysis is based on the 
##I want to take the data with both the BP genes and the CC genes
data0<-NeuronalData[row.names(NeuronalData) %in% SynGenes,] ####

###Prefiltering----
##Filter out genes expressed in less than 10% of the cells
total = ncol(data0)
min = total* 0.1
genes <- c()
for (gene in row.names(data0)){
  n = 0
  for (value in data0[gene,]){
    if (value != 0){
      n = n+1
  }}
  if (n >= min){
    genes <- c(genes,gene)
  }}
rm(gene,min,n,total,value)

#Transpose the expression data for further analysis.
data = as.data.frame(t(data0))

##Check for excessive missing values. Shouldn't be any since both genes and samples were filtered before
gsg = goodSamplesGenes(data,verbose = 4)
gsg$allOK
#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(data)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(data)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  data = data[gsg$goodSamples, gsg$goodGenes]
}
###Network construction and module detection----
##Analysis of network topology (Choice of soft-thresholding power)
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#Call the network topology analysis function
sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(10, 5)
par(mfrow = c(1,2))
cex1 = 0.9
#Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
#This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
#Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h = 1, col = "red")
##Co-expression similarity and adjacency
softPower = 4
adjacency = adjacency(data, power = softPower)
##Topological Overlap Matrix
#Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
##Clustering using TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
geneTree$labels<-row.names(adjacency)
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Synaptic genes clustering",
     labels = NULL, cex = 0.4, hang = 0.04)

#Module identification using dynamic tree cut-----
dyn_modules1 = cutreeDynamic(dendro = geneTree, method = 'tree', cutHeight = '0.9', minClusterSize = 2)
table(dyn_modules1)
Colors1 = labels2colors(dyn_modules1)
table(Colors1)
dyn_modules2 = cutreeDynamic(dendro = geneTree, method = 'tree',cutHeight = '0.95', minClusterSize = 2)
table(dyn_modules2)
Colors2 = labels2colors(dyn_modules2)
table(Colors2)
dyn_modules3 = cutreeDynamic(dendro = geneTree, method = 'tree',cutHeight = '0.98', minClusterSize = 2)
table(dyn_modules3)
Colors3 = labels2colors(dyn_modules3)
table(Colors3)
dyn_modules4 = cutreeDynamic(dendro = geneTree, method = 'tree', cutHeight = '0.99', minClusterSize = 2)
table(dyn_modules4)
Colors4 = labels2colors(dyn_modules4)
table(Colors4)

#Plot the dendrogram and colors-------
plotDendroAndColors(dendro = geneTree , cbind(Colors1, Colors2, Colors3, Colors4),
                    c("Mods1-0.90","Mods2-0.95", "Mods3-0.98", "Mods4-0.99"),
                    dendroLabels = NULL, cex.dendroLabels = 0.25,hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Synaptic genes - Modules")


#Obtain labels and genes from modules
labels <- rev(geneTree$labels[geneTree$order])
id1 <- rev(Colors1[geneTree$order])
genesmod1 <- labels[id1 != "grey"]
id2 <- rev(Colors2[geneTree$order])
genesmod2 <- labels[id2 != "grey"]
id3 <- rev(Colors3[geneTree$order])
genesmod3 <- labels[id3 != "grey"]
id4 <- rev(Colors4[geneTree$order])
genesmod4 <- labels[id4 != "grey"]

###Write modules content in a file----
genes_in_modules <- c()
##Modules in the order they appear in the dendogram (color order)
cols <- c()
for(color in Colors4[geneTree$order]){  ###
  if( color %in% cols == F & color != "grey"){
    cols = c(cols,color)
  }
}
##Add each gene to the column of the table corresponding to its module. Right to left order
#for each module
for(color in cols){
  modulegenes <- c(labels[id4 == color])####
  modulegenes <- c(modulegenes, rep(NA, max(table(Colors4[Colors4!="grey"]))###
                                    - length(modulegenes)))
  genes_in_modules <- cbind(genes_in_modules, modulegenes)
}
colnames(genes_in_modules) <- paste(cols, 4, sep = ".") ##
genes_in_modules -> genes_in_modules4 ###
rm(cols, color, modulegenes, genes_in_modules)
write.table(genes_in_modules4,      ######
            file = "TASIC/results/wgcna/genes_in_modules4.txt", ###
            row.names = FALSE, col.names = T, sep = "\t", quote = FALSE, na="")

save(geneTree, labels, Colors1, Colors2,Colors3, Colors4,
     genes_in_modules1,  genes_in_modules2,genes_in_modules3, genes_in_modules4,file = "TASIC/results/wgcna/syn_gene_modules.RData")

### Curated final wgcna modules and hierarchical clustering tree
#modules obtained with dynamic tree cut are reviewd and selected according to the convergece of the 3 tested tree cuts
mods <- as.character(read.table('TASIC0/results/wgcna/curatedgenemodules.txt')[,1])
final.mods <- c()
for(module in mods){
  print(module)
  n = unlist(strsplit(module,split = '.', fixed = T))[2]
  module <- unlist(strsplit(module,split = '.', fixed = T))[1]
  #retrieve module
  mod <- get(paste('genes_in_modules', n, sep = ''))[,module]
  #adjust length - 96 is the maximum possible length found in all modules saved
  mod<- c(mod, rep(NA, 96-length(mod)))
  final.mods<- cbind(final.mods, mod)
}
colnames(final.mods)<- mods
write.table(final.mods,      ######
            file = "TASIC0/results/wgcna/curated_gene_mods.tsv", ###
            row.names = FALSE, col.names = T, sep = "\t", quote = FALSE, na="")

#new colors for curated modules
melt.final.mods <- melt(final.mods)
melt.final.mods <- melt.final.mods[is.na(melt.final.mods$value)==F, ]
melt.final.mods$Var1 <- as.numeric(melt.final.mods$Var2)

num.final.mods <- geneTree$labels
num.final.mods[num.final.mods%in%melt.final.mods$value == F] <- '0grey'
for(gene in num.final.mods[num.final.mods!= '0grey']){
  num.final.mods[num.final.mods%in%gene] = melt.final.mods$Var2[melt.final.mods$value == gene]
}
#num.final.mods[melt.final.mods$value] <- as.character(melt.final.mods$Var2[melt.final.mods$value==])
num.final.mods<- as.numeric(as.factor(num.final.mods)) -1

final.colors <- labels2colors(num.final.mods, zeroIsGrey = T )

table.final.colors <- data.frame(table(final.colors[final.colors!='grey']))%>%arrange(Freq)
table.final.colors<- cbind(table.final.colors, new.c = rev(standardColors(33)))

new.final.colors <- rep(NA, 1049)
for(org.col in unique(final.colors[final.colors!='grey'])){
  print(paste(org.col, 'changes to', table.final.colors$new.c[table.final.colors$Var1 == org.col]))
  new.final.colors[final.colors%in%org.col] <- table.final.colors$new.c[table.final.colors$Var1 == org.col]
}
new.final.colors[is.na(new.final.colors)] <- 'grey'

for(color in unique(new.final.colors)[unique(new.final.colors)!='grey']){
  modulegenes <- c(labels[new.final.colors == color])####
  modulegenes <- c(modulegenes, rep(NA, max(table(new.final.colors[new.final.colors!="grey"]))###
                                    - length(modulegenes)))
  final.mods.new <- cbind(final.mods.new, modulegenes)
}

#opt brown 
new.final.colors[Colors3 == 'blue'] <- 'brown4'
new.final.colors[Colors2 == 'yellow'] <- 'brown4'

  #final dendrogram
plotDendroAndColors(dendro = geneTree , cbind(new.final.colors),
                    c(""),
                    dendroLabels = NULL, cex.dendroLabels = 0.25,hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Synaptic gene modules")

#save file with modules
genes_in_modules.merged <- c()
colnam <- c()
for (modname in unique(new.final.colors)){
  if(modname != 'grey'){
    #genes en una lista de ordenadors por order
    modulegenes <- c(geneTree$labels[new.final.colors== modname])####
    modulegenes <- c(modulegenes, rep(NA, max(table(new.final.colors[new.final.colors!="grey"]))- length(modulegenes)))
    genes_in_modules.merged <- cbind(genes_in_modules.merged, modulegenes)
    colnam<- c(colnam, modname)
  }
}
colnames(genes_in_modules.fin) <- colnam
write.table(genes_in_modules.fin, 'TASIC0/results/wgcna/moduletable.csv', sep =',', row.names = F, na = '')

#save all produced objects
save.image(file = 'TASIC0/results/wgcna/finalmergeed/20220125_final_merged_modules.RData')

##Plot tSNE color coded with WGCNA clusters average expression----
#load Seurat object with data to plot tSNE
load("TASIC0/results/02_clustering/04_clustering_Syn.RData")
object <- UpdateSeuratObject(Sobject)

i = 0
for (module in colnames(genes_in_modules.fin)) {
  module = as.character(module)
  nmod = unlist(strsplit(module, split =  '.', fixed = T))[2]
  i = i+1
  # Select genes
  gene.set <- genes_in_modules.brown[,module]
  gene.set <- gene.set[is.na(gene.set) == F]
  gene.set = gene.set[gene.set %in% row.names(Sobject@assays[["RNA"]]@data) == TRUE]
  #Normalize
  norm.exp <- Sobject@assays[["RNA"]]@data[gene.set, ]/max(Sobject@assays[["RNA"]]@data[gene.set, ])
  max(norm.exp)
  # Get mean expression of genes of interest per cell
  mean.exp <- colMeans(x = as.matrix(norm.exp[gene.set, ]), na.rm = TRUE)
  # Add mean expression values in 'object@meta.data$gene.set.score'   DO THIS FOR EACH MODULE SI THAT IT CAN BE RETRIEVED LATER (assign)
  modname = module#paste(module, mods, sep = ".")
  Sobject@meta.data$gene.set.score <- mean.exp
  # Plot mean expression using Seurat::FeaturePlot()
  plota <- FeaturePlot(object = Sobject, features = "gene.set.score", cols = c("light grey","red"))
  #plotb <- plota + labs(title = modname)
  plotb <- plota + labs(title = paste (modname, sep = ' '))
  assign(paste("plot",i, sep="" ), plotb)
  rm(plota, plotb)
}
png('TASIC0/results/wgcna/finalmergeed/tsne1.png', width = 900, height = 1200)
cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5,  plot6, plot7, plot8, plot9, plot10, plot11, plot12, ncol = 3)
dev.off()

png('TASIC0/results/wgcna/finalmergeed/tsne2.png', width = 900, height = 1200)
cowplot::plot_grid(plot13, plot14, plot15, plot16, plot17,  plot18, plot19, plot20, plot21, plot22, plot23, plot24, ncol = 3)
dev.off()
