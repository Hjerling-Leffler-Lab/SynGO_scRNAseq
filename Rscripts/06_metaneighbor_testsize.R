##########################
#description: Test and visuailse differences between SynGO and random sets performance. Generates the figures in Figure 2 and Suppementary Figure 2 A,C
#author: Amparo Roig Adam
#last updated: 2022-09-13
#input: 
# - ranSynAUR, syngoAUR, ranAllAUR = AUROC scores of random synaptic gene sets, SynGO terms and random all-genes gene sets (from 05_metaneighbor_randomsets.R and 04_metaneighbor_synGOsets.R)
# - SynGO annotations and all SynGO genes
#output:
# - Scatter plot of all AUROC scores from randomly generated synaptic gene sets
# - Calculated empirical p values for each SynGO term compared to 10000 random synaptic gene sets of the same size
# - Scatter plot of AUROC vs gene set size relationship showing significant SynGO terms that perform different than random genes
##########################

library(reshape2)
library(ggplot2)
library(ggridges)
library(dplyr)
library(ggrepel)

# Read files with AUROC scores ---------------
#Random sets AUROC
ranSynAUR <- read.csv(file = "TASIC0/results/01_metaneighbor/20200228_bootstrap_subclass/200204_metaneighFast_subclass_randomsetAUROC.csv",
                      header = T, row.names = 1)
#SynGO AUROC
syngoAUR <- read.csv(file = "TASIC0/results/01_metaneighbor/20200228_bootstrap_subclass/200123_metaneighFast_subclass_synGOauroc.csv",
                     header = T, row.names = 1)
# Random sets all genes
ranAllAUR <- read.csv(file = "TASIC0/results/01_metaneighbor/200228_aurocNOsynRandomgenes_subclass.csv", 
                      header = T, row.names = 1)

#Remove too small groups (<5)
ranSynAUR <- na.omit(ranSynAUR[ranSynAUR[,1] >= 5,])
syngoAUR <- na.omit(syngoAUR[syngoAUR[,1] >= 5,])
ranAllAUR <- na.omit(ranAllAUR[ranAllAUR[,1] >= 5,])

rand <- melt(ranSynAUR, id.vars="X.1")
#rand$X.1 <- as.factor(rand$X.1)
syn <- melt(syngoAUR, id.vars="X.1")
allrand <- melt(ranAllAUR, id.vars="X.1")
#allrand$X.1 <- as.factor(allrand$X.1)

#Plot all random synaptic and SynGO  AUROCs----------
pdf(file = "TASIC0/figures/allscores.pdf", width = 6, height = 4,bg = "transparent")
ggplot() +
  geom_point(data = rand, aes(X.1, value, color = 'Random SynGO gene sets'), stroke=0,size = 1, alpha = 1) +
  #Different colours for different celltypes
  geom_point(data = syn, aes(X.1, value, color = 'SynGO categories'),stroke = 0, size = 2, alpha = 1) +
  scale_color_manual(values = c('Random SynGO gene sets' = 'black','SynGO categories' = 'red'))+
  #scale_fill_manual(values = c(NA, NA))+
  theme_minimal() + xlab("Gene set size (log10)") + ylab("AUROC score") + scale_x_log10()+
  theme(axis.line = element_line(size = 0.25), axis.ticks = element_line(size = 0.25), 
        text = element_text(size = 10))
dev.off()

#read calculated pvals
# Compute p values ----------
##pval = sum(s >= s0)/N 
  #sum(s >= s0) = number of times that the bootstrapped gene sets performed better than the SynGO gene sets
  # N = number of bootstrap samples
bootstrap.pval <- function(s0, len, N, celltype, bootstrap.values) {
  
  #Doublecheck that the number of random gene sets with a valid score (not NA) is equal to N
  if (nrow(bootstrap.values[bootstrap.values[,1] == len & !is.na(bootstrap.values[,2]),]) != N){
    stop ('Number of gene sets is not N')
  }
  pval = sum(bootstrap.values[bootstrap.values[,1] == len ,celltype] <= s0)/N #change direction of sign to do the opposite test - lower cell type discriminability than expected randomly
  return (as.numeric(pval))
}

## P value all cell types together
#Sum score in all cell types
ranAUR.sum <- as.data.frame(rowSums(ranSynAUR[,2:17]))
ranAUR.sum <- cbind(ranSynAUR[,1], rowSums(ranSynAUR[,2:17]))
colnames(ranAUR.sum) <- c("size", "sumAUR")

syngoAUR.sum <- as.data.frame(rowSums(syngoAUR[,2:17]))
colnames(syngoAUR.sum) <- "sumAUR"

#Empirical p-value calculation for SynGO terms vs random synaptic gene sets
test.pvals.sum <- sapply(row.names(syngoAUR.sum), function(go) 
{bootstrap.pval(s0 = syngoAUR.sum[go,], len = syngoAUR[go,1], N = 10000, celltype = "sumAUR", bootstrap.values = ranAUR.sum)})

#Write result table
write.csv(test.pvals.sum, file = "TASIC0/results/01_metaneighbor/lessthan_pval_bootstrap.csv", quote = F, row.names = T)

#write pval filtered for sunburst plot in SynGO platform
test.pvals.sum[test.pvals.sum$test.pvals.sum>0.05,]<-1
test.pvals.sum[test.pvals.sum$test.pvals.sum<0.05&test.pvals.sum$test.pvals.sum>0.01,]<-0.05
test.pvals.sum[test.pvals.sum$test.pvals.sum<0.01,]<-0.01
write.csv(test.pvals.sum, file = "TASIC0/results/01_metaneighbor/lessthan_pval_bootstrap_sunburstdata.csv", quote = F, row.names = T)

#Calculate mean AUROC values per group -------
mean_allrand <- allrand %>%
  group_by(X.1) %>%
  summarise(value = mean(value))
mean_rand <- rand %>%
  group_by(X.1) %>%
  summarise(value = mean(value))
mean_syn <- syn %>%
  group_by(X.1) %>%
  summarise(value = mean(value))
#Get all BP and CC annotations to generate independent plots
BP <- c()
for(i in BP_ontology){
  BP<- rbind(BP,as.data.frame(i)[,c('id', 'name')])
}
CC<-c()
for(i in CC_ontology){
  CC<- rbind(CC,as.data.frame(i)[,c('id', 'name')])
}
#
mean_syn <- as.data.frame(cbind(X.1 = syngoAUR$X.1,value = rowMeans(syngoAUR[,2:17])))
mean_syn <- mean_syn[row.names(mean_syn)%in%BP$id,]


#Plot results - mean AUROC values of SynGO terms+ color code by empirical p value; line of Mean values of randomly generated gene sets ------
pvals = data.frame(x = test.pvals.sum)

ggplot() + #to generate CC annotations as separate plots, change all occurrences of 'BP' below
  geom_line(data = mean_rand, aes(X.1, value, color = 'Random SynGO gene sets'), size = 0.9) +
  geom_line(data = mean_allrand, aes(X.1, value, color = 'Random gene sets \n expressed in the dataset'),size = 0.9) +
  theme_minimal() + xlab("Gene set size (log10)") + ylab("Mean AUROC score") + scale_x_log10() + 
  geom_point(data = mean_syn[row.names(mean_syn)%in%BP$id,], aes(X.1, value, 
        color = row.names(mean_syn[row.names(mean_syn)%in%BP$id,])%in%row.names(pvals)[pvals$x<0.05]), 
        size = 2, alpha = 1)+scale_color_manual(values = c('#619CFF','grey60' , 'black', 'red'))+ 
        theme(axis.line = element_line(size = 1), axis.text = element_text(size = 8, colour = 'black'),
        axis.title = element_text(size = 10), legend.position = c(0.7,0.25))

ggsave('TASIC0/figures/lessthanran_BP_labels.svg', width = 4, height = 4)

#Test random synaptic genes AUROC > random all genes AUROC
wilcox.test(x = mean_rand$value, y = mean_allrand$value, alternative = 'greater')






