##########################
#description: SynGO enrichment of WGCNA gene modules
#author: Amparo Roig Adam
#last updated: 2022-01-25
#input: 
# - table of WGCNA modules from  04_wgcna.R
# - NeuronalData and NeuronalMetadata = Raw count data previously filtered in 00_read+filter_NeuronalData
# - SynGO gene annotations
#output:
# - SynGO enrichment of each individual WGCNA gene module
# - combinedSynGOenrichment.txt for each cluster of WGCNA modules
##########################

library("Category")

#Data
#Modules wgcna - after merging similar modules
load("~/Desktop/SynGO/TASIC0/results/wgcna/finalmergeed/20220125_final_merged_modules.RData")

# GO term enrichment: Hypergeometrical test
#read syngo
BP_annotations <- fromJSON(file = "data/SynGO_genesets_2018-6-21/BP_annotations.json")
BP_ontology <- fromJSON(file = "data/SynGO_genesets_2018-6-21/BP_ontology.json")
CC_annotations <- fromJSON(file = "data/SynGO_genesets_2018-6-21/CC_annotations.json")
CC_ontology <- fromJSON(file = "data/SynGO_genesets_2018-6-21/CC_ontology.json")

#Background data -------
## SynGO terms - list of GOs that contain a list of characters(genes)
SynGOsets <- c(BP_annotations$aggregate$hgnc_symbol, CC_annotations$aggregate$hgnc_symbol)
SynGOsets <- Filter(function(x) length(x) >= 1, SynGOsets) #eliminate empty GO terms
#Ontology
SynGO_ontology <- c(BP_ontology, CC_ontology)
rm(BP_ontology, CC_ontology)

#Filter out genes from GO terms that are not in the dataset (NeuronalData)
load("TASIC0/data/00_Neuronal_Data.RData")

SynGOsetsDB <- SynGOsets 
for (setnumber in 1:length(SynGOsetsDB)){
  SynGOsetsDB[[setnumber]] <- Filter(function(x) x%in%row.names(NeuronalData) == T, SynGOsetsDB[[setnumber]]) #NeuronalData needs to be changed to HumanNeuronalData if running human dataset
}  
rm(NeuronalData, NeuronalMetadata)

#Universe: all genes in SynGO detected in the dataset
universe <- unique(unlist(SynGOsetsDB, use.names= F ))

#write(universe, file = "universe.txt", sep = ' ')

#genes in the module  #FROM HERE IT SHOULD BE DONE FOR ALL MODULES/all selected modules.. ---------

i = 0
for (module in colnames(genes_in_modules.merged)){###
  i = i+1
  module = as.character(module)
  nmod = unlist(strsplit(module, split =  '.', fixed = T))[2]
  
  siggenes <- genes_in_modules.merged[,module]
  siggenes <- siggenes[is.na(siggenes) == F] #remove NAs
  
  #GOs to which the genes in the module belong
  sigsets <- Map(function(x, y) x[x %in% y], SynGOsetsDB, MoreArgs=list(y=siggenes))
  
  #HYPERGEOMETRIC TEST--------
  set.seed(123) 
  go_enrich <- hyperg(SynGOsetsDB, sigsets, universe)
  #Filter results
  go_enrich <- go_enrich[go_enrich$significant != 0,]
  
  #Order by p-val
  go_enrich <- go_enrich[order(unlist(go_enrich$p)),]
  #Add corrected p-val
  correctedPval <- p.adjust(go_enrich$p, method = 'bonferroni')
  go_enrich <- cbind(go_enrich, correctedPval)
  
  #Add column with Domain (BP or CC) to go_enrich--------------
  domain <- c()
  go_name <- c()
  
  for(go in row.names(go_enrich)){   
    if (go %in% names(BP_annotations$aggregate$hgnc_symbol) == T){
      domain <- c(domain, "BP")
    }
    if (go %in% names(CC_annotations$aggregate$hgnc_symbol) == T){
      domain <- c(domain, "CC")
    }
    go_name <- c(go_name,SynGO_ontology[[which(lapply( SynGO_ontology, function(x) x$id) %in% go)]]$name)
    rm(go)}
  go_enrich <- cbind(go_enrich, domain, go_name)
  
  #Generate Output file for SynGO sunburst colour-coding tool ------------
  tofile <- as.data.frame(go_enrich[,'correctedPval'],row.names = row.names(go_enrich))
  colnames(tofile)<- 'correctedPval'
  write.csv(tofile,file = paste('TASIC0/results/wgcna/finalmergeed/GOenrichment/',module,'.txt', sep = ''),      ###  
            quote = F, sep = ",", row.names = T, col.names = F)
   
 table <-  go_enrich[1:10,c('correctedPval', 'go_name')]
 write.csv(table, file = paste('TASIC0/results/wgcna/finalmergeed/GOenrichment/',module, '.csv', sep =''),sep = ',', quote = F) 
}


##############
#Summarised enrichment for each module type as classified in 09_gene_module_clustering.R
#module classification results
load('TASIC0/results/wgcna/finalmergeed/TASIC0/results/wgcna/finalmergeed/module_classification.RData')
genmods.kmeans$cluster

#for each type
for(type in c(1:3)){
  sumGO <- data.frame()
  print(type)
  #for each module in the module group collect syngo enrichment
  for(module in names(genmods.kmeans$cluster[genmods.kmeans$cluster==type])){
    mod.enrich <- read.table(paste('TASIC0/results/wgcna/finalmergeed/GOenrichment/',module,'.txt',sep = ''), header = T, sep = ',')
    colnames(mod.enrich)<- c('GOid', 'correctedPval')
    for(GO in mod.enrich$GOid){
      #add new row if GO has not appeared before
      if(GO%in%sumGO == F){
        #score 0 if it was not significant and 1 if significant
        if(mod.enrich[mod.enrich$GOid==GO, 'correctedPval']<0.05){
          sumGO <- rbind(sumGO, cbind(GO = GO, n = 1))
        }else{
          #print('else1')
          sumGO <- rbind(sumGO, cbind(GO = GO, n = 0.1))} #score can't be 0 when visualising on syngo portal
      }
      #if it is already in the sumGO dataframe, +1 score if significant
      if(GO%in%sumGO == T){
        print('else2')
        if(mod.enrich[mod.enrich$GOid==GO, 'correctedPval']<0.05){
          sumGO[sumGO$GO == GO,]$n <- as.numeric(sumGO[sumGO$GO == GO,]$n)+1
        }}}}
  #save final table as input for sunburst colour coding tool in syngo.com
  write.table(sumGO, file = paste('TASIC0/results/wgcna/finalmergeed/GOenrichment/modulesckuster',type,'_combinedSynGOenrichment.txt', sep = ''), col.names = F, row.names = F, quote = F, sep = ',')
}
