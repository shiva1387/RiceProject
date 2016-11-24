#####################################################
# Data Analysis in R-Metabolomics data Prof Prakash #
#####################################################
# Author(s): Shiv
# Version: 11102016
# Categorizing plant gsea results

##Libraries

library(VennDiagram)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(data.table)
library(pheatmap)
library(ggplot2)
library(preprocessCore)
library(gplots)
library(reshape2)

#### Genes

##Reading in the files
rm(list=ls())
setwd('/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/ricemanuscript/Data/PlantGsea_101016/')

########## Fuctions

normalizeData<-function(data_matrix){
  data_matrix<-as.matrix(log1p(data_matrix))
  processed_data<-normalize.quantiles(as.matrix(data_matrix),copy=TRUE)
  colnames(processed_data)<-colnames(data_matrix)
  rownames(processed_data)<-rownames(data_matrix)
  return(processed_data)
}
gg_color_hue <- function(n) {
  hues = seq(10, 380, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

CountUpDownReg<-function(x,datasetName) #Returns number of up-rgulated, down-regulated and log2ratio of number of up-reg/down-reg genes
{
 #x is a vector of loc id  
 #
 x<-x[!is.na(x)]
 DiffgenesInCategory<-datasetName[datasetName$LOC %in% x,]
 DiffgenesInCategory_v1<-DiffgenesInCategory[,c(grep("LOC|logFC_WT",colnames(DiffgenesInCategory)))]
 #DiffgenesInCategory_v1<-DiffgenesInCategory_v1 %>% group_by (LOC) %>% summarise(logFc=mean(logFC_WTEC))
 DiffgenesInCategory_v1<-DiffgenesInCategory_v1 %>% group_by (LOC) %>% summarise(logFc=mean(logFC_WTMU))
 DiffgenesInCategory_v2<-DiffgenesInCategory_v1 %>% summarise(upReginPerturbation=length(which(logFc<0)),
                                                              downReginPerturbation=length(which(logFc>0)))
 RegulationRatio<-log2(DiffgenesInCategory_v2$upReginPerturbation/DiffgenesInCategory_v2$downReginPerturbation) #If more than 50% of genes are up-regulated, then log2 will be positive, else negative
 names(RegulationRatio)<-"RegulationRatio"
 #Getting only the regulation
 return(list(unlist(DiffgenesInCategory_v2),RegulationRatio))
}
 
CategorySpecificRatio_WTEC<-function(data1) {
  #testData,WTEC_diffGeneTable_loc
  df<-apply(data1,1,function(x) {
  #datasetName1<-substitute(datasetName)  
  GenesinCategory<-unlist(strsplit(as.character(x["GeneList"]),"\\s"))
  tryCatch(CountUpDownReg(GenesinCategory,WTEC_diffGeneTable_loc), error=function(x) NA)
  })
  return(df)
}

CategorySpecificRatio_WTMU<-function(data1) {
  #testData,WTEC_diffGeneTable_loc
  df<-apply(data1,1,function(x) {
    #datasetName1<-substitute(datasetName)  
    GenesinCategory<-unlist(strsplit(as.character(x["GeneList"]),"\\s"))
    tryCatch(CountUpDownReg(GenesinCategory,WTMU_diffGeneTable_loc), error=function(x) NA)
  })
  return(df)
}

Conversion.list.df<-function(CategorySpecificRatioOutput,dataset) {
  output1<-data.frame(matrix(unlist(CategorySpecificRatioOutput), nrow=nrow(dataset), byrow=T),
                      stringsAsFactors=FALSE)
  return(output1)
}

#Calculate from the genes in overlap, how many are up-reg (% differential genes in that category) and how many down-reg (% differential genes in that category)
#Category represents what % of differential genes
#Category size, number of genes in overlap

######################Reading in the transcriptomics data

diffGeneTable<-read.table("/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/Transcriptomics/data/differentialGeneTable.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
diffGeneTable$GeneId<-rownames(diffGeneTable)
diffGeneTable_V1<-data.matrix(cbind(diffGeneTable$GeneId,normalizeData(as.matrix(diffGeneTable[,c(1:6)]))))
colnames(diffGeneTable_V1)[1]<-"GeneId"

WTEC<-read.csv("/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject//Transcriptomics/data/toptags_tt_WTEC_edgeR.csv",header=TRUE)
WTMU<-read.csv("/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/Transcriptomics/data/toptags_tt_WTMU_edgeR.csv",header=TRUE)
rownames(WTEC)<-WTEC$X
rownames(WTMU)<-WTMU$X
## CONVERT TO NCBI GENE IDS
GeneDesc_GeneId<-read.table("../../../Metabolomics/data/EnrichmentAnalysis/GeneDescrip_GeneId.txt",header=TRUE,check.names=FALSE,stringsAsFactors = FALSE,sep="\t")
GeneDesc_GeneId$KEGG_geneId<-gsub('osa:','',GeneDesc_GeneId$KEGG_geneId)
colnames(GeneDesc_GeneId)<-c('GeneId','KEGG_geneId')
GeneDesc_GeneId <- mutate_each(GeneDesc_GeneId, funs(toupper)) #Converting all names to upper case to match

diffGeneTable_Kegg<-merge(diffGeneTable,GeneDesc_GeneId,all=FALSE, by='GeneId', all.x= TRUE)
diffGeneTable_Kegg<-diffGeneTable_Kegg[complete.cases(diffGeneTable_Kegg$KEGG_geneId),]
write.table(diffGeneTable_Kegg,"GeneExpressionTable_withKegg.txt",sep="\t",quote=FALSE,col.names=NA)

## Reading a list of locus ids
LocId<-read.table("../../../Transcriptomics/GenId_LOCLink.txt",header=TRUE,check.names=FALSE,stringsAsFactors = FALSE,sep="\t")
colnames(LocId)<-c("GeneId","LOC")
LocId$GeneId <- toupper(LocId$GeneId) #Converting geneId to upper case

##Creating a list of differential genes (fold change > or < 2 and p-values(FDR) <0.05)
#WTEC_diff<-WTEC[WTEC$logFC >= 1 | WTEC$logFC <= -1 & WTEC$FDR <=0.05,c("X","logFC","FDR")] 
WTEC_diff<-WTEC[(WTEC$logFC >= 1 | WTEC$logFC <= -1) & WTEC$FDR <=0.05,c("X","logFC","FDR")] #Modified on 11-Aug 2016 as the previous code will select genes if it satisfies either condition!
colnames(WTEC_diff)<-c("GeneId","logFC_WTEC","FDR(p-value)_WTEC")
# nrow(WTEC_diff)
# [1] 3705
# [1] 3669 #Modified on 11-Aug 2016 as the previous code will select genes if it satisfies either condition! #There are now 36 genes less

nrow(WTEC_diff[WTEC_diff$logFC_WTEC>=0,])
#Up-regulated in EC 
#1632 #11-Aug
nrow(WTEC_diff[WTEC_diff$logFC_WTEC<=0,])
#Down-regulated in EC 
#2037
WTEC_diff<-merge(WTEC_diff,GeneDesc_GeneId,all=FALSE, by='GeneId', all.x= TRUE)
###NOTE: All logFC_WTEC regulation arewith respect to WT. As we are interested in the over-expression or mutant
WTEC_diffGeneTable<-merge(WTEC_diff,diffGeneTable,all=FALSE, by='GeneId', all.x= TRUE)
##Adding gene id, note, one gene id may may to multiple locus id
WTEC_diffGeneTable_loc<-merge(WTEC_diffGeneTable,LocId,all=FALSE, by='GeneId', all.x= TRUE) 
# write.table(WTEC_diffGeneTable_loc,"WTEC_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

WTMU_diff<-WTMU[(WTMU$logFC >= 1 | WTMU$logFC <= -1) & WTMU$FDR <=0.05,c("X","logFC","FDR")]
colnames(WTMU_diff)<-c("GeneId","logFC_WTMU","FDR(p-value)_WTMU")
# nrow(WTMU_diff)
#1089 #11-Aug
nrow(WTMU_diff[WTMU_diff$logFC_WTMU>=0,])
#Up-regulated in MU 
#285
nrow(WTMU_diff[WTMU_diff$logFC_WTMU<=0,])
#Down-regulated in MU 
#804
WTMU_diff<-merge(WTMU_diff,GeneDesc_GeneId,all=FALSE, by='GeneId', all.x= TRUE)
WTMU_diffGeneTable<-merge(WTMU_diff,diffGeneTable,all=FALSE, by='GeneId', all.x= TRUE) 

##Adding gene id, note, one gene id may may to multiple locus id
WTMU_diffGeneTable_loc<-merge(WTMU_diffGeneTable,LocId,all=FALSE, by='GeneId', all.x= TRUE) 

#### Analyzing categories
########### WT EC
Gsea_output<-read.table("PLANTGSEA_MotifGeneSets_WTMU_101016.txt",header = FALSE,sep="\t",stringsAsFactors = FALSE, fill = NA)
colnames(Gsea_output)<-c("GeneSetName","NumberOfGenes","Class","GenesinOverlap","pValue","FDR","GeneList") #"Category"

listWithRatio<-CategorySpecificRatio_WTMU(Gsea_output)
dfWithRatio<-Conversion.list.df(listWithRatio,Gsea_output)
colnames(dfWithRatio)<-c("NumUpRegulatedGenes","NumDownRegulatedGenes","log2(Ratio_UpDown)")
Gsea_output <- cbind(Gsea_output,dfWithRatio)
Gsea_output_sig<-Gsea_output[Gsea_output$FDR<=0.05,]

write.table(Gsea_output_sig,"PLANTGSEA_MotifGeneSets_WTMU_101016_sig.txt",sep="\t",quote=FALSE,col.names=NA)
  
