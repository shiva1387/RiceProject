#####################################################
# Data Analysis in R-Metabolomics data Prof Prakash #
#####################################################
# Author(s): Shiv
# Version: 28072015
# Analysis inverse expression trends in genes and metabolites

##Libraries

library(VennDiagram)
library(d3heatmap) #interactive heatmaps
library(RColorBrewer)
library(dplyr)

#### Genes

##Reading in the files

setwd('../../Transcriptomics/data/')

WTEC<-read.csv("toptags_tt_WTEC_edgeR.csv",header=TRUE)
WTMU<-read.csv("toptags_tt_WTMU_edgeR.csv",header=TRUE)

## CONVERT TO NCBI GENE IDS
GeneDesc_GeneId<-read.table("../../Metabolomics/data/EnrichmentAnalysis/GeneDescrip_GeneId.txt",header=TRUE,check.names=FALSE,stringsAsFactors = FALSE,sep="\t")
GeneDesc_GeneId$KEGG_geneId<-gsub('osa:','',GeneDesc_GeneId$KEGG_geneId)
colnames(GeneDesc_GeneId)<-c('GeneId','KEGG_geneId')
GeneDesc_GeneId <- mutate_each(GeneDesc_GeneId, funs(toupper)) #Converting all names to upper case to match

##Creating a list of differential genes (fold change > or < 2 and p-values(FDR) <0.05)
WTEC_diff<-WTEC[WTEC$logFC >= 1 | WTEC$logFC <= -1 & WTEC$FDR <=0.05,c("X","logFC","FDR")]
colnames(WTEC_diff)<-c("GeneId","logFC_WTEC","FDR(p-value)_WTEC")
WTEC_diff<-merge(WTEC_diff,GeneDesc_GeneId,all=FALSE, by='GeneId', all.x= TRUE)
###NOTE: All logFC_WTEC regulation arewith respect to WT. As we are interested in the over-expression or mutant
### i have reversed the colors. Therefore, those 


WTEC_diff$color
# nrow(WTEC_diff)
# [1] 3705
nrow(WTEC_diff[WTEC_diff$logFC_WTEC>=0,])
#Up-regulated in WT 
#406
nrow(WTEC_diff[WTEC_diff$logFC_WTEC<=0,])
#Down-regulated in WT 
#2037
write.table(WTEC_diff,"WTEC_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

WTMU_diff<-WTMU[WTMU$logFC >= 1 | WTMU$logFC <= -1 & WTMU$FDR <=0.05,c("X","logFC","FDR")]
colnames(WTMU_diff)<-c("GeneId","logFC_WTMU","FDR(p-value)_WTMU")
# nrow(WTMU_diff)
# [1] 1210
nrow(WTMU_diff[WTMU_diff$logFC_WTMU>=0,])
#Up-regulated in WT 
#406
nrow(WTMU_diff[WTMU_diff$logFC_WTMU<=0,])
#Down-regulated in WT 
#804
# write.table(WTMU_diff,"WTMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

#### Plotting venn diagram of differential genes
venn.plot<-venn.diagram(list(as.character(WTEC_diff$GeneId),as.character(WTMU_diff$GeneId)), filename=NULL, margin=0.1, scaled=TRUE, fill=c("red", "green"), 
                        alpha=c(0.5,0.5), cex = 2, cat.dist=0.1, cat.fontface=4, category.names=c("WT vs EC", "WT vs MU"))
pdf("DiffentialGenes_pval2fc2.pdf")
grid.draw(venn.plot)
dev.off()

#Merging the 2 comparisons to look at trends

merge_WTECMU<-merge(WTEC_diff,WTMU_diff,all=FALSE)
rownames(merge_WTECMU)<-merge_WTECMU$GeneId
d3heatmap(merge_WTECMU[,c(2,4)], labRow=NULL, colors = "RdYlBu",scale="none", dendogram="both")

#Extracting Genes Up-reg in WT(EC) and Down-reg in WT(MU)
upWTEC_downWTMU<-merge_WTECMU[merge_WTECMU$logFC_WTEC>=0 & merge_WTECMU$logFC_WTMU<=0,]
nrow(upWTEC_downWTMU)
#25
#write.table(upWTEC_downWTMU,"upWTEC_downWTMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

#Extracting Genes Down-reg in WT(EC) and Up-reg in WT(MU)
downWTEC_upWTMU<-merge_WTECMU[merge_WTECMU$logFC_WTEC<=0 & merge_WTECMU$logFC_WTMU>=0,]
nrow(downWTEC_upWTMU)
#29
#write.table(downWTEC_upWTMU,"downWTEC_upWTMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

#Extracting Genes Up-reg in WT(EC) and Up-reg in WT(MU)
upWTEC_upWTMU<-merge_WTECMU[merge_WTECMU$logFC_WTEC>=0 & merge_WTECMU$logFC_WTMU>=0,]
nrow(upWTEC_upWTMU)
#253
#write.table(upWTEC_upWTMU,"upWTEC_upWTMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

#Extracting Genes Down-reg in WT(EC) and Down-reg in WT(MU)
downWTEC_downWTMU<-merge_WTECMU[merge_WTECMU$logFC_WTEC<=0 & merge_WTECMU$logFC_WTMU<=0,]
nrow(downWTEC_downWTMU)
#600
write.table(downWTEC_downWTMU,"downWTEC_downWTMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

##Genes differential only in WTEC
onlyWTEC_notWTMU<-anti_join(WTEC_diff, WTMU_diff, by="GeneId")
nrow(onlyWTEC_notWTMU)
#2798
write.table(onlyWTEC_notWTMU,"onlyWTEC_notWTMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

##Genes differential only in WTMU
onlyWTMU_notWTEC<-anti_join(WTMU_diff, WTEC_diff, by="GeneId")
nrow(onlyWTMU_notWTEC)
#303
write.table(onlyWTMU_notWTEC,"onlyWTMU_notWTEC_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)


