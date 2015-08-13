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
library(tidyr)

#### Genes

##Reading in the files

setwd('../../Transcriptomics/data/')

WTEC<-read.csv("toptags_tt_WTEC_edgeR.csv",header=TRUE)
WTMU<-read.csv("toptags_tt_WTMU_edgeR.csv",header=TRUE)
rownames(WTEC)<-WTEC$X
rownames(WTMU)<-WTMU$X
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
WTEC_diff<- mutate(WTEC_diff, KEGGcolor = ifelse(logFC_WTEC <= 0 , "red","blue")) 
#by default, it should be red for up-regulation and blue for down, however here the 
#colors are based on expression patterns in overexpression, and hence, reversed
# nrow(WTEC_diff)
# [1] 3705
nrow(WTEC_diff[WTEC_diff$logFC_WTEC>=0,])
#Up-regulated in WT 
#1668
nrow(WTEC_diff[WTEC_diff$logFC_WTEC<=0,])
#Down-regulated in WT 
#2037
# write.table(WTEC_diff,"WTEC_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

WTMU_diff<-WTMU[WTMU$logFC >= 1 | WTMU$logFC <= -1 & WTMU$FDR <=0.05,c("X","logFC","FDR")]
colnames(WTMU_diff)<-c("GeneId","logFC_WTMU","FDR(p-value)_WTMU")
WTMU_diff<-merge(WTMU_diff,GeneDesc_GeneId,all=FALSE, by='GeneId', all.x= TRUE)
WTMU_diff<- mutate(WTMU_diff, KEGGcolor = ifelse(logFC_WTMU <= 0 , "red","blue")) 
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
merge_WTECMU<-merge(WTEC_diff,WTMU_diff,by='GeneId', all=FALSE)
rownames(merge_WTECMU)<-merge_WTECMU$GeneId
write.table(merge_WTECMU,"Common_WTEC-WTMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)
#d3heatmap(merge_WTECMU[,c(2,4)], labRow=NULL, colors = "RdYlBu",scale="none", dendogram="both")

#Extracting Genes Up-reg in EC and Down-reg in MU
upWTEC_downWTMU<-merge_WTECMU[merge_WTECMU$logFC_WTEC<=0 & merge_WTECMU$logFC_WTMU>=0,]
nrow(upWTEC_downWTMU)
#29
#write.table(upWTEC_downWTMU,"upEC_downMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

#Extracting Genes Down-reg in EC and Up-reg in MU
downWTEC_upWTMU<-merge_WTECMU[merge_WTECMU$logFC_WTEC>=0 & merge_WTECMU$logFC_WTMU<=0,]
nrow(downWTEC_upWTMU)
#25
#write.table(downWTEC_upWTMU,"downEC_upMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

#Extracting Genes Up-reg in EC and Up-reg in MU
upWTEC_upWTMU<-merge_WTECMU[merge_WTECMU$logFC_WTEC<=0 & merge_WTECMU$logFC_WTMU<=0,]
nrow(upWTEC_upWTMU)
#600
#write.table(upWTEC_upWTMU,"upEC_upMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

#Extracting Genes Down-reg in EC and Down-reg in MU
downWTEC_downWTMU<-merge_WTECMU[merge_WTECMU$logFC_WTEC>=0 & merge_WTECMU$logFC_WTMU>=0,]
nrow(downWTEC_downWTMU)
#253
write.table(downWTEC_downWTMU,"downEC_downMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

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

#### Formatting gene id and pathway membership table
#For EC
WTEC_genePathwayInfo<-read.table("../../Metabolomics/data/EnrichmentAnalysis/gFC2mFC2pv05/WTEC_genePathway.txt",header=TRUE,
                                 sep="\t")

WTEC_geneRegulation<-WTEC_genePathwayInfo[,c(3,4)]
##Splitting the annMolecule list column which contains gene id's separated by ';'
##this list will contain only genes that have been assigned to a pathway in KEGG
WTEC_genePathwayInfo_long<-setDT(WTEC_genePathwayInfo)[, strsplit(as.character(annMoleculeList),";",fixed=TRUE), by=list(pathwayName,annMoleculeList)][,.(pathwayName,annMoleculeList=V1)]
setnames(WTEC_genePathwayInfo_long,"annMoleculeList","GeneId")
WTEC_genePathwayInfo_long<-as.data.frame(WTEC_genePathwayInfo_long)
WTEC_genePathwayInfo_longReg<-merge(WTEC_genePathwayInfo_long,WTEC_geneRegulation,by='GeneId', all.x= TRUE)

#Count and summarise using dplyr

WTEC_genePathwayInfo_longReg_count<-WTEC_genePathwayInfo_longReg %>% 
  group_by (pathwayName, Regulation) %>% 
  tally(sort=TRUE)
write.table(WTEC_genePathwayInfo_longReg_count,"WTEC_genePathwayInfo_longReg_count.txt",sep="\t",quote=FALSE,row.names=FALSE)

WTEC_genePathwayInfo_longReg_count$sample<-rep("EC",nrow(WTEC_genePathwayInfo_longReg_count))

#For MU
WTMU_genePathwayInfo<-read.table("../../Metabolomics/data/EnrichmentAnalysis/gFC2mFC2pv05/WTMU_genePathway.txt",header=TRUE,
                                 sep="\t")

WTMU_geneRegulation<-WTMU_genePathwayInfo[,c(3,4)]
##Splitting the annMolecule list column which contains gene id's separated by ';'
##this list will contain only genes that have been assigned to a pathway in KEGG
WTMU_genePathwayInfo_long<-setDT(WTMU_genePathwayInfo)[, strsplit(as.character(annMoleculeList),";",fixed=TRUE), by=list(pathwayName,annMoleculeList)][,.(pathwayName,annMoleculeList=V1)]
setnames(WTMU_genePathwayInfo_long,"annMoleculeList","GeneId")
WTMU_genePathwayInfo_long<-as.data.frame(WTMU_genePathwayInfo_long)
WTMU_genePathwayInfo_longReg<-merge(WTMU_genePathwayInfo_long,WTMU_geneRegulation,by='GeneId', all.x= TRUE)

#Count and summarise using dplyr

WTMU_genePathwayInfo_longReg_count<-WTMU_genePathwayInfo_longReg %>% 
  group_by (pathwayName, Regulation) %>% 
  tally(sort=TRUE)

write.table(WTMU_genePathwayInfo_longReg_count,"WTMU_genePathwayInfo_longReg_count.txt",sep="\t",quote=FALSE,row.names=FALSE)

WTMU_genePathwayInfo_longReg_count$sample<-rep("MU",nrow(WTMU_genePathwayInfo_longReg_count))

WTECMU_pathway_trends<-rbind(WTEC_genePathwayInfo_longReg_count,WTMU_genePathwayInfo_longReg_count)

WTECMU_pathway_trends1<-WTECMU_pathway_trends %>% 
  group_by(pathwayName,sample) %>% 
  mutate(EC=cumsum(ifelse(sample=="EC",n,0))) %>% 
  mutate(MU=cumsum(ifelse(sample=="MU",n,0)))



##Plotting pathway trends
WTECMU_pathway_trends$Regulation<-factor(WTECMU_pathway_trends$Regulation, levels=c('up','down'))
  
plot1<-ggplot(data=WTECMU_pathway_trends, aes(x=pathwayName, y=n, fill= factor(Regulation))) +
  geom_bar(stat="identity", position="dodge") + facet_grid(sample~.) +
  theme_bw() + theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=12),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     #panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black')) + ggtitle("Pathway trends") + ylab("No of Diff genes")

ggsave("PathwayTrends.pdf",plot1)

#####Making long to wide

WTECMU_pathway_trends_wide <-dcast(WTECMU_pathway_trends, pathwayName + Regulation ~ sample, value.var="n")
write.table(WTECMU_pathway_trends_wide,"WTECMU_pathway_trends.txt",sep="\t",quote=FALSE,row.names=FALSE)
