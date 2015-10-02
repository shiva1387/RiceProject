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
library(data.table)
library(pheatmap)
library(ggplot2)
library(preprocessCore)
library(gplots)

#### Genes

##Reading in the files

setwd('../../Transcriptomics/data/')

########## Fuctions

normalizeData<-function(data_matrix){
  data_matrix<-as.matrix(log1p(data_matrix))
  processed_data<-normalize.quantiles(as.matrix(data_matrix),copy=TRUE)
  colnames(processed_data)<-colnames(data_matrix)
  rownames(processed_data)<-rownames(data_matrix)
  return(processed_data)
}
###########

diffGeneTable<-read.table("differentialGeneTable.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
diffGeneTable$GeneId<-rownames(diffGeneTable)

## Reading a list pathway memberships for each gene
GenePathway<-read.table("osa_pathway.list",header=TRUE,check.names=FALSE,stringsAsFactors = FALSE,sep="\t")
colnames(GenePathway)<-c("KEGG_geneId","PathwayId")
GenePathway$KEGG_geneId <- gsub("osa:","",GenePathway$KEGG_geneId) #Removing osa:
GenePathway$PathwayId <- gsub("osa","",GenePathway$PathwayId) #Removing osa:

###################### Plotting PCA for differential genes

#### PCA

dataset<-normalizeData(as.matrix(diffGeneTable[,c(1:6)])) #ignoring GeneId column
pca_results <- princomp(dataset,cor=F,scores=T) ### IMP: choose quantile normalized or scaled data

residual_variance<-pca_results$sdev^2/sum(pca_results$sdev^2)

plot(pca_results$loadings[,2]~pca_results$loadings[,1])
text(pca_results$loadings[,2]~pca_results$loadings[,1], labels = colnames(dataset), cex=0.6, pos=4)

SampleGroups<-sapply(colnames(dataset), function(x) strsplit(x,"\\.")[[1]][1])
SampleGroups<-as.vector(SampleGroups)


### PCA ggplot
forPlot<-data.frame(PCaxisA = pca_results$loadings[,1],PCaxisB = pca_results$loadings[,2], 
                    dataType=SampleGroups) #Subset SampleGroups as the last 9 are blanks

##plotting
plot1<- ggplot(data=forPlot, aes(x=PCaxisA, y=PCaxisB, colour= factor(SampleGroups), shape= factor(SampleGroups))) + geom_point(size=4) #dataType
plot2<- plot1 +  #geom_text(data = forPlot, aes(x=PCaxisA, y=PCaxisB, label= colnames(dataForPlotting)), hjust = 1) +
  theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1,size=12),axis.text.y=element_text(size=12),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black'))+ 
  xlab(paste0("PC 1 loadings","\n","Variation exp= ",round(residual_variance[1]*100,2),"%")) + 
  ylab(paste0("PC 2 loadings","\n","Variation exp= ",round(residual_variance[2]*100,2),"%")) +
  ggtitle("PCA RNA-seq")
# - (metabolomics data for samples with rnaseq + without outliers)

pdf("pca_pc12_rnaseq.pdf",height=12,width=12)
plot2
dev.off()

######################

WTEC<-read.csv("toptags_tt_WTEC_edgeR.csv",header=TRUE)
WTMU<-read.csv("toptags_tt_WTMU_edgeR.csv",header=TRUE)
rownames(WTEC)<-WTEC$X
rownames(WTMU)<-WTMU$X
## CONVERT TO NCBI GENE IDS
GeneDesc_GeneId<-read.table("../../Metabolomics/data/EnrichmentAnalysis/GeneDescrip_GeneId.txt",header=TRUE,check.names=FALSE,stringsAsFactors = FALSE,sep="\t")
GeneDesc_GeneId$KEGG_geneId<-gsub('osa:','',GeneDesc_GeneId$KEGG_geneId)
colnames(GeneDesc_GeneId)<-c('GeneId','KEGG_geneId')
GeneDesc_GeneId <- mutate_each(GeneDesc_GeneId, funs(toupper)) #Converting all names to upper case to match

## Reading a list of locus ids
LocId<-read.table("../GenId_LOCLink.txt",header=TRUE,check.names=FALSE,stringsAsFactors = FALSE,sep="\t")
colnames(LocId)<-c("GeneId","LOC")
LocId$GeneId <- toupper(LocId$GeneId) #Converting geneId to upper case

##Creating a list of differential genes (fold change > or < 2 and p-values(FDR) <0.05)
WTEC_diff<-WTEC[WTEC$logFC >= 1 | WTEC$logFC <= -1 & WTEC$FDR <=0.05,c("X","logFC","FDR")]
colnames(WTEC_diff)<-c("GeneId","logFC_WTEC","FDR(p-value)_WTEC")
# nrow(WTEC_diff)
# [1] 3705
nrow(WTEC_diff[WTEC_diff$logFC_WTEC>=0,])
#Up-regulated in WT 
#1668
nrow(WTEC_diff[WTEC_diff$logFC_WTEC<=0,])
#Down-regulated in WT 
#2037
WTEC_diff<-merge(WTEC_diff,GeneDesc_GeneId,all=FALSE, by='GeneId', all.x= TRUE)
###NOTE: All logFC_WTEC regulation arewith respect to WT. As we are interested in the over-expression or mutant
### i have reversed the colors. Therefore, those 
WTEC_diff<- mutate(WTEC_diff, KEGGcolor = ifelse(logFC_WTEC <= 0 , "red","blue")) 
#by default, it should be red for up-regulation and blue for down, however here the 
#colors are based on expression patterns in overexpression, and hence, reversed
WTEC_diffGeneTable<-merge(WTEC_diff,diffGeneTable,all=FALSE, by='GeneId', all.x= TRUE)
##Adding gene id, note, one gene id may may to multiple locus id
WTEC_diffGeneTable_loc<-merge(WTEC_diffGeneTable,LocId,all=FALSE, by='GeneId', all.x= TRUE) 
# write.table(WTEC_diffGeneTable_loc,"WTEC_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

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
WTMU_diff<-merge(WTMU_diff,GeneDesc_GeneId,all=FALSE, by='GeneId', all.x= TRUE)
WTMU_diff<- mutate(WTMU_diff, KEGGcolor = ifelse(logFC_WTMU <= 0 , "red","blue")) 
WTMU_diffGeneTable<-merge(WTMU_diff,diffGeneTable,all=FALSE, by='GeneId', all.x= TRUE) 

##Adding gene id, note, one gene id may may to multiple locus id
WTMU_diffGeneTable_loc<-merge(WTMU_diffGeneTable,LocId,all=FALSE, by='GeneId', all.x= TRUE) 
# write.table(WTMU_diffGeneTable_loc,"WTMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

#pdf("WTMU_diffGenes.pdf",onefile=FALSE)
pheatmap(temp, scale="row",cluster_rows = TRUE, cluster_cols = FALSE,#labels_row = as.character(WTMU_diffGeneTable[,1]),
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=8)
#dev.off()


#### Plotting venn diagram of differential genes
venn.plot<-venn.diagram(list(as.character(WTEC_diff$GeneId),as.character(WTMU_diff$GeneId)), filename=NULL, margin=0.1, scaled=TRUE, fill=c("red", "green"), 
                        alpha=c(0.5,0.5), cex = 2, cat.dist=0.1, cat.fontface=4, category.names=c("WT vs EC", "WT vs MU"))
pdf("DiffentialGenes_pval2fc2.pdf")
grid.draw(venn.plot)
dev.off()

#Merging the 2 comparisons to look at trends
merge_WTECMU<-merge(WTEC_diff,WTMU_diff,by='GeneId', all=FALSE)
rownames(merge_WTECMU)<-merge_WTECMU$GeneId
#write.table(merge_WTECMU,"Common_WTEC-WTMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)
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

##Getting genes which have inverse trends
InverseTrends<-rbind(upWTEC_downWTMU,downWTEC_upWTMU)
InverseTrends_loc<-merge(InverseTrends,LocId,all=FALSE, by='GeneId', all.x= TRUE) 
InverseTrends_diffGeneTable_loc<-merge(InverseTrends_loc[,c(1,4)],diffGeneTable,all=FALSE, by='GeneId', all.x= TRUE) 
rownames(InverseTrends_diffGeneTable_loc)<-as.character(InverseTrends_diffGeneTable_loc[,1])

pdf("InverseTrends.pdf",onefile=FALSE)
pheatmap(InverseTrends_diffGeneTable_loc[,c(3:8)], scale="row",cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=8)
dev.off()

# write.table(InverseTrends_loc,"InverseTrends_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

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

#### Formatting gene id and pathway membership table #obtained from gene_enrichment results using isubpathway_miner.R
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

##########################################################################################
##################### Plotting Heatmaps according to pathway memberships #################
##### The pathway mmberships are in a similar concept as above
##### The only difference being that the pathways are grouped in different categories according to PathwayClassificationHeatmaps.xlsx
##### All the genes in that particular category is then extracted and plotted as heatmaps

##Getting amino acid  associated genes
aminoacids<-c("path:00400", "path:00270","path:00260","path:00280","path:00330","path:00360","path:00380",
"path:00350","path:00250","path:00410","path:00460")
GenePathway_aminoacids<-GenePathway[GenePathway$PathwayId %in% aminoacids,]
Genes_aminoacids<-unique(as.numeric(GenePathway_aminoacids$KEGG_geneId))

#EC
WTEC_diffGeneTable_aa<-merge(WTEC_diffGeneTable,GenePathway_aminoacids,all=FALSE, by='KEGG_geneId') 
pdf("WTEC_aa.pdf",onefile=FALSE)
pheatmap(WTEC_diffGeneTable_aa[,c(6:11)], scale="row", labels_row= as.character(WTEC_diffGeneTable_aa$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs EC Amino acids ")
dev.off()

#MU
WTMU_diffGeneTable_aa<-merge(WTMU_diffGeneTable,GenePathway_aminoacids,all=FALSE, by='KEGG_geneId') 
pdf("WTMU_aa.pdf",onefile=FALSE)
pheatmap(WTMU_diffGeneTable_aa[,c(6:11)], scale="row", labels_row= as.character(WTMU_diffGeneTable_aa$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs MU Amino acids ")
dev.off()

#ECMU
WTECMU_diffGeneTable_aa_GeneID<-unique(c(as.character(WTEC_diffGeneTable_aa$GeneId),as.character(WTMU_diffGeneTable_aa$GeneId)))
WTECMU_diffGeneTable_aa<-diffGeneTable[diffGeneTable$GeneId%in% WTECMU_diffGeneTable_aa_GeneID,]

pdf("WTECMU_aa.pdf",onefile=FALSE)
pheatmap(WTECMU_diffGeneTable_aa[,c(1:6)], scale="row", cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs EC or MU, Amino acids ")
dev.off()

##Getting fatty acid  associated genes
fattyacids<-c("path:01212,path:00592,path:00590,path:01040,path:00061,path:00071,path:00591,path:00062")
fattyacids<-unlist(strsplit(fattyacids,","))
GenePathway_fattyacids<-GenePathway[GenePathway$PathwayId %in% fattyacids,]
Genes_fattyacids<-unique(as.numeric(GenePathway_fattyacids$KEGG_geneId))

#EC
WTEC_diffGeneTable_fa<-merge(WTEC_diffGeneTable,GenePathway_fattyacids,all=FALSE, by='KEGG_geneId') 
pdf("WTEC_fa.pdf",onefile=FALSE)
pheatmap(WTEC_diffGeneTable_fa[,c(6:11)], scale="row", labels_row= as.character(WTEC_diffGeneTable_fa$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs EC Fatty acids ")
dev.off()

#MU
WTMU_diffGeneTable_fa<-merge(WTMU_diffGeneTable,GenePathway_fattyacids,all=FALSE, by='KEGG_geneId') 
pdf("WTMU_fa.pdf",onefile=FALSE)
pheatmap(WTMU_diffGeneTable_fa[,c(6:11)], scale="row", labels_row= as.character(WTMU_diffGeneTable_fa$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs MU Fatty acids ")
dev.off()

#ECMU
WTECMU_diffGeneTable_fa_GeneID<-unique(c(as.character(WTEC_diffGeneTable_fa$GeneId),as.character(WTMU_diffGeneTable_fa$GeneId)))
WTECMU_diffGeneTable_fa<-diffGeneTable[diffGeneTable$GeneId%in% WTECMU_diffGeneTable_fa_GeneID,]

pdf("WTECMU_fa.pdf",onefile=FALSE)
pheatmap(WTECMU_diffGeneTable_fa[,c(1:6)], scale="row", cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs EC or MU, Fatty acids")
dev.off()

##Getting phenylpropanoid  associated genes
phenylprop<-c("path:00940,path:00941,path:00944")
phenylprop<-unlist(strsplit(phenylprop,","))
GenePathway_phenylprop<-GenePathway[GenePathway$PathwayId %in% phenylprop,]
Genes_phenylprop<-unique(as.numeric(GenePathway_phenylprop$KEGG_geneId))

#EC
WTEC_diffGeneTable_phenylprop<-merge(WTEC_diffGeneTable,GenePathway_phenylprop,all=FALSE, by='KEGG_geneId') 
pdf("WTEC_phenylprop.pdf",onefile=FALSE)
pheatmap(WTEC_diffGeneTable_phenylprop[,c(6:11)], scale="row", labels_row= as.character(WTEC_diffGeneTable_phenylprop$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs EC Phenylpropanoids")
dev.off()

#MU
WTMU_diffGeneTable_phenylprop<-merge(WTMU_diffGeneTable,GenePathway_phenylprop,all=FALSE, by='KEGG_geneId') 
pdf("WTMU_phenylprop.pdf",onefile=FALSE)
pheatmap(WTMU_diffGeneTable_phenylprop[,c(6:11)], scale="row", labels_row= as.character(WTMU_diffGeneTable_phenylprop$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs MU Phenylpropanoids")
dev.off()


#ECMU
WTECMU_diffGeneTable_phenylprop_GeneID<-unique(c(as.character(WTEC_diffGeneTable_phenylprop$GeneId),as.character(WTMU_diffGeneTable_phenylprop$GeneId)))
WTECMU_diffGeneTable_phenylprop<-diffGeneTable[diffGeneTable$GeneId%in% WTECMU_diffGeneTable_phenylprop_GeneID,]

pdf("WTECMU_phenylprop.pdf",onefile=FALSE)
pheatmap(WTECMU_diffGeneTable_phenylprop[,c(1:6)], scale="row", cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs EC or MU, Phenylpropanoids")
dev.off()


##Getting sugar and nucleotide  associated genes
sunu<-c("path:00030,path:00520,path:00500,path:01210,path:00620,path:00053,path:00051,path:00052,
path:00480,path:00010,path:00630,path:00910,path:00190,path:00040,path:00230,path:00240,
path:01200")
sunu<-unlist(strsplit(sunu,","))
GenePathway_sunu<-GenePathway[GenePathway$PathwayId %in% sunu,]
Genes_sunu<-unique(as.numeric(GenePathway_sunu$KEGG_geneId))

#EC
WTEC_diffGeneTable_sunu<-merge(WTEC_diffGeneTable,GenePathway_sunu,all=FALSE, by='KEGG_geneId') 
pdf("WTEC_sunu.pdf",onefile=FALSE)
pheatmap(WTEC_diffGeneTable_sunu[,c(6:11)], scale="row", labels_row= as.character(WTEC_diffGeneTable_sunu$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs EC Sugars and Nucleotides")
dev.off()

#MU
WTMU_diffGeneTable_sunu<-merge(WTMU_diffGeneTable,GenePathway_sunu,all=FALSE, by='KEGG_geneId') 
pdf("WTMU_sunu.pdf",onefile=FALSE)
pheatmap(WTMU_diffGeneTable_sunu[,c(6:11)], scale="row", labels_row= as.character(WTMU_diffGeneTable_sunu$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs MU Sugars and Nucleotides")
dev.off()

#ECMU
WTECMU_diffGeneTable_sunu_GeneID<-unique(c(as.character(WTEC_diffGeneTable_sunu$GeneId),as.character(WTMU_diffGeneTable_sunu$GeneId)))
WTECMU_diffGeneTable_sunu<-diffGeneTable[diffGeneTable$GeneId%in% WTECMU_diffGeneTable_sunu_GeneID,]

pdf("WTECMU_sunu.pdf",onefile=FALSE)
pheatmap(WTECMU_diffGeneTable_sunu[,c(1:6)], scale="row", cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs EC or MU, Sugars and Nucleotides")
dev.off()

##Getting terpenoid associated genes
terp<-c("path:00909,path:00906,path:00904,path:00900")
terp<-unlist(strsplit(terp,","))
GenePathway_terp<-GenePathway[GenePathway$PathwayId %in% terp,]
Genes_terp<-unique(as.numeric(GenePathway_terp$KEGG_geneId))

#EC
WTEC_diffGeneTable_terp<-merge(WTEC_diffGeneTable,GenePathway_terp,all=FALSE, by='KEGG_geneId') 
pdf("WTEC_terp.pdf",onefile=FALSE)
pheatmap(WTEC_diffGeneTable_terp[,c(6:11)], scale="row", labels_row= as.character(WTEC_diffGeneTable_terp$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs EC Terpenoids")
dev.off()

#MU
WTMU_diffGeneTable_terp<-merge(WTMU_diffGeneTable,GenePathway_terp,all=FALSE, by='KEGG_geneId') 
pdf("WTMU_terp.pdf",onefile=FALSE)
pheatmap(WTMU_diffGeneTable_terp[,c(6:11)], scale="row", labels_row= as.character(WTMU_diffGeneTable_terp$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs MU Terpenoids")
dev.off()

#ECMU
WTECMU_diffGeneTable_terp_GeneID<-unique(c(as.character(WTEC_diffGeneTable_terp$GeneId),as.character(WTMU_diffGeneTable_terp$GeneId)))
WTECMU_diffGeneTable_terp<-diffGeneTable[diffGeneTable$GeneId%in% WTECMU_diffGeneTable_terp_GeneID,]

pdf("WTECMU_terp.pdf",onefile=FALSE)
pheatmap(WTECMU_diffGeneTable_terp[,c(1:6)], scale="row", cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs EC or MU, Terpenoids")
dev.off()
