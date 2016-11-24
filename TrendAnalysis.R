#####################################################
# Data Analysis in R-Metabolomics data Prof Prakash #
#####################################################
# Author(s): Shiv
# Version: 11082016
# Analysis inverse expression trends in genes and metabolites

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
setwd('../../Transcriptomics/data/')

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

compute_linearModel<-function(gene_matrix,dependent.factor) { #dependent.factor is either RunDay_day4 or 12 (or) SampleGroup_day4 or 12 
  lm_model<-apply(gene_matrix,1, function(x) {
    lm_val<-lm(x~ as.factor(dependent.factor)) 
    lm_r2<-summary(lm_val)
    p.val<-anova(lm_val)$'Pr(>F)'[1]
    return(list(lm_r2$r.squared,p.val))
  })
}

#function to extract r2 value from list containing r2 and p.val returned from linear model
compute.r2.pval<-function(linearmodel_list,r2.pval) {
  if(r2.pval=="r2") { #WARNING:code implictly assumes r2 is in the first column and p.val in the second
    return (sapply(linearmodel_list, function(x){as.numeric(x[1])}))
  } else{
    return (sapply(linearmodel_list, function(x){as.numeric(x[2])}))
  }
}

###########

diffGeneTable<-read.table("differentialGeneTable.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
diffGeneTable$GeneId<-rownames(diffGeneTable)
diffGeneTable_V1<-data.matrix(cbind(diffGeneTable$GeneId,normalizeData(as.matrix(diffGeneTable[,c(1:6)]))))
colnames(diffGeneTable_V1)[1]<-"GeneId"

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

rownames(forPlot)<-gsub("EX","EX",rownames(forPlot))
forPlot$dataType<-gsub("EX","EX",forPlot$dataType)

##plotting
plot1<- ggplot(data=forPlot, aes(x=PCaxisA, y=PCaxisB, colour= factor(SampleGroups), shape= factor(SampleGroups))) + geom_point(size=3, alpha = 0.5) #dataType
plot2<- plot1 + # geom_text(data = forPlot, aes(x=PCaxisA, y=PCaxisB, label= rownames(forPlot)), hjust = 1, vjust=1, size=2) +
  theme_bw() + theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     #panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black'),legend.position="none")+ 
  xlab(paste0("PC 1 loadings","\n",round(residual_variance[1]*100,2),"% of total variation")) + 
  ylab(paste0("PC 2 loadings","\n",round(residual_variance[2]*100,2),"% of total variation")) +
  ggtitle("PCA RNA-seq")

# ggsave("pca_pc12_rnaseq.pdf",useDingbats=FALSE, h=4,w=6,units="cm",scale = 2)
# plot2
# dev.off()

# - (metabolomics data for samples with rnaseq + without outliers)
# 
# pdf("pca_pc12_rnaseq.pdf",height=4,width=4)
# plot2
# dev.off()

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
#WTEC_diff<-WTEC[WTEC$logFC >= 1 | WTEC$logFC <= -1 & WTEC$FDR <=0.05,c("X","logFC","FDR")] 
WTEC_diff<-WTEC[(WTEC$logFC >= 1 | WTEC$logFC <= -1) & WTEC$FDR <=0.05,c("X","logFC","FDR")] #Modified on 11-Aug 2016 as the previous code will select genes if it satisfies either condition!
colnames(WTEC_diff)<-c("GeneId","logFC_WTEC","FDR(p-value)_WTEC")
# nrow(WTEC_diff)
# [1] 3705
# [1] 3669 #Modified on 11-Aug 2016 as the previous code will select genes if it satisfies either condition! #There are now 36 genes less

nrow(WTEC_diff[WTEC_diff$logFC_WTEC>=0,])
#Up-regulated in EC 
#1668 #prev
#1632 #11-Aug
nrow(WTEC_diff[WTEC_diff$logFC_WTEC<=0,])
#Down-regulated in EC 
#2037
WTEC_diff<-merge(WTEC_diff,GeneDesc_GeneId,all=FALSE, by='GeneId', all.x= TRUE)
###NOTE: All logFC_WTEC regulation arewith respect to WT. As we are interested in the over-expression or mutant
### i have reversed the colors. Therefore, those 
WTEC_diff<- mutate(WTEC_diff, KEGGcolor = ifelse(logFC_WTEC <= 0 , "blue","red")) 
#It is red for up-regulation and blue for down, colors are based on expression patterns in overexpression
WTEC_diffGeneTable<-merge(WTEC_diff,diffGeneTable,all=FALSE, by='GeneId', all.x= TRUE)
##Adding gene id, note, one gene id may may to multiple locus id
WTEC_diffGeneTable_loc<-merge(WTEC_diffGeneTable,LocId,all=FALSE, by='GeneId', all.x= TRUE) 
# write.table(WTEC_diffGeneTable_loc,"WTEC_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

WTMU_diff<-WTMU[(WTMU$logFC >= 1 | WTMU$logFC <= -1) & WTMU$FDR <=0.05,c("X","logFC","FDR")]
colnames(WTMU_diff)<-c("GeneId","logFC_WTMU","FDR(p-value)_WTMU")
# nrow(WTMU_diff)
# [1] 1210 #prev
#1089 #11-Aug
nrow(WTMU_diff[WTMU_diff$logFC_WTMU>=0,])
#Up-regulated in MU 
#285
nrow(WTMU_diff[WTMU_diff$logFC_WTMU<=0,])
#Down-regulated in MU 
#804
WTMU_diff<-merge(WTMU_diff,GeneDesc_GeneId,all=FALSE, by='GeneId', all.x= TRUE)
WTMU_diff<- mutate(WTMU_diff, KEGGcolor = ifelse(logFC_WTMU <= 0 , "blue","red")) 
WTMU_diffGeneTable<-merge(WTMU_diff,diffGeneTable,all=FALSE, by='GeneId', all.x= TRUE) 

##Adding gene id, note, one gene id may may to multiple locus id
WTMU_diffGeneTable_loc<-merge(WTMU_diffGeneTable,LocId,all=FALSE, by='GeneId', all.x= TRUE) 
# write.table(WTMU_diffGeneTable_loc,"WTMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)
# 
# #pdf("WTMU_diffGenes.pdf",onefile=FALSE)
# pheatmap(temp, scale="row",cluster_rows = TRUE, cluster_cols = FALSE,#labels_row = as.character(WTMU_diffGeneTable[,1]),
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=8)
# #dev.off()

# 
# # #### Plotting venn diagram of differential genes
# venn.plot<-venn.diagram(list(as.character(WTEC_diff$GeneId),as.character(WTMU_diff$GeneId)), filename=NULL, margin=0.1, scaled=TRUE, fill=c("red", "green"),
#                         alpha=c(0.5,0.5), cex = 2, cat.dist=0.1, cat.fontface=4, category.names=c("WT vs EX", "WT vs Mutant"))
# pdf("DiffentialGenes_pval2fc2.pdf")
# grid.draw(venn.plot)
# dev.off()

#Merging the 2 comparisons to look at trends
merge_WTECMU<-merge(WTEC_diff,WTMU_diff,by='GeneId', all=FALSE)
rownames(merge_WTECMU)<-merge_WTECMU$GeneId
#write.table(merge_WTECMU,"Common_WTEC-WTMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)
#d3heatmap(merge_WTECMU[,c(2,4)], labRow=NULL, colors = "RdYlBu",scale="none", dendogram="both")

#Extracting Genes Up-reg in EC and Down-reg in MU
upWTEC_downWTMU<-merge_WTECMU[merge_WTECMU$logFC_WTEC>=0 & merge_WTECMU$logFC_WTMU<=0,]
nrow(upWTEC_downWTMU)
#22
#write.table(upWTEC_downWTMU,"upEC_downMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

#Extracting Genes Down-reg in EC and Up-reg in MU
downWTEC_upWTMU<-merge_WTECMU[merge_WTECMU$logFC_WTEC<=0 & merge_WTECMU$logFC_WTMU>=0,]
nrow(downWTEC_upWTMU)
#27
#write.table(downWTEC_upWTMU,"downEC_upMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

##Getting genes which have inverse trends
InverseTrends<-rbind(upWTEC_downWTMU,downWTEC_upWTMU)
InverseTrends_loc<-merge(InverseTrends,LocId,all=FALSE, by='GeneId', all.x= TRUE) 
InverseTrends_diffGeneTable_loc<-merge(InverseTrends_loc[,c(1,4)],diffGeneTable,all=FALSE, by='GeneId', all.x= TRUE) 
rownames(InverseTrends_diffGeneTable_loc)<-as.character(InverseTrends_diffGeneTable_loc[,1])
# 
# pdf("InverseTrends.pdf",onefile=FALSE)
# pheatmap(InverseTrends_diffGeneTable_loc[,c(3:8)], scale="row",cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=8)
# dev.off()
# 
# write.table(InverseTrends_loc,"InverseTrends_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

#Extracting Genes Up-reg in EC and Up-reg in MU
upWTEC_upWTMU<-merge_WTECMU[merge_WTECMU$logFC_WTEC>=0 & merge_WTECMU$logFC_WTMU>=0,]
nrow(upWTEC_upWTMU)
#164
upWTEC_upWTMU_V1<-merge(upWTEC_upWTMU[,c(1,4)],diffGeneTable_V1,all=FALSE, by='GeneId', all.x= TRUE) 
rownames(upWTEC_upWTMU_V1)<-as.character(upWTEC_upWTMU_V1[,1])
upWTEC_upWTMU_V1<-upWTEC_upWTMU_V1[,-c(1:2)]
upWTEC_upWTMU_V1<-data.matrix(upWTEC_upWTMU_V1[,order(colnames(upWTEC_upWTMU_V1),decreasing = TRUE)])

PlantConditions<-as.factor(unlist(sapply(colnames(upWTEC_upWTMU_V1), function(x) strsplit(x,"\\.")[[1]][1])))
PlantConditions<-relevel(PlantConditions,"WT")
PlantConditions_V1<-as.numeric(PlantConditions)-1
# GeneExpression<-as.numeric(upWTEC_upWTMU_V1[1,])
# #b<-runif(6,0.5,3)
# names(GeneExpression)<-PlantConditions; plot(a)
# #mt.maxT(upWTEC_upWTMU_V1,as.factor(PlantConditions),test = "f")
# b<-lm(upWTEC_upWTMU_V1[1,]~PlantConditions)
# summary(b)
# anova(b)
# b1<-mt.sample.teststat(upWTEC_upWTMU_V1[1,],PlantConditions_V1,test = "f",B=1000)
upWTEC_upWTMU_model<-mt.maxT(upWTEC_upWTMU_V1,classlabel = PlantConditions_V1,test = "f",B=1000)
length(which(upWTEC_upWTMU_model$adjp<0.05))
p.adj1<-p.adjust(upWTEC_upWTMU_model$rawp,method = "BH")
length(which(p.adj1<0.05))

#0
#write.table(upWTEC_upWTMU,"upEC_upMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

#Extracting Genes Down-reg in EC and Down-reg in MU
downWTEC_downWTMU<-merge_WTECMU[merge_WTECMU$logFC_WTEC<=0 & merge_WTECMU$logFC_WTMU<=0,]
nrow(downWTEC_downWTMU)
#600
downWTEC_downWTMU_V1<-merge(downWTEC_downWTMU[,c(1,4)],diffGeneTable_V1,all=FALSE, by='GeneId', all.x= TRUE) 
rownames(downWTEC_downWTMU_V1)<-as.character(downWTEC_downWTMU_V1[,1])
downWTEC_downWTMU_V1<-downWTEC_downWTMU_V1[,-c(1:2)]
downWTEC_downWTMU_V1<-downWTEC_downWTMU_V1[,order(colnames(downWTEC_downWTMU_V1),decreasing = TRUE)]
genes_downWTEC_downWTMU_V1<-as.character(rownames(downWTEC_downWTMU_V1))
downWTEC_downWTMU_V1<-apply(downWTEC_downWTMU_V1, 2, as.numeric)
rownames(downWTEC_downWTMU_V1)<-genes_downWTEC_downWTMU_V1
head(downWTEC_downWTMU_V1)
class(downWTEC_downWTMU_V1)

PlantConditions<-as.factor(unlist(sapply(colnames(downWTEC_downWTMU_V1), function(x) strsplit(x,"\\.")[[1]][1])))
PlantConditions<-relevel(PlantConditions,"WT")
PlantConditions_V1<-as.numeric(PlantConditions)-1
# GeneExpression<-as.numeric(upWTEC_upWTMU_V1[1,])
# #b<-runif(6,0.5,3)
# names(GeneExpression)<-PlantConditions; plot(a)
# #mt.maxT(upWTEC_upWTMU_V1,as.factor(PlantConditions),test = "f")
# b<-lm(upWTEC_upWTMU_V1[1,]~PlantConditions)
# summary(b)
# anova(b)
# b1<-mt.sample.teststat(upWTEC_upWTMU_V1[1,],PlantConditions_V1,test = "f",B=1000)
downWTEC_downWTMU_lm_model<-compute_linearModel(downWTEC_downWTMU_V1,PlantConditions)
downWTEC_downWTMU_lm_r.sq<-compute.r2.pval(downWTEC_downWTMU_lm_model,"r2")
downWTEC_downWTMU_lm_pval<-compute.r2.pval(downWTEC_downWTMU_lm_model,"pval")

downWTEC_downWTMU_lm_pval_v1<-p.adjust(downWTEC_downWTMU_lm_pval,method = "bonferroni")
length(which(downWTEC_downWTMU_lm_pval_v1<0.05))

downWTEC_downWTMU_model<-mt.maxT(downWTEC_downWTMU_V1,classlabel = PlantConditions_V1,test = "f",B=100)
p.adj1<-p.adjust(downWTEC_downWTMU_model$rawp,method = "BH")
length(which(p.adj1<0.05))
length(which(downWTEC_downWTMU_model$adjp<0.05))
rownames(downWTEC_downWTMU_model[which(downWTEC_downWTMU_model$adjp<0.05),])
#1
#"OS09G0344500"

# > downWTEC_downWTMU[downWTEC_downWTMU$GeneId==rownames(downWTEC_downWTMU_model[which(downWTEC_downWTMU_model$adjp<0.05),]),]
# GeneId logFC_WTEC FDR(p-value)_WTEC KEGG_geneId.x KEGGcolor.x logFC_WTMU FDR(p-value)_WTMU KEGG_geneId.y KEGGcolor.y
# OS09G0344500 OS09G0344500  -7.096603           6.4e-71       4346795        blue  -5.452687          2.73e-56       4346795        blue

# > downWTEC_downWTMU_V1[rownames(downWTEC_downWTMU_V1)==rownames(downWTEC_downWTMU_model[which(downWTEC_downWTMU_model$adjp<0.05),]),]
# WT.2      WT.1      MU.2      MU.1      EC.2      EC.1 
# 2.9897242 3.0406705 0.3740199 0.3740199 0.2189517 0.2024276 

#write.table(downWTEC_downWTMU,"downEC_downMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

##Genes differential only in WTEC
onlyWTEC_notWTMU<-anti_join(WTEC_diff, WTMU_diff, by="GeneId")
nrow(onlyWTEC_notWTMU)
#2856
#write.table(onlyWTEC_notWTMU,"onlyWTEC_notWTMU_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

##Genes differential only in WTMU
onlyWTMU_notWTEC<-anti_join(WTMU_diff, WTEC_diff, by="GeneId")
nrow(onlyWTMU_notWTEC)
#276
#write.table(onlyWTMU_notWTEC,"onlyWTMU_notWTEC_differentialGenes.txt",sep="\t",quote=FALSE,row.names=FALSE)

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
#write.table(WTEC_genePathwayInfo_longReg_count,"WTEC_genePathwayInfo_longReg_count.txt",sep="\t",quote=FALSE,row.names=FALSE)

WTEC_genePathwayInfo_longReg_count$sample<-rep("EX",nrow(WTEC_genePathwayInfo_longReg_count))

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

#write.table(WTMU_genePathwayInfo_longReg_count,"WTMU_genePathwayInfo_longReg_count.txt",sep="\t",quote=FALSE,row.names=FALSE)

WTMU_genePathwayInfo_longReg_count$sample<-rep("MU",nrow(WTMU_genePathwayInfo_longReg_count))

WTECMU_pathway_trends<-rbind(WTEC_genePathwayInfo_longReg_count,WTMU_genePathwayInfo_longReg_count)
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
# pdf("WTEC_aa.pdf",onefile=FALSE)
# pheatmap(WTEC_diffGeneTable_aa[,c(6:11)], scale="row", labels_row= as.character(WTEC_diffGeneTable_aa$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs EC Amino acids ")
# dev.off()

#MU
WTMU_diffGeneTable_aa<-merge(WTMU_diffGeneTable,GenePathway_aminoacids,all=FALSE, by='KEGG_geneId') 
# pdf("WTMU_aa.pdf",onefile=FALSE)
# pheatmap(WTMU_diffGeneTable_aa[,c(6:11)], scale="row", labels_row= as.character(WTMU_diffGeneTable_aa$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs MU Amino acids ")
# dev.off()

#ECMU
WTECMU_diffGeneTable_aa_GeneID<-unique(c(as.character(WTEC_diffGeneTable_aa$GeneId),as.character(WTMU_diffGeneTable_aa$GeneId)))
WTECMU_diffGeneTable_aa<-diffGeneTable[diffGeneTable$GeneId%in% WTECMU_diffGeneTable_aa_GeneID,]

# pdf("WTECMU_aa.pdf",onefile=FALSE)
# pheatmap(WTECMU_diffGeneTable_aa[,c(1:6)], scale="row", cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs EC or MU, Amino acids ")
# dev.off()

#### boxplot
a<-melt(WTECMU_diffGeneTable_aa)
a$Genotype<-as.character(unlist(sapply(as.character(a$variable), function(x) strsplit(x,"\\.")[[1]][1])))
a$Genotype<-gsub("EX","EX",a$Genotype)
colour_type<-gg_color_hue(length(unique(a$Genotype)))
b<-a %>% group_by(Genotype,GeneId) %>% mutate(avgGE=mean(value)) %>% as.data.frame
b<-b[,c(1,2,4,5)]
  
plot1<-ggplot(b,aes(x=Genotype,y= avgGE)) + geom_boxplot(aes(alpha=0.05, fill = Genotype),outlier.shape=NA) + #facet_grid(CompoundCategory~., drop = TRUE) + 
  scale_fill_manual(name = 'Genotype', values = colour_type) + scale_y_continuous(limits = c(0,110), breaks=seq(0,110,by=20)) +#geom_jitter(alpha=0.5,size=1) +
  xlab("Genotype") +
  ylab("log(gene expression)")  + 
  theme_bw() + theme(axis.text.x=element_text(size=8,angle=45, vjust=0.5),axis.text.y=element_text(size=8),strip.text.y = element_text(size = 6),
                     axis.line = element_line(size=1, colour = "black"), 
                     panel.grid.major = element_line(colour = "#d3d3d3"),
                     #panel.grid.minor = element_blank(), # to x remove gridlines
                     panel.background = element_blank(),
                     legend.position = "none") + ggtitle("WT Vs EX or Mu, Amino acids")

ggsave("WTECMU_aa_boxplot.pdf",plot1,w=3,h=3)

##Getting fatty acid  associated genes
fattyacids<-c("path:01212,path:00592,path:00590,path:01040,path:00061,path:00071,path:00591,path:00062")
fattyacids<-unlist(strsplit(fattyacids,","))
GenePathway_fattyacids<-GenePathway[GenePathway$PathwayId %in% fattyacids,]
Genes_fattyacids<-unique(as.numeric(GenePathway_fattyacids$KEGG_geneId))

#EC
WTEC_diffGeneTable_fa<-merge(WTEC_diffGeneTable,GenePathway_fattyacids,all=FALSE, by='KEGG_geneId') 
# pdf("WTEC_fa.pdf",onefile=FALSE)
# pheatmap(WTEC_diffGeneTable_fa[,c(6:11)], scale="row", labels_row= as.character(WTEC_diffGeneTable_fa$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs EC Fatty acids ")
# dev.off()

#MU
WTMU_diffGeneTable_fa<-merge(WTMU_diffGeneTable,GenePathway_fattyacids,all=FALSE, by='KEGG_geneId') 
# pdf("WTMU_fa.pdf",onefile=FALSE)
# pheatmap(WTMU_diffGeneTable_fa[,c(6:11)], scale="row", labels_row= as.character(WTMU_diffGeneTable_fa$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs MU Fatty acids ")
# dev.off()

#ECMU
WTECMU_diffGeneTable_fa_GeneID<-unique(c(as.character(WTEC_diffGeneTable_fa$GeneId),as.character(WTMU_diffGeneTable_fa$GeneId)))
WTECMU_diffGeneTable_fa<-diffGeneTable[diffGeneTable$GeneId%in% WTECMU_diffGeneTable_fa_GeneID,]

# pdf("WTECMU_fa.pdf",onefile=FALSE)
# pheatmap(WTECMU_diffGeneTable_fa[,c(1:6)], scale="row", cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs EC or MU, Fatty acids")
# dev.off()

#### boxplot
a<-melt(WTECMU_diffGeneTable_fa)
a$Genotype<-as.character(unlist(sapply(as.character(a$variable), function(x) strsplit(x,"\\.")[[1]][1])))
a$Genotype<-gsub("EX","EX",a$Genotype)
colour_type<-gg_color_hue(length(unique(a$Genotype)))
b<-a %>% group_by(Genotype,GeneId) %>% mutate(avgGE=mean(value)) %>% as.data.frame
b<-b[,c(1,2,4,5)]

plot1<-ggplot(b,aes(x=Genotype,y= avgGE)) + geom_boxplot(aes(alpha=0.05, fill = Genotype),outlier.shape=NA) + #facet_grid(CompoundCategory~., drop = TRUE) + 
  scale_fill_manual(name = 'Genotype', values = colour_type) + scale_y_continuous(limits = c(0,110), breaks=seq(0,110,by=20)) +
  xlab("Genotype") +
  ylab("log(gene expression)")  + 
  theme_bw() + theme(axis.text.x=element_text(size=8,angle=45, vjust=0.5),axis.text.y=element_text(size=8),strip.text.y = element_text(size = 6),
                     axis.line = element_line(size=1, colour = "black"), 
                     panel.grid.major = element_line(colour = "#d3d3d3"),
                     #panel.grid.minor = element_blank(), # to x remove gridlines
                     panel.background = element_blank(),
                     legend.position = "none") + ggtitle("WT Vs EX or Mu, Fatty acids")

ggsave("WTECMU_fa_boxplot.pdf",plot1,w=3,h=3)

##Getting phenylpropanoid  associated genes
phenylprop<-c("path:00940,path:00941,path:00944")
phenylprop<-unlist(strsplit(phenylprop,","))
GenePathway_phenylprop<-GenePathway[GenePathway$PathwayId %in% phenylprop,]
Genes_phenylprop<-unique(as.numeric(GenePathway_phenylprop$KEGG_geneId))

#EC
WTEC_diffGeneTable_phenylprop<-merge(WTEC_diffGeneTable,GenePathway_phenylprop,all=FALSE, by='KEGG_geneId') 
# pdf("WTEC_phenylprop.pdf",onefile=FALSE)
# pheatmap(WTEC_diffGeneTable_phenylprop[,c(6:11)], scale="row", labels_row= as.character(WTEC_diffGeneTable_phenylprop$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs EC Phenylpropanoids")
# dev.off()

#MU
WTMU_diffGeneTable_phenylprop<-merge(WTMU_diffGeneTable,GenePathway_phenylprop,all=FALSE, by='KEGG_geneId') 
# pdf("WTMU_phenylprop.pdf",onefile=FALSE)
# pheatmap(WTMU_diffGeneTable_phenylprop[,c(6:11)], scale="row", labels_row= as.character(WTMU_diffGeneTable_phenylprop$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs MU Phenylpropanoids")
# dev.off()


#ECMU
WTECMU_diffGeneTable_phenylprop_GeneID<-unique(c(as.character(WTEC_diffGeneTable_phenylprop$GeneId),as.character(WTMU_diffGeneTable_phenylprop$GeneId)))
WTECMU_diffGeneTable_phenylprop<-diffGeneTable[diffGeneTable$GeneId%in% WTECMU_diffGeneTable_phenylprop_GeneID,]

# pdf("WTECMU_phenylprop.pdf",onefile=FALSE)
# pheatmap(WTECMU_diffGeneTable_phenylprop[,c(1:6)], scale="row", cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs EC or MU, Phenylpropanoids")
# dev.off()

#### boxplot
a<-melt(WTECMU_diffGeneTable_phenylprop)
a$Genotype<-as.character(unlist(sapply(as.character(a$variable), function(x) strsplit(x,"\\.")[[1]][1])))
a$Genotype<-gsub("EX","EX",a$Genotype)
colour_type<-gg_color_hue(length(unique(a$Genotype)))
b<-a %>% group_by(Genotype,GeneId) %>% mutate(avgGE=mean(value)) %>% as.data.frame
b<-b[,c(1,2,4,5)]

plot1<-ggplot(b,aes(x=Genotype,y= avgGE)) + geom_boxplot(aes(alpha=0.05, fill = Genotype),outlier.shape=NA) + #facet_grid(CompoundCategory~., drop = TRUE) + 
  scale_fill_manual(name = 'Genotype', values = colour_type) + scale_y_continuous(limits = c(0,110), breaks=seq(0,110,by=20)) +
  xlab("Genotype") +
  ylab("log(gene expression)")  + 
  theme_bw() + theme(axis.text.x=element_text(size=8,angle=45, vjust=0.5),axis.text.y=element_text(size=8),strip.text.y = element_text(size = 6),
                     axis.line = element_line(size=1, colour = "black"), 
                     panel.grid.major = element_line(colour = "#d3d3d3"),
                     #panel.grid.minor = element_blank(), # to x remove gridlines
                     panel.background = element_blank(),
                     legend.position = "none") + ggtitle("WT Vs EX or Mu, Phenylpropanoids")

ggsave("WTECMU_phenylprop_boxplot.pdf",plot1,w=3,h=3)


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
# pheatmap(WTEC_diffGeneTable_sunu[,c(6:11)], scale="row", labels_row= as.character(WTEC_diffGeneTable_sunu$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs EC Sugars and Nucleotides")
# dev.off()

#MU
WTMU_diffGeneTable_sunu<-merge(WTMU_diffGeneTable,GenePathway_sunu,all=FALSE, by='KEGG_geneId') 
pdf("WTMU_sunu.pdf",onefile=FALSE)
# pheatmap(WTMU_diffGeneTable_sunu[,c(6:11)], scale="row", labels_row= as.character(WTMU_diffGeneTable_sunu$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs MU Sugars and Nucleotides")
# dev.off()

#ECMU
WTECMU_diffGeneTable_sunu_GeneID<-unique(c(as.character(WTEC_diffGeneTable_sunu$GeneId),as.character(WTMU_diffGeneTable_sunu$GeneId)))
WTECMU_diffGeneTable_sunu<-diffGeneTable[diffGeneTable$GeneId%in% WTECMU_diffGeneTable_sunu_GeneID,]

# pdf("WTECMU_sunu.pdf",onefile=FALSE)
# pheatmap(WTECMU_diffGeneTable_sunu[,c(1:6)], scale="row", cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs EC or MU, Sugars and Nucleotides")
# dev.off()


#### boxplot
a<-melt(WTECMU_diffGeneTable_sunu)
a$Genotype<-as.character(unlist(sapply(as.character(a$variable), function(x) strsplit(x,"\\.")[[1]][1])))
a$Genotype<-gsub("EX","EX",a$Genotype)
colour_type<-gg_color_hue(length(unique(a$Genotype)))
b<-a %>% group_by(Genotype,GeneId) %>% mutate(avgGE=mean(value)) %>% as.data.frame
b<-b[,c(1,2,4,5)]

plot1<-ggplot(b,aes(x=Genotype,y= avgGE)) + geom_boxplot(aes(alpha=0.05, fill = Genotype),outlier.shape=NA) + #facet_grid(CompoundCategory~., drop = TRUE) + 
  scale_fill_manual(name = 'Genotype', values = colour_type) + scale_y_continuous(limits = c(0,110), breaks=seq(0,110,by=20)) +
  xlab("Genotype") +
  ylab("log(gene expression)")  + 
  theme_bw() + theme(axis.text.x=element_text(size=8,angle=45, vjust=0.5),axis.text.y=element_text(size=8),strip.text.y = element_text(size = 6),
                     axis.line = element_line(size=1, colour = "black"), 
                     panel.grid.major = element_line(colour = "#d3d3d3"),
                     #panel.grid.minor = element_blank(), # to x remove gridlines
                     panel.background = element_blank(),
                     legend.position = "none") + ggtitle("WT Vs EX or Mu, Sugars and Nucleotides")

ggsave("WTECMU_sunu_boxplot.pdf",plot1,w=3,h=3)


##Getting terpenoid associated genes
terp<-c("path:00909,path:00906,path:00904,path:00900")
terp<-unlist(strsplit(terp,","))
GenePathway_terp<-GenePathway[GenePathway$PathwayId %in% terp,]
Genes_terp<-unique(as.numeric(GenePathway_terp$KEGG_geneId))

#EC
WTEC_diffGeneTable_terp<-merge(WTEC_diffGeneTable,GenePathway_terp,all=FALSE, by='KEGG_geneId') 
# pdf("WTEC_terp.pdf",onefile=FALSE)
# pheatmap(WTEC_diffGeneTable_terp[,c(6:11)], scale="row", labels_row= as.character(WTEC_diffGeneTable_terp$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs EC Terpenoids")
# dev.off()

#MU
WTMU_diffGeneTable_terp<-merge(WTMU_diffGeneTable,GenePathway_terp,all=FALSE, by='KEGG_geneId') 
# pdf("WTMU_terp.pdf",onefile=FALSE)
# pheatmap(WTMU_diffGeneTable_terp[,c(6:11)], scale="row", labels_row= as.character(WTMU_diffGeneTable_terp$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs MU Terpenoids")
# dev.off()

#ECMU
WTECMU_diffGeneTable_terp_GeneID<-unique(c(as.character(WTEC_diffGeneTable_terp$GeneId),as.character(WTMU_diffGeneTable_terp$GeneId)))
WTECMU_diffGeneTable_terp<-diffGeneTable[diffGeneTable$GeneId%in% WTECMU_diffGeneTable_terp_GeneID,]

# pdf("WTECMU_terp.pdf",onefile=FALSE)
# pheatmap(WTECMU_diffGeneTable_terp[,c(1:6)], scale="row", cluster_rows = TRUE, cluster_cols = FALSE,
#          clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
#          main= "WT Vs EC or MU, Terpenoids")
# dev.off()

#### boxplot
a<-melt(WTECMU_diffGeneTable_terp)
a$Genotype<-as.character(unlist(sapply(as.character(a$variable), function(x) strsplit(x,"\\.")[[1]][1])))
a$Genotype<-gsub("EX","EX",a$Genotype)
colour_type<-gg_color_hue(length(unique(a$Genotype)))
b<-a %>% group_by(Genotype,GeneId) %>% mutate(avgGE=mean(value)) %>% as.data.frame
b<-b[,c(1,2,4,5)]

plot1<-ggplot(b,aes(x=Genotype,y= avgGE)) + geom_boxplot(aes(alpha=0.05, fill = Genotype),outlier.shape=NA) + scale_y_continuous(limits = c(0,110), breaks=seq(0,110,by=20)) +
  scale_fill_manual(name = 'Genotype', values = colour_type) +
  xlab("Genotype") +
  ylab("log(gene expression)")  + 
  theme_bw() + theme(axis.text.x=element_text(size=8,angle=45, vjust=0.5),axis.text.y=element_text(size=8),strip.text.y = element_text(size = 6),
                     axis.line = element_line(size=1, colour = "black"), 
                     panel.grid.major = element_line(colour = "#d3d3d3"),
                     #panel.grid.minor = element_blank(), # to x remove gridlines
                     panel.background = element_blank(),
                     legend.position = "none") + ggtitle("WT Vs EX or Mu,Terpenoids")

ggsave("WTECMU_terp_boxplot.pdf",plot1,w=3,h=3)

##Getting genes associated with other pathways
others<-c("path:00950,path:00564,path:00196,path:00860,path:00100,path:00960,path:00908,
path:00780,path:00710,path:00073,path:01220,path:00561,path:00562")
others<-unlist(strsplit(others,","))
GenePathway_others<-GenePathway[GenePathway$PathwayId %in% others,]
Genes_others<-unique(as.numeric(GenePathway_others$KEGG_geneId))

#EC
WTEC_diffGeneTable_others<-merge(WTEC_diffGeneTable,GenePathway_others,all=FALSE, by='KEGG_geneId') 
pdf("WTEC_others.pdf",onefile=FALSE)
pheatmap(WTEC_diffGeneTable_others[,c(6:11)], scale="row", labels_row= as.character(WTEC_diffGeneTable_others$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs EC other significant pathways")
dev.off()

#MU
WTMU_diffGeneTable_others<-merge(WTMU_diffGeneTable,GenePathway_others,all=FALSE, by='KEGG_geneId') 
pdf("WTMU_others.pdf",onefile=FALSE)
pheatmap(WTMU_diffGeneTable_others[,c(6:11)], scale="row", labels_row= as.character(WTMU_diffGeneTable_others$GeneId) ,cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs MU other significant pathways")
dev.off()


#ECMU
WTECMU_diffGeneTable_others_GeneID<-unique(c(as.character(WTEC_diffGeneTable_others$GeneId),as.character(WTMU_diffGeneTable_others$GeneId)))
WTECMU_diffGeneTable_others<-diffGeneTable[diffGeneTable$GeneId%in% WTECMU_diffGeneTable_others_GeneID,]

pdf("WTECMU_others.pdf",onefile=FALSE)
pheatmap(WTECMU_diffGeneTable_others[,c(1:6)], scale="row", cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",clustering_method="average",color = greenred(15),  fontsize_row=4,
         main= "WT Vs EC or MU, other significant pathways")
dev.off()

#### boxplot
a<-melt(WTECMU_diffGeneTable_others)
a$Genotype<-as.character(unlist(sapply(as.character(a$variable), function(x) strsplit(x,"\\.")[[1]][1])))
a$Genotype<-gsub("EX","EX",a$Genotype)
colour_type<-gg_color_hue(length(unique(a$Genotype)))
b<-a %>% group_by(Genotype,GeneId) %>% mutate(avgGE=mean(value)) %>% as.data.frame
b<-b[,c(1,2,4,5)]

plot1<-ggplot(b,aes(x=Genotype,y= avgGE)) + geom_boxplot(aes(alpha=0.05, fill = Genotype),outlier.shape=NA) + #facet_grid(CompoundCategory~., drop = TRUE) + 
  scale_fill_manual(name = 'Genotype', values = colour_type) + scale_y_continuous(limits = quantile(b$avgGE, c(0.1, 0.9))) +#geom_jitter(alpha=0.5,size=1) +
  xlab("Genotype") +
  ylab("log(gene expression)")  + 
  theme_bw() + theme(axis.text.x=element_text(size=8,angle=45, vjust=0.5),axis.text.y=element_text(size=8),strip.text.y = element_text(size = 6),
                     axis.line = element_line(size=1, colour = "black"), 
                     panel.grid.major = element_line(colour = "#d3d3d3"),
                     #panel.grid.minor = element_blank(), # to x remove gridlines
                     panel.background = element_blank(),
                     legend.position = "none") + ggtitle("WT Vs EX or Mu,other significant pathways")

ggsave("WTECMU_others_boxplot.pdf",plot1)


########## RT-PCR #########
#Plotting trends for RT-PCR

#Reading in the file
rtpcr<-read.table("../../ricemanuscript/Data/rtpcrResults.txt",header=TRUE,na.strings = "",
                  sep="\t",check.names=FALSE,stringsAsFactors = FALSE)
rtpcr_geneid<-read.table("../../ricemanuscript/Data/PrimerNo_GeneId.txt",header=TRUE,na.strings = "",
                         sep="\t",check.names=FALSE,stringsAsFactors = FALSE)
  
rtpcrNew<-rtpcr[,c("Line","Bioref1","Bioref2","PrimerNo")] %>% mutate(PrimerNoNew=zoo::na.locf(PrimerNo))
rtpcrNew<-rtpcrNew[,-c(4)]#Removing the old PrimerNo column

rtpcr_long <- melt(rtpcrNew, id=c("Line","PrimerNoNew"))
rtpcr_long$variable = factor(rtpcr_long$variable,levels=c("Bioref1","Bioref2"))
#rtpcr_long$value<-as.numeric(rtpcr_long$value)
#rtpcr_long<-rtpcr_long %>% group_by(Line,PrimerNoNew) %>% summarise_each(funs(mean), mean_rt=value)

#plot_data_long$location = factor(plot_data_long$location,levels=c("Top","Bottom","Single culture"))

rtpcr_long <- summarySE(rtpcr_long, measurevar="value", groupvars=c("Line","PrimerNoNew"))


rtpcr_long<-data.frame(rtpcr_long[,-c(3)])
rtpcr_wide<-dcast(rtpcr_long,PrimerNoNew~Line)
rtpcr_wide$Mu<-rtpcr_wide$Mu/rtpcr_wide$WT #Dividing by the expression level observed in wildType
rtpcr_wide$ECE<-rtpcr_wide$ECE/rtpcr_wide$WT#Dividing by the expression level observed in wildType
rtpcr_long<-melt(data.frame(rtpcr_wide[,c(1:3)]))

rtpcr_long<-inner_join(rtpcr_long,rtpcr_geneid)
rtpcr_long$Pathway<-factor(rtpcr_long$Pathway, levels=c("Terpenoids","FattyAcids","PhenylPropanoids","Others"))
rtpcr_long$variable<-gsub("ECE","EX",rtpcr_long$variable)
rtpcr_long$variable<-factor(rtpcr_long$variable, levels=c("EX", "Mu"))
rtpcr_long$GeneId<-factor(rtpcr_long$GeneId, levels= rtpcr_long[order(rtpcr_long$value), "GeneId"])

#Adapted from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/ [Accessed Oct 09 2015]

plot1<-ggplot(rtpcr_long, aes(x=GeneId, y=value, fill=variable)) + geom_bar(position="dodge",stat = "identity") +facet_wrap(~Pathway, scales="free_x") + geom_errorbar(aes(ymin=value - sd, ymax=value + sd),  size=1, width=.2) 
  scale_fill_manual(values=c("grey15","grey45"),name = "Perturbation") + ylab("Gene expression relative to WT") + geom_hline(yintercept=1, linetype="dashed") +                        # Expand y range
  theme_bw() +   theme(axis.text.x=element_text(size=6, angle=45, vjust = 1, hjust=1), axis.text.y=element_text(size=8),
                     strip.text.x = element_text(size = 12),axis.title.y = element_text(vjust=1),
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black'))

ggsave("rtPcr.pdf",height=2, width=2,plot1, scale=4, useDingbats=FALSE)

##Plotting GSEA

#Reading in the file
#Wtmu
wtmu_gsea<-read.table("../../ricemanuscript/Data/WTMU_GSEA.txt",header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors = FALSE)
wtec_gsea<-read.table("../../ricemanuscript/Data/WTEC_GSEA.txt",header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors = FALSE)

###For merging
wtmu_gsea_filter<-wtmu_gsea %>% filter(FDR<=0.05)
wtmu_gsea_filter<-wtmu_gsea_filter[,c("GeneSetName","Class","Category","FDR","GenesinOverlap","GeneList")]
colnames(wtmu_gsea_filter)[4:6]<-c("FDR.MU","GenesinOverlap.WTMU","GeneList.WTMU")
wtec_gsea_filter<-wtec_gsea %>% filter(FDR<=0.05)
wtec_gsea_filter<-wtec_gsea_filter[,c("GeneSetName","Class","Category","FDR","GenesinOverlap","GeneList")]
colnames(wtec_gsea_filter)[4:6]<-c("FDR.EC","GenesinOverlap.WTEC","GeneList.WTEC")

list.of.data.frames<-list(wtec_gsea_filter,wtmu_gsea_filter)
merged.data.frame <- merge_recurse(list.of.data.frames)

write.table(merged.data.frame,"GSEA_analysis_sig.txt",sep="\t",quote=FALSE,row.names = FALSE)

########
wtmu_gsea_sig_GO<-filter(wtmu_gsea, FDR<=0.05 & grepl("GO_",Category))

wtmu_gsea_lit<-filter(wtmu_gsea, FDR<=0.05 & Category=="LIT")
wtmu_gsea_lit_exp<- wtmu_gsea_lit %>%  mutate(GeneList = strsplit(as.character(GeneList), " ")) %>%  unnest(GeneList)
colnames(wtmu_gsea_lit_exp)[8]<-"LOC"
wtmu_gsea_lit_exp_loc<-inner_join(WTMU_diffGeneTable_loc[,c("logFC_WTMU","LOC")],wtmu_gsea_lit_exp,all=FALSE, by='LOC', all.x= TRUE)
wtmu_gsea_lit_exp_loc$logFC_WTMU<-wtmu_gsea_lit_exp_loc$logFC_WTMU*-1 #As we want fold changes w.r.t MU
wtmu_gsea_lit_exp_loc<- wtmu_gsea_lit_exp_loc %>% group_by(GeneSetName) %>% mutate(avgFC=mean(logFC_WTMU))

###for plotting using GOplot
wtmu_gsea_lit_exp_loc_forplot1<-wtmu_gsea_sig_GO[,c("Category","Class","GeneList","FDR")]
#Formatting acoording to the data frame required by GOplot
wtmu_gsea_lit_exp_loc_forplot1$Category<-gsub("GO_","",wtmu_gsea_lit_exp_loc_forplot1$Category)
#Splitting Class
wtmu_gsea_lit_exp_loc_forplot1$ID<-as.character(unlist(sapply(
  as.character(unlist(sapply(wtmu_gsea_lit_exp_loc_forplot1$Class, function(x) strsplit(x," ")[[1]][1]))),
  function(x) strsplit(x,",")[[1]][1])))
wtmu_gsea_lit_exp_loc_forplot1$Term<-gsub(pattern = "^\\S+\\d\\s+(.+)\\,\\s+GOslim\\S+$", replacement = "\\1", x = wtmu_gsea_lit_exp_loc_forplot1$Class)

wtmu_gsea_lit_exp_loc_forplot1$GeneList<-gsub(" ",", ",wtmu_gsea_lit_exp_loc_forplot1$GeneList)
wtmu_gsea_lit_exp_loc_forplot1<-wtmu_gsea_lit_exp_loc_forplot1[,c("Category","ID","Term","GeneList","FDR")]
colnames(wtmu_gsea_lit_exp_loc_forplot1 )[c(4,5)]<-c("Genes","adj_pval")

wtmu_gsea_lit_exp_loc_forplot2<-WTMU_diffGeneTable_loc[complete.cases(WTMU_diffGeneTable_loc[,c("logFC_WTMU","LOC")]),c("LOC","logFC_WTMU")]
colnames(wtmu_gsea_lit_exp_loc_forplot2)<-c("ID","logFC")
circ <- circle_dat(wtmu_gsea_lit_exp_loc_forplot1, wtmu_gsea_lit_exp_loc_forplot2)

# Generate a simple barplot
GOBar(subset(circ, category == 'BP'))
GOBar(circ, display = 'multiple')
GOBubble(circ, labels = 3)
GOCircle(circ)
GOCluster(circ, EC$process, clust.by = 'logFC', term.width = 2)
process<-c("phosphorylation","glycosylation")
  
chord <- chord_dat(circ, wtmu_gsea_lit_exp_loc_forplot2, EC$process)
                   

##########

ggplot(data=wtmu_gsea_lit_exp_loc, aes(x=reorder(GeneSetName, avgFC), y=logFC_WTMU,fill=avgFC)) + geom_boxplot(aes(alpha=0.5)) +
         geom_jitter(position = position_jitter(width = 0.2), aes(alpha=0.5)) + scale_fill_gradient2(low="darkgreen", high="red") + 
theme_bw() + theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     #panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black')) + ggtitle("WTMU Lit")

wtmu_gsea_gfam<-filter(wtmu_gsea, FDR<=0.05 & Category=="GFam")
wtmu_gsea_gfam_exp<- wtmu_gsea_gfam %>%  mutate(GeneList = strsplit(as.character(GeneList), " ")) %>%  unnest(GeneList)
colnames(wtmu_gsea_gfam_exp)[8]<-"LOC"
wtmu_gsea_gfam_exp_loc<-inner_join(WTMU_diffGeneTable_loc[,c("logFC_WTMU","LOC")],wtmu_gsea_gfam_exp,all=FALSE, by='LOC', all.x= TRUE)
wtmu_gsea_gfam_exp_loc$logFC_WTMU<-wtmu_gsea_gfam_exp_loc$logFC_WTMU*-1 #As we want fold changes w.r.t MU
wtmu_gsea_gfam_exp_loc<- wtmu_gsea_gfam_exp_loc %>% group_by(GeneSetName) %>% mutate(avgFC=mean(logFC_WTMU))

plot1<-ggplot(data=wtmu_gsea_gfam_exp_loc, aes(x=reorder(GeneSetName, avgFC), y=logFC_WTMU,fill=avgFC)) + geom_boxplot(aes(alpha=0.5)) +
  geom_jitter(size=2, position = position_jitter(width = 0.2)) + scale_fill_gradient2(low="darkgreen", high="red") +
  theme_bw() + theme(axis.text.x=element_text(size=8, angle=45),axis.text.y=element_text(size=12),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     #panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black')) + ggtitle("WTMU GFam")
ggsave("wtmu_gfam.png",plot1)

#EC
wtec_gsea_gfam<-filter(wtec_gsea, FDR<=0.05 & Category=="GFam")
wtec_gsea_gfam_exp<- wtec_gsea_gfam %>%  mutate(GeneList = strsplit(as.character(GeneList), " ")) %>%  unnest(GeneList)
colnames(wtec_gsea_gfam_exp)[8]<-"LOC"
wtec_gsea_gfam_exp_loc<-inner_join(WTEC_diffGeneTable_loc[,c("logFC_WTEC","LOC")],wtec_gsea_gfam_exp,all=FALSE, by='LOC', all.x= TRUE)
wtec_gsea_gfam_exp_loc$logFC_WTEC<-wtec_gsea_gfam_exp_loc$logFC_WTEC*-1 #As we want fold changes w.r.t EC
wtec_gsea_gfam_exp_loc<- wtec_gsea_gfam_exp_loc %>% group_by(GeneSetName) %>% mutate(avgFC=mean(logFC_WTEC))

plot1<-ggplot(data=wtec_gsea_gfam_exp_loc, aes(x=reorder(GeneSetName, avgFC), y=logFC_WTEC,fill=avgFC)) + geom_boxplot(aes(alpha=0.5)) +
  geom_jitter(size=2, position = position_jitter(width = 0.2)) + scale_fill_gradient2(low="darkgreen", high="red") +
  theme_bw() + theme(axis.text.x=element_text(size=8, angle=45),axis.text.y=element_text(size=12),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     #panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black')) + ggtitle("WTEC GFam")
                     
ggsave("wtec_gfam.png",plot1)
