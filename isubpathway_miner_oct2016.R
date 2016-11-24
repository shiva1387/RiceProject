################################################
# Metabolomics Pathway enrichment in Rice data #
################################################
# Author(s): Shiv
# Version: 23082016 
# Input: ".tsv" file from XCMS 
# Software: XCMS
# Modified By :Shivshankar Umashankar 
##############################################################
install.packages("../../Metabolomics/scripts//iSubpathwayMiner_3.0.tar.gz", repos = NULL, type = "source")

##Impostant: Reanalyze the sub-pathways based on improved annotation performed in Jan2015 using MatchCompoundIdKEGG.R
# source("http://bioconductor.org/biocLite.R")
# biocLite(c('graph', 'RBGL', 'igraph', 'XML'))

rm(list = ls())

library("iSubpathwayMiner")

####### Updating pathways
setwd("/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/KEGG_2016/osa_pathwayminerVersion")
path<-paste("/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/KEGG_2016/osa_pathwayminerVersion")

updateOrgAndIdType("osa","ncbi-geneid","/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/KEGG_2016/osa_pathwayminerVersion")
getOrgAndIdType()

#####getting list of significant genes and metabolites

setwd("/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/Metabolomics/scripts/")

sig_genes<-read.table('/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/Transcriptomics/data/WTMU_differentialGenes.txt',header=TRUE,sep='\t')
sig_metab<-read.table('/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/Metabolomics/data/differential_AR_WTvsMU_metabolites.txt',header=FALSE,sep='\t',fill = TRUE,check.names = FALSE, stringsAsFactors = FALSE)
gene_data<-sig_genes$KEGG_geneId
gene_data<-as.integer(na.omit(gene_data))

metab_data<-sig_metab[,1]
metab_data<-as.character(metab_data[metab_data!=""])

moleculeList1<-as.character(c(gene_data,metab_data))

##########identify metabolic subpathways by using the function SubpathwayGM###########
#get a set of interesting genes and metabolites
#moleculeList<-getExample(geneNumber=100,compoundNumber=50)
#annotate gene and metabolite sets to metabolic subpathways 
#and identify subpathways
reGM<-SubpathwayGM(moleculeList1,n=2,s=5)
#convert ann to data.frame
result<-printGraph(reGM$ann)
#print the results to screen
result[1:10,]

#Each row of the result (data.frame) is a subpathway. Its columns include pathwayId, 
#pathwayName, annMoleculeRatio, annBgRatio, pvalue, and fdr. They correspond to subpathway
#identifier, pathway name, the ratio of the annotated interesting molecules, the ratio 
#of the annotated background, p-value of the hypergeometric test, and Benjamini-Hochberg
#fdr values. For annMoleculeRatio, 29/1100 means that 29 molecules in 1100 interesting molecules
#are annotated to the subpathway. For annBgRatio, 67/25051 means that 67 molecules 
#in 25051 background molecules are annotated to the subpathway. 
result1<-printGraph(reGM$ann,detail=TRUE)
##write the annotation results to tab delimited file. 
write.table(result1,file="../../Metabolomics/data/EnrichmentAnalysis/Aug2016/WTvsMU_subpathway-result-pv05_fc2_meta_n2s5.txt",row.names=FALSE,sep="\t",quote=FALSE)

#visualize subpathways in R and KEGG web site
str(result1)
result2<-result[which(result[,"pvalue"]<0.05),]
for(i in 1:nrow(result1))
{
  pathway_id<-result1[i,1]
  plotAnnGraph(pathway_id,reGM$subGraphList,reGM$ann,displayInR=TRUE,gotoKEGG=TRUE)
}

#plotAnnGraph("path:00941",reGM$subGraphList,reGM$ann,displayInR=TRUE,gotoKEGG=TRUE,multipleCell=TRUE,match=TRUE)

#Differential genes analysis

##########get the list of pathway graphs#########
##Convert all metabolic pathways to graphs.
#convert pathways to a list in R
#setwd('../../KEGG_2016/kgml/metabolic/osa/')
setwd("/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/KEGG_2016/kgml/metabolic/osa/")
a<-list.files()
p<-getPathway("/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/KEGG_2016/kgml/metabolic/osa/",a)
gm<-getMetabolicGraph(p)
gm[[1]]$title
plotGraph(gm[[1]])
summary(gm[[1]])

graphList<-c(gm)

##get a set of genes annotate gene sets to graphs and identify significant graphs
ann<-identifyGraph(as.character(gene_data),graphList,type="gene")
#convert ann to data.frame
result<-printGraph(ann,detail=TRUE)
##write the annotation results to tab delimited file.
write.table(result,file="../../../../Metabolomics/data/EnrichmentAnalysis/Aug2016/WTvsMU_enrichment_result_pv05fc2.txt",row.names=FALSE,sep="\t",quote=FALSE)

##########annotate metabolite sets and identify enriched metabolic pathways###########

ann<-identifyGraph(metab_data,graphList,type="compound")
#convert ann to data.frame
result<-printGraph(ann,detail=TRUE)
result[which(result[,"pvalue"]<0.05),]
write.table(result,file="../../../../Metabolomics/data/EnrichmentAnalysis/Aug2016/WTvsMU_metabolite_enrichment_result_pv05.txt",row.names=FALSE,sep="\t",quote=FALSE)

##########annotate gene and compound sets and identify enriched metabolic pathways###########
#annotate gene and compound sets to metabolic graphs and identify significant graphs
ann<-identifyGraph(moleculeList1,graphList,type="gene_compound")
#convert ann to data.frame
result<-printGraph(ann)
result[,1:5]
write.table(result,file="../../../../Metabolomics/data/EnrichmentAnalysis/Aug2016/WTvsMU_gene_pv05fc2metabolite_enrichment_result_pv05.txt",row.names=FALSE,sep="\t",quote=FALSE)
