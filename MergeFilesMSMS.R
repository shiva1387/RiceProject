# Merging files for getting a combined list for MS/MS #
#######################################################
# Author: Shiv
# Input: *_selected.txt files generated from MatchCompoundIDKEGG.R
#######################################################

#Clearing the environment
rm(list=ls())
graphics.off()

#### Loading the required libraries
#library(KEGGREST)
library(dplyr)
library(tidyr)
library(reshape)

## Functions

## Reading in the data

setwd('/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/Metabolomics/data/MetaboliteIdentification/Jan2015/')
  
WT.Vs.EC<-read.table("differential_PB_WTvsEC_forId_selected.txt",header=TRUE,sep="\t",check.names=FALSE,row.names=1, quote="")
WT.Vs.MU<-read.table("differential_PB_WTvsMU_forId_selected.txt",header=TRUE,sep="\t",check.names=FALSE,row.names=1, quote="")

#Selecting columns to merge
WT.Vs.EC_s<-cbind(WT.Vs.EC[,c("mz","rt.minute","mz.rt","Compounds","Annotation")],WT.Vs.EC[,grep('FoldChange',colnames(WT.Vs.EC))])
WT.Vs.MU_s<-cbind(WT.Vs.MU[,c("mz","rt.minute","mz.rt","Compounds","Annotation")],WT.Vs.MU[,grep('FoldChange',colnames(WT.Vs.MU))])

list.of.data.frames<-list(WT.Vs.EC_s,WT.Vs.MU_s)
merged.data.frame <- merge_recurse(list.of.data.frames)
colnames(merged.data.frame)[6:7]<-c("FoldChange WT/EC","FoldChange WT/MU")
#colnames(merged.data.frame)[6:9]<-c("FoldChange PAO/Mixed","FoldChange PAO/PF","FoldChange PF/KP","FoldChange PAO/KP")

##Adding raw data to check the absolute intensity levels for each mode
ms_data_camera<-read.csv("../../PB_240315_xcms1_4_camera.csv",header=TRUE,stringsAsFactors = FALSE,row.names=1)  
#For samples where the input is in the form of tsv files
#ms_data_camera<-read.table("../../../Oct2015/sar_131015_neg-exo-ace_xcms1_42_hier_camera.csv",header=TRUE,stringsAsFactors = FALSE,sep="\t",row.names=1)
ms_data_camera$mz.rt<-paste0(ms_data_camera$mz,"@",ms_data_camera$rt)
drops <- c("mz","rt")
ms_data_camera<-ms_data_camera[,!(names(ms_data_camera) %in% drops)]#As mz is already a column in merged.data.frame
#ms_data_camera<-ms_data_camera[,c("rtmin","rtmax","mz.rt","isotopes","adduct","pcgroup")]#retaining only the compounds of interest
#For samples which could not be processed by CAMERA
#ms_data_camera<-ms_data_camera[,c("rtmin","rtmax","mz.rt")]
ms_data2map<-left_join(merged.data.frame,ms_data_camera,by="mz.rt") #Only used for differential feature table
colnames(ms_data2map)
write.table(ms_data2map,"forMSMS_Rice_20Jan2016.txt",sep="\t",quote=FALSE,row.names = FALSE)
