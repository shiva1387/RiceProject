### Merging annotation and plotting using heatmaps with MS/MS data for getting a combined list for MS/MS #
#######################################################
## Rice metabolome
# Author: Shiv
# Input: MS/MS files generated using QTOF and converted using MD-DIAL VER 2.06
##The database used was 
#MassBank_MSMS_Pos_Rev173_vs1.msp
#Respect_20120925_ESI_Positive_MSMS.msp
# In Ms-Dial, all replicated from MU have the same class id, OE have same class id and so on...
# If the samples are blank, MU_, OE_X, US3_X the class Id's will be 1,2,3 and 4 respectively.
# All single charge adducts with 1M were selected for identification
# Version: 13052016
# Output: A merged list containg identified compounts and their corresponding m/z

# Merge these to create a unified non-redundant list of mz/rt with annotations and be passed onto Analytical chemist for verification using Mass Hunter

#This table can then be used for matching with differential features identified in MS1
#For the next step, run the mgf files created by MS-DIAL, in SIRIUS to obtain possible molecular formula, then create a similar list as above

####HeatMap plots obtained from HeatMapAnalysis_V1.R version 23032016

# Plot expression trends in metabolites
#Clearing the environment
rm(list=ls())
graphics.off()

#### Loading the required libraries
#library(KEGGREST)
library(dplyr)
library(tidyr)
library(reshape)
library(stringi)
library(stringr)

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

extract_mzrt<-function(dataset,sampleMode){
  mode_mzrt_temp <- strsplit(rownames(dataset), "\\@")
  mz<-as.numeric(sapply(mode_mzrt_temp , function (x) if(length(x) == 2) x[1] else as.character(NA)))
  mz<-round(mz,5)
  rt<-as.numeric(sapply(mode_mzrt_temp , function (x) if(length(x) == 2) x[2] else as.character(NA)))
  rt<-round(rt/60,2)
  mode_mzrt<-data.frame(rownames(dataset),mz,rt)
  colnames(mode_mzrt)<-c("mz@rt","mz","rt.minute")
  return(mode_mzrt)
}

####################################################################################################################################
## Important section of code! Mapping MS/MS annotations onto heatmaps
##################################################################################################################################
## Reading in the MS/MS data
setwd('/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/ricemanuscript/MSMS_2_Identified/')

PositiveMode_MB_msms<-read.table("20165131426_peakID_0_MassBank.txt",header=FALSE,sep="\t",check.names=FALSE, quote="", fill=TRUE)
PositiveMode_Rs_msms<-read.table("20165131442_peakID_1_ResPEct.txt",header=FALSE,sep="\t",check.names=FALSE, quote="", fill=TRUE)

###Before ANALYZING the the peakid files, remove the first 3 rows which contain Class, Type and Injection order information 
PositiveMode_MB_msms<-PositiveMode_MB_msms[-c(1:3),]
colnames(PositiveMode_MB_msms)<-as.character(unlist(PositiveMode_MB_msms[1,]))
PositiveMode_MB_msms<-PositiveMode_MB_msms[-c(1),]

PositiveMode_Rs_msms<-PositiveMode_Rs_msms[-c(1:3),]
colnames(PositiveMode_Rs_msms)<-as.character(unlist(PositiveMode_Rs_msms[1,]))
PositiveMode_Rs_msms<-PositiveMode_Rs_msms[-c(1),]

#Selecting columns to merge
PositiveMode_MB_msms_id<-PositiveMode_MB_msms[,c("Average Rt(min)","Average Mz","Metabolite name","LINK")]
PositiveMode_Rs_msms_id<-PositiveMode_Rs_msms[,c("Average Rt(min)","Average Mz","Metabolite name","LINK")]

colnames(PositiveMode_MB_msms_id)<-c("RT","MZ","MetaboliteName","LINK")
colnames(PositiveMode_Rs_msms_id)<-c("RT","MZ","MetaboliteName","LINK")

#Removing features with missing values
PositiveMode_MB_msms_id<-PositiveMode_MB_msms_id[PositiveMode_MB_msms_id$MetaboliteName!="Unknown",]
PositiveMode_Rs_msms_id<-PositiveMode_Rs_msms_id[PositiveMode_Rs_msms_id$MetaboliteName!="Unknown",]

PositiveMode_MB_msms_id$MZ<-as.numeric(as.character(PositiveMode_MB_msms_id$MZ))
PositiveMode_Rs_msms_id$MZ<-as.numeric(as.character(PositiveMode_Rs_msms_id$MZ))
PositiveMode_MB_msms_id$Kegg<-gsub("KEGG ","*",str_extract(PositiveMode_MB_msms_id$LINK, "KEGG\\s+(C\\d{5})"))
PositiveMode_Rs_msms_id$Kegg<-gsub("KEGG ","*",str_extract(PositiveMode_Rs_msms_id$LINK, "KEGG\\s+(C\\d{5})"))

PositiveMode_msms_annotation<-rbind(PositiveMode_MB_msms_id,PositiveMode_Rs_msms_id)

PositiveMode_msms_annotation$Identifier<-gsub("w/o MS2:","^",paste0(as.character(sapply(as.character(PositiveMode_msms_annotation$MetaboliteName), 
                                                                                        function(x) strsplit(x,";")[[1]][1])),";",PositiveMode_msms_annotation$Kegg))

write.table(PositiveMode_msms_annotation,"PositiveMode_msms_annotation.txt",sep="\t",quote=FALSE,col.names=NA)
# write.table(NegativeMode_msms_annotation,"NegativeMode_msms_annotation.txt",sep="\t",quote=FALSE,col.names=NA)

#### Write a function to match the mz rt from ms/ms data with the differential ms1 data in the respective modes
###############################################################
###Mapping mz onto compound id
#functions to get matching mz and rt
# mz_tolerance=5 #ppm
# 
# ###Positive mode
# ms_data2map_Pos$MSMS = ms_data2map_Pos$mz %>% sapply(function(mz){
#   ###!! Bells and Whistles!!!! For positive mode, its mz - Hydrogen. For negative mode, its mz + Hydrogen
#   ppm_tolerance<-mz*mz_tolerance/10^6
#   filter(PositiveMode_msms_annotation, abs(mz - MZ) <=ppm_tolerance) %>% #Mapping mz onto KEGG ExactMass
#     summarise(msms_ids=paste(Kegg,collapse=";"))
# })
# ms_data2map_Pos$MSMS<-as.character(sapply(as.character(unlist(ms_data2map_Pos$MSMS)),function(x) {compound_ids=unlist(strsplit(x,";"));
# paste(unique(compound_ids),collapse=";")}))
# 
# 
# ###Negative mode
# ms_data2map_Neg$MSMS = ms_data2map_Neg$mz %>% sapply(function(mz){
#   ###!! Bells and Whistles!!!! For Negative mode, its mz - Hydrogen. For negative mode, its mz + Hydrogen
#   ppm_tolerance<-mz*mz_tolerance/10^6
#   filter(NegativeMode_msms_annotation, abs(mz - MZ) <=ppm_tolerance) %>% #Mapping mz onto KEGG ExactMass
#     summarise(msms_ids=paste(Kegg,collapse=";"))
# })
# ms_data2map_Neg$MSMS<-as.character(sapply(as.character(unlist(ms_data2map_Neg$MSMS)),function(x) {compound_ids=unlist(strsplit(x,";"));
# paste(unique(compound_ids),collapse=";")}))
# 
# 
#### Write a function to match the mz from ms/ms data with the mz of ms1 data in the respective modes
###############################################################
###Mapping mz onto compound id
#functions to get matching mz and rt
# mz_tolerance=5 #ppm
# 
# ###Positive mode
# ms_data_raw_Pos_mzrt$MSMS = ms_data_raw_Pos_mzrt$mz %>% sapply(function(mz){
#   ###!! Bells and Whistles!!!! For positive mode, its mz - Hydrogen. For negative mode, its mz + Hydrogen
#   ppm_tolerance<-mz*mz_tolerance/10^6
#   filter(PositiveMode_msms_annotation, abs(mz - MZ) <=ppm_tolerance) %>% #Mapping mz onto KEGG ExactMass
#     summarise(msms_ids=paste(Identifier,collapse=";"))
# })
# ms_data_raw_Pos_mzrt$MSMS<-as.character(sapply(as.character(unlist(ms_data_raw_Pos_mzrt$MSMS)),function(x) {compound_ids=unlist(strsplit(x,";"));
# paste(unique(compound_ids),collapse=";")}))
# 
# ###Negative mode
# ms_data_raw_Neg_mzrt$MSMS = ms_data_raw_Neg_mzrt$mz %>% sapply(function(mz){
#   ###!! Bells and Whistles!!!! For Negative mode, its mz - Hydrogen. For negative mode, its mz + Hydrogen
#   ppm_tolerance<-mz*mz_tolerance/10^6
#   filter(NegativeMode_msms_annotation, abs(mz - MZ) <=ppm_tolerance) %>% #Mapping mz onto KEGG ExactMass
#     summarise(msms_ids=paste(Identifier,collapse=";"))
# })
# ms_data_raw_Neg_mzrt$MSMS<-as.character(sapply(as.character(unlist(ms_data_raw_Neg_mzrt$MSMS)),function(x) {compound_ids=unlist(strsplit(x,";"));
# paste(unique(compound_ids),collapse=";")}))

#write.table(ms_data_raw_Pos_mzrt,"PositiveMode_ms1WithMSMS_annotation.txt",sep="\t",quote=FALSE,col.names=NA)

#### Write a function to match the mz from ms/ms data with the mz of for MSMS list sent to Peter for data in the respective modes
###############################################################
forMsMS_data_raw_Pos<-read.table("../../Metabolomics/data/MetaboliteIdentification/Jan2016/forMSMS_Rice_20Jan2016.txt",header=TRUE,stringsAsFactors = FALSE,check.names=FALSE,fill =TRUE )  
colnames(forMsMS_data_raw_Pos)

###Mapping mz onto compound id
#functions to get matching mz and rt
mz_tolerance=5 #ppm

###Positive mode
forMsMS_data_raw_Pos$MSMS = forMsMS_data_raw_Pos$mz %>% sapply(function(mz){
  ###!! Bells and Whistles!!!! For positive mode, its mz - Hydrogen. For negative mode, its mz + Hydrogen
  ppm_tolerance<-mz*mz_tolerance/10^6
  filter(PositiveMode_msms_annotation, abs(mz - MZ) <=ppm_tolerance) %>% #Mapping mz onto KEGG ExactMass
    summarise(msms_ids=paste(Identifier,collapse=";"))
})
forMsMS_data_raw_Pos$MSMS<-as.character(sapply(as.character(unlist(forMsMS_data_raw_Pos$MSMS)),function(x) {compound_ids=unlist(strsplit(x,";"));
paste(unique(compound_ids),collapse=";")}))

write.table(forMsMS_data_raw_Pos,"forMSMS_Peat_pos_11Feb2016_WithMSMS_annotation.txt",sep="\t",quote=FALSE,col.names=NA)
