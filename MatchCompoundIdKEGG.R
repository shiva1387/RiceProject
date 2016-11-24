# Compound identification using KEGG Compound table #
# Mapping mz values extracted from XCMS data format #
#####################################################
# Author(s): Shiv, Wesley Goi (wrote the KEGG parser for extracting from KEGG Compound)
# Version: 20012016 
# Input: ".tsv" file from XCMS Metabolomics data, ".txt" file containing the compound data from KEGG
# Compound file format is
# CPD  NAME	ExactMass	MolecularMass
# Upon comparing with Agilent PCDL, we observe that the compound id is obtained by mapping 
# mz value to the ExactMass after adding (in negative mode) or removing hydrogen mass (in positive mode)
# For example, 76.038 is the mz in positive mode. It maps to Glycine (ExactMass: 75.03203 in KEGG). Thus the diff is -1.007 (mass of hydrogen)

###Criteria for selection
# 1: Peak width > 1 sec (to remove noise/artefacts)
# 2: For each pairwise comparison, the feature should be present in atleast one of the samples and in a majority of the replicates (>50%) within that species
# 3: P value (FDR) <0.05 and Fold change > + (or) - 2
# 4: Annotated (multiple, unique to species or..) to the common pool of metabolites

###################################

#Clearing the environment
rm(list=ls())
graphics.off()

#### Loading the required libraries
#library(KEGGREST)
library(dplyr)
library(VennDiagram)
library(tidyr)
library(preprocessCore)
library(data.table)

####Loading KEGG compound tables
#Kindly note : the version used is Nov2015
##Imp: This is not an organism spMUific table. This contains a list of all compounds in KEGG
setwd('/Users/shiv/Dropbox-lsius/Dropbox/Projects/NUS/RiceProject/Metabolomics/scripts/')
kegg.compound<-read.table("Compound_Aug2016.txt",header=FALSE,sep="\t",check.names=FALSE,row.names=NULL, quote="")
kegg.compound<-unique(kegg.compound) #removing duplicate rows if any
colnames(kegg.compound)<-c("CpdId","Name","ExactMass","MolecularMass")
# dim(kegg.compound)
# [1] 15676     4

# Oryza sativa japonica (Japanese rice), osa, http://www.genome.jp/kegg-bin/show_organism?org=osa

## Reading in compounds extracted from the pathway file
## Downloaded from KEGG FTP path /kegg/pathway/organisms/osa
## The file is *.list (here * is the organism id)
#OSA
#Obtained 20 January 2016
OrgSpCpd_osa_list<-read.table("../../KEGG_2016/osa_pathwayminerVersion/osa/osa.list",
                          header=FALSE,sep="\t",check.names=FALSE,row.names=NULL, quote="")
OrgSpCpd_osa<-filter(OrgSpCpd_osa_list, grepl('cpd:', OrgSpCpd_osa_list$V2)) #Extracting the column containing only compound id
OrgSpCpd_osa_cpd<-unique(as.character(OrgSpCpd_osa$V2))
length(OrgSpCpd_osa_cpd)
#3735

kegg.compound_osa<-filter(kegg.compound, CpdId %in% OrgSpCpd_osa_cpd)

#OSA
# length(unique(as.character(kegg.compound_osa$CpdId)))
# 3355

##formatting cpd id
kegg.compound_osa$CpdId<-gsub("cpd:","",kegg.compound_osa$CpdId)

# write.table(kegg.compound_osa,"TheoreticalMetabolites_OSA.txt",sep="\t",quote=FALSE,row.names=FALSE)

## Functions
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
normalizeData<-function(data_matrix){
  data_matrix<-as.matrix(log1p(data_matrix))
  processed_data<-normalize.quantiles(data_matrix,copy=TRUE)
  colnames(processed_data)<-colnames(data_matrix)
  rownames(processed_data)<-rownames(data_matrix)
  return(processed_data)
}
zero_percent_species<-function(x) {                              ### Counting number of features observed for each m/z feature in each group
  no_of_zero<-sum(x < 1e-3) 
  zero_percent<-round((no_of_zero/length(x)),2)
  return(zero_percent)
}
#loading the dataset to be mapped
#loading the full dataset
#
filename<-"../data/AR_230816_xcms1_48.tsv"
# ms_data2map<-read.csv(filename,header=TRUE)
# #For samples where the input is in the form of tsv files
# #ms_data2map<-read.table(filename,header=TRUE,stringsAsFactors = FALSE,sep="\t",row.names=1)
# ms_data2map<-ms_data2map[,c("mz","rt","isotopes","adduct","pcgroup","rtmin","rtmax")]
# #For samples which could not be processed by CAMERA
# #ms_data2map<-ms_data2map[,c("mz","rt","rtmin","rtmax")]
# rownames(ms_data2map)<-paste0(ms_data2map$mz,"@",ms_data2map$rt)
# colnames(ms_data2map)[grep('\\bmz\\b',colnames(ms_data2map))]<-"mz.orig"#replacing mz in xcms table with mz.orig as later on "mz" column
# # will contain values that have been rounded to 5 dec places
#reading differential features
ms_data2map<-read.table(filename,header=TRUE,stringsAsFactors = FALSE,sep="\t",row.names=1)
# #Creating a dataframe by splitting rownames composed of mz@rt into mz@rt,round(mz,5),round(rt/60,2)
head(ms_data2map)
colnames(ms_data2map)
ms_data2map_mzrt<-extract_mzrt(ms_data2map)

ms_data2map<-cbind(ms_data2map,ms_data2map_mzrt[,c("mz","rt.minute")])
ms_data2map$mz.rt<-rownames(ms_data2map)
colnames(ms_data2map)
#loading the dataset containing the isotope, pcgroup information from CAMERA
# Only used for differential metabolite lists
##The following section not used for AR dataset as it already contains a filtered list that has features whose
##retention time is atleast 1 second and had a abs(fold change)>=2 and FDR <= 0.05
# if(exists("ms_data_camera") ==FALSE) # Written so that if i use a diff differential list (after rea), then i dont need to read
# {
# ms_data_camera<-read.csv("../../PB_240315_xcms1_4_camera.csv",header=TRUE,stringsAsFactors = FALSE,row.names=1)  
# #For samples where the input is in the form of tsv files
# #ms_data_camera<-read.table("data/Oct2015/sar_131015_pos-exo-meoh_xcms1_42_hier_camera.csv",header=TRUE,stringsAsFactors = FALSE,sep="\t",row.names=1)
# ms_data_camera$mz.rt<-paste0(ms_data_camera$mz,"@",ms_data_camera$rt)
# ms_data_camera<-ms_data_camera[,c("rtmin","rtmax","mz.rt","isotopes","adduct","pcgroup")]#retaining only the compounds of interest
# #For samples which could not be processed by CAMERA
# #ms_data_camera<-ms_data_camera[,c("rtmin","rtmax","mz.rt")]
# }
# ms_data2map<-left_join(ms_data2map,ms_data_camera,by="mz.rt") #Only used for differential feature table
# colnames(ms_data2map)
###############################################################
###Mapping mz onto compound id
#functions to get matching mz and rt
mz_tolerance=10 #ppm
Hydrogen<-1.007

ms_data2map$Compounds = ms_data2map$mz %>% sapply(function(mz){
  ###!! Bells and Whistles!!!! For posiive mode, its mz - Hydrogen. For negative mode, its mz + Hydrogen
  ppm_tolerance<-trunc(mz - Hydrogen)*mz_tolerance/10^6
  filter(kegg.compound_osa, abs((mz - Hydrogen)-ExactMass) <=ppm_tolerance) %>% #Mapping mz onto KEGG ExactMass
    summarise(compounds=paste(CpdId,collapse=";"))
})
ms_data2map$Compounds<-as.character(unlist(ms_data2map$Compounds))
#Adding NA to features which were not annotated to any of the compounds
ms_data2map_cpds<-ms_data2map %>% mutate(Compounds = replace(Compounds, which(Compounds<0), NA))
ms_data2map_cpds$Annotation<-ms_data2map_cpds$Compounds
ms_data2map_cpds<-ms_data2map_cpds %>% mutate(Annotation = ifelse(grepl("C",Annotation),ifelse(grepl("\\;",Annotation), "Conflicting","Single"),"NA"))

#Total number of annotated features
dim(ms_data2map_cpds[complete.cases(ms_data2map_cpds[,"Compounds"]),])
#79 31 for WT vs EC
#14 31 for WT vs MU
###Writing the data out
#write.table(ms_data2map_cpds,"AR_230816_xcms1_48_annotated.txt",sep="\t",quote=FALSE,row.names=FALSE)

##The following section not used for AR dataset as it already contains a filtered list that has features whose
##retention time is atleast 1 second and had a abs(fold change)>=2 and FDR <= 0.05

### Filtering based on the above criteria
###Criteria for selection
# 1: Peak width > 1 sec
# 2: For comparison, majority of the replicate (>50%), Detected/Not detected in majority of replicates within that species
# 3: Fold change, highest & lowest > 2, Sing (3) group Vs Mixed,highest or lowest in Sing Vs Mixed
# 4: Annotated (multiple, unique to species or..) to the common pool of metabolites

#Step 1
#ms_data2map_cpds_filter<-ms_data2map_cpds %>% filter(rtmax-rtmin >=1)
#Step2 #Look at group wise presence ### Used only for pairwise comparison
#Only remove features that have majority of zeroes in both groups
# if(grep("pos-endo",filename)) { #As the sample group name differ based on extraction

#### Manually check to see if dim(zero_percent_mz_feature)==length(zero_percent_mz_feature_selection)
#grep samples using the below text, based on the files used
# #Endo_(or) AcetateExo_ (or)
# ms_data2map_cpds_filter_2<-ms_data2map_cpds_filter[, grep('MeOHExo_', colnames(ms_data2map_cpds_filter))] 
# #ms_data2map_cpds_filter_2<-ms_data2map_cpds_filter_2[1:100,]
# #Always check the groups
# Species_samples<-as.character(unlist(sapply(
#   as.character(unlist(sapply(colnames(ms_data2map_cpds_filter_2), function(x) strsplit(x,"_")[[1]][2]))),#function(x) strsplit(x,"_")[[1]][3]))), 
#   function(x) strsplit(x,"\\.")[[1]][1])))#Removing .day5 and after from the timepoint
# ms_data_table<-data.table(cbind(t(normalizeData(ms_data2map_cpds_filter_2)),as.factor(Species_samples)))
# colnames(ms_data_table)[dim(ms_data_table)[2]]<-"Species"
# zero_percent_mz_feature<-t(ms_data_table[,lapply(.SD,zero_percent_species),by=Species]) #last column had the species name
# zero_percent_mz_feature<-as.data.frame(zero_percent_mz_feature[2:nrow(zero_percent_mz_feature),])
# colnames(zero_percent_mz_feature)<-unique(Species_samples)
# dim(zero_percent_mz_feature)
# head(zero_percent_mz_feature)
# #Only remove features that have majority of zeroes in both groups
# #Check which groups are present in the dataframe and comment those which are not present
# zero_percent_mz_feature_selection<-which(
#                                            # zero_percent_mz_feature$Mixed < 0.6
#                                              zero_percent_mz_feature$KP < 0.6
#                                            # & zero_percent_mz_feature$PAO1 < 0.6
#                                            & zero_percent_mz_feature$PF < 0.6
#                                          )
# length(zero_percent_mz_feature_selection)
#Now we have a table where each column has a value expressing the % of features that are zeroes 
#(within that species, i.e replicates) for that mz
# rm(ms_data2map_cpds_filter_2,Species_samples,ms_data_table,zero_percent_mz_feature,zero_percent_mz_feature_selection)

#colnames(ms_data2map_cpds_filter)
#Step3 (Only for pairwise comparisons, as others wont have fold change)
#grep("FoldChange.",colnames(ms_data2map_cpds_filter),value = TRUE)
#The exact name of the fold change column in the dataset. Use this name in the following filter command
####$$$$####Check of the fold change is calculated on log transformed or raw value
## If log transformed, then fold change if you require fold change above 2, the condition should be abs(fold change)>1)
#ms_data2map_cpds_filter<-ms_data2map_cpds_filter %>% filter(abs(FoldChange.WT.MU.)>2)
####$$$$####
#Step4
#dim(ms_data2map_cpds_filter)
#Only used for cases where Fold Change is not available
#Selecting only features that were annotated to mz
#ms_data2map_cpds_filter_complete<-ms_data2map_cpds_filter[complete.cases(ms_data2map_cpds_filter[,"Compounds"]),]

ms_data2map_cpds_filter_complete<-ms_data2map_cpds[complete.cases(ms_data2map_cpds[,"Compounds"]),]

#Compounds detected, separating multiple hits into unique rows
compounds_detected<- ms_data2map_cpds_filter_complete %>%  mutate(Compounds = strsplit(as.character(Compounds), ";")) %>%  unnest(Compounds)
dim(compounds_detected)
#[1] 214  31 WT vs EC
#22 31 WT vs MU

# ##Added on 12-Oct-2016
# forSubPathwayAnalysis<-data.frame(compounds_detected)
# colnames(forSubPathwayAnalysis)<-gsub("\\.+2w_A\\.","",colnames(forSubPathwayAnalysis))
# rownames(forSubPathwayAnalysis)<-as.character(make.unique(paste0(forSubPathwayAnalysis$mz,"@",forSubPathwayAnalysis$rt)))
# nonRNAseq <- "R|CO|BL|Blank|rt"
# forSubPathwayAnalysis<-forSubPathwayAnalysis[, -grep(nonRNAseq, colnames(forSubPathwayAnalysis))] #removing non-rnaseq samples
# forSubPathwayAnalysis<-forSubPathwayAnalysis[, -c(1:9)] #removing non-rnaseq samples
# 
# highestCorrelations<-"E2|E3|WT1|WT2|M1|M2|Compounds"
# withinSampleCorrelations<-forSubPathwayAnalysis[,grep(highestCorrelations,colnames(forSubPathwayAnalysis))]
# 
# withinSampleCorrelations<-withinSampleCorrelations %>% group_by(Compounds) %>% summarise_each(funs(mean))
# 
# write.table(withinSampleCorrelations,"KeggCompoundMassSpec_AR_230816_xcms1_48.txt",sep="\t",quote=FALSE,col.names=NA)

#######

unique_compounds<-unique(as.character(compounds_detected$Compounds))
unique_compounds <- kegg.compound_osa[kegg.compound_osa$CpdId %in% unique_compounds,c("CpdId","Name")]
dim(unique_compounds)
#172 2 WT vs EC
#22 2 WT vs MU
write.table(unique_compounds,"differential_AR_WTvsMU_metabolites.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names = FALSE)
