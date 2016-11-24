############################################################################
# Extraction of .mzxml data using XCMS in R-Metabolomics data Prof Prakash #
############################################################################
# Author(s): Shiv
# Version: 30082016
# Running Amit's 4 WEEK data from AB SCiex with latest xcms

# Software: XCMS
## SCRIPT for Rice metabolomics data from Ram (NUS-DBS)
#MS analysis performed by Peter Benke and Amit Rai -There are 2 datasets
#MS date: 2014 
#Data location 1 (statistician) (both raw and processed) : F:\Projects\NUS\RiceProject\Metabolomics\data
#Sample extraction was performed by Ram/Peter (NERI), MS analysis by Peter Benke (NERI)
## Data analysis for positive mode

#To identify a set of metabolites characteristic of gain of function and loss of function mutant in Rice
#There are 4 biological replicates and two technical replicate for each biological replicate- 
#Total number of sample = 6*4*2= 48 + 9 Blanks --> 57 columns

#The profiles are generated for  #(first X denotes biorep and second X techrep number)
#Lines                  GivenName       ModifiedName(SampleGroup)
#WT                     - WTXTX         WT
#Mutant                 - MXTX          MU
#RNAi neutral           - RNXTX         RN
#RNAi lines             - RXTX          RL
#Complementation strain - COXTX         CO
#Blank                  - BlankX        BL


###Important

# Raw data from Amit was converted by Peter to mzxml format
# The raw data is stored  in the smbl server
# 
# /data1/Projects/peter/SMBL_MS_data/Amit/Rice/2012-Nov/2-week
# /data1/Projects/peter/SMBL_MS_data/Amit/Rice/2012-Nov/4-week
# 
# The converted mzxml files for the same are stored at 
# 
# /data1/Projects/peter/SMBL_MS_data/Amit/Rice/mzxml/2 Week AR
# /data1/Projects/peter/SMBL_MS_data/Amit/Rice/mzxml/4 Week AR 
# 
# A subset of the same mzxml files for only 2-week old overexpression, mutant and wild-type lines are located at
# /data1/Projects/peter/SMBL_MS_data/Amit/Rice/mzxml/Kumar
# 
# One more copy of the same set of files are stored at
# /data1/Projects/peter/SMBL_MS_data/Amit/Rice/mzxml/26th (2week rice)
# 
# Teh 2 week run was performed on 26-April, while the 4-week run was carried out on 12-November

############


#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 
library(xcms)
library(CAMERA)

# sessionInfo()
# R version 3.3.1 (2016-06-21)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] CAMERA_1.28.0        xcms_1.48.0          Biobase_2.32.0       ProtGenerics_1.4.0   BiocGenerics_0.18.0  mzR_2.6.3           
# [7] Rcpp_0.12.5          BiocInstaller_1.22.3 packrat_0.4.3       
# 
# loaded via a namespace (and not attached):
#   [1] igraph_1.0.1        graph_1.50.0        Formula_1.2-1       magrittr_1.5        cluster_2.0.4       splines_3.3.1       munsell_0.4.3      
# [8] colorspace_1.2-6    lattice_0.20-33     plyr_1.8.4          tools_3.3.1         nnet_7.3-12         grid_3.3.1          data.table_1.9.6   
# [15] gtable_0.2.0        latticeExtra_0.6-28 survival_2.39-5     RBGL_1.48.1         Matrix_1.2-6        gridExtra_2.2.1     RColorBrewer_1.1-2 
# [22] ggplot2_2.1.0       acepack_1.3-3.3     codetools_0.2-14    rpart_4.1-10        scales_0.4.0        Hmisc_3.17-4        stats4_3.3.1       
# [29] chron_2.3-47        foreign_0.8-66  
#############
# User      #
# Specific  #
# Variables #
#############
#Raw data stored in external hard disk
setwd('C:/Users/lsius/Desktop/AmitData/4 week AR')

## Parameters for AB SCIEX

set1<-xcmsSet(nSlaves=12,method='centWave',ppm=15,peakwidth=c(10,60), prefilter=c(0,0),snthresh=6,mzdiff=0.01,)
save(set1,file="set1_AR_WK4_300816_xcms1_48.rda")
set2 <- group(set1,bw=5,mzwid=0.015,minfrac=0.66) 
set3 <- retcor(set2,method="obiwarp",plottype="none")
set4 <- group(set3,bw=5,mzwid=0.015,minfrac=0.66)
save(set4,file="set4_AR_WK4_300816_xcms1_48.rda")

set5 <- fillPeaks(set4) 
save(set5,file="set5_AR_WK4_300816_xcms1_48.rda")

peaklist<-peakTable(set5,filebase="AR_WK4_300816_xcms1_48")


xsa<-xsAnnotate(set5)
xsaF<-groupFWHM(xsa, perfwhm=0.6)
xsaC<-groupCorr(xsaF)
xsaFI<-findIsotopes(xsaC)
xsaFA<-findAdducts(xsaFI, polarity='positive')
peaklist<-getPeaklist(xsaFA)
write.csv(peaklist, file='AR_WK4_300816_xcms1_48_camera.csv')
