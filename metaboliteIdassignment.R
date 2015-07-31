#####################################################
# Data Analysis in R-Metabolomics data Prof Prakash #
#####################################################
# Author(s): Shiv
# Version: 17062015
# Assinging KEGG compound id's to mass features
# Using PCDL in positive mode, the m/z feature is subtracted by 
# Mass feature (m/z) in the dataset    227.10204	
# Mass submitted in PCDL(after correcting for positive mode) 226.0947635	
# Difference in masses 1.0072765
# Input file a combined list of metabolites identified for differential features when WT is compared to EC and MU
# In the input file 'MetaboliteMass' is the mass of the metabolite, 'MassSubmittedPCDL' is the m/z value submitted by PCDL
# after correcting for positive mode (i.e subtracting the original m/z obtained by xcms, by 1.0072)
######################################################

#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 

#############
# User      #
# Specific  #
# Variables #
#############

##############################################
# Reading in the combined metabolite id data #
##############################################

setwd('MetaboliteIdentification/Positive/KEGGRice_database/')

### RUN THE USER CREATED FUNCTIONS FIRST

metaboliteListFileName<-"Combined_EC-MUVsWT_DifferentialMetabolite.txt"
metaboliteList<-read.table(metaboliteListFileName,sep="\t",header=T,check.names=FALSE)

##Obtain list of differential mass features from Analysis_scripts.R
#WTvsEC_mass.list and WTvsMU_mass.list
apply(WTvsEC_mass.list,2,class)
WTvsEC_mass.list_V1 <- apply(WTvsEC_mass.list,2,as.numeric)
WTvsMU_mass.list_V1 <- apply(WTvsMU_mass.list,2,as.numeric)
##Adding version number to avoid making accidental changes to the original matrix

closestPCDLMass<-function(massList) {
    list_MassCpd<-sapply(massList, function(x)
                          {ind<-which(as.numeric(metaboliteList$MassSubmittedPCDL) > (x-1.0073) & 
                                          as.numeric(metaboliteList$MassSubmittedPCDL) < (x-1.00715));
                            unique_mass<-paste0(unique(metaboliteList$MassSubmittedPCDL[ind]),collapse=';')
                            unique_cpds<-paste0(unique(metaboliteList$KEGG[ind]),collapse=';')
                            cpd_mass<-c(x,"@",unique_mass,"@",unique_cpds)
                           return(cpd_mass)
                          }
                        )
      return(list_MassCpd)}

#WTvsEC_mass.list_V1<-cbind(WTvsEC_mass.list_V1,closestPCDLMass)

a<-closestPCDLMass(as.numeric(WTvsEC_mass.list_V1[,c(1)]))
write.table(t(a),"mapped_metabolite.txt",quote=FALSE,sep="\t",row.names=FALSE)

  which(as.numeric(metaboliteList$MassSubmittedPCDL) < (130.0494-1.0071) & 
        unique(metaboliteList$MassSubmittedPCDL[1:3]) > (130.0494-1.00735))
  


#################mapping metabolites Modified by shiv to include the actual mass also 14-Nov-2014

mappingMetabolites_mod<-function(massSpecDataMatrix,metaboliteData){  #massSpecDataMatrix is the mz data, metaboliteDat contains cpd id 
  #and mass
  mappedData<-list()
  sigfeat_bc<-splitMzRt(massSpecDataMatrix) #Change to d4 or d12 according to analysis
  for (i in 1:nrow(metaboliteData))
  {
    # STEP1:Obtain values within the range
    interMappedData<- sigfeat_bc[sigfeat_bc$sigfeat_bc_mz_s >= metaboliteData$lowerRange[i] & 
                                   sigfeat_bc$sigfeat_bc_mz_s <= metaboliteData$upperRange[i],]
    #print(paste(i,length(interMappedData)))
    if(!is.na(interMappedData[1,1]))
    {mappedData[[i]]<-cbind(metaboliteData$actualMass[i],metaboliteData$Mass[i],metaboliteData$CompoundId[i],interMappedData)}
    # STEP2: For values which do not map within the tolerance limits,
    #        obtain the mz which is closest to the value
    if(is.na(interMappedData[1,1])) {
      match_mz<-cmpMZmin(sigfeat_bc$sigfeat_bc_mz_s,metaboliteData$actualMass[i])
      interMappedData<-sigfeat_bc[match_mz,]
    }
    if(!is.na(interMappedData[1,1]))
    {mappedData[[i]]<-cbind(metaboliteData$actualMass[i],metaboliteData$Mass[i],metaboliteData$CompoundId[i],interMappedData)}
    # STEP3: Check if any mz is left out
    if(is.na(interMappedData[1,1]))
    { print(paste(metaboliteData$lowerRange[i],metaboliteData$upperRange[i]))}
  }
  mappedData1<-do.call(rbind,mappedData) #Complete list of mapped metabolites and corresponding features
  # Will contain mmultiple hits to each metabolite id
  return (mappedData1)
}


