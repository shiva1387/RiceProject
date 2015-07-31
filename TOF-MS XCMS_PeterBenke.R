# Centwave option, gives more peaks
library(xcms)

set <- xcmsSet(nSlaves=7,method='centWave',ppm=15,mzdiff=0.01,peakwidth=c(10,60),snthresh=6,prefilter=c(0,0))
set1 <- group(set)
set2 <- retcor(set1, method="obiwarp",profStep=0.1,plottype="deviation")
set3 <- group(set2, bw=5,mzwid=0.015,minfrac=0.25,minsamp=1)
set4 <- fillPeaks(set3)
peaklist<-peakTable(set4, filebase="amit 2wk result") 
write.table(peaklist,file="amit result.tsv",sep="\t",quote=FALSE,row.names=FALSE)

png("amit 2wk mzRT.png")
plot(peaklist$mz,peaklist$rt,pch=1,main="result",ylab="m/z",xlab="Rt[sec]")
dev.off()

library(CAMERA)
xsa<-xsAnnotate(set4)
xsaF<-groupFWHM(xsa, perfwhm=0.6)
xsaC<-groupCorr(xsaF)
xsaFI<-findIsotopes(xsaC)
xsaFA<-findAdducts(xsaFI, polarity='positive')
peaklist<-getPeaklist(xsaFA)
write.csv(peaklist, file='deisotoped_results.csv')

#reporttab <- diffreport(set4, "high", "low", "Diffreport", 1000)
#mzfilename<-"Diffreport.tsv"
#ms_data_total<-read.table(mzfilename,sep="\t",header=T,check.names=FALSE)
png("mzRT.png")
plot(peaklist$mz,peaklist$rt,pch=1,main="result",ylab="m/z",xlab="Rt[sec]")
dev.off()


# Matched Filter option, gives more peaks
library(xcms)

set <- xcmsSet(nSlaves=8,method='matchedFilter',fwhm=10,snthresh=1,step=0.1)
set1 <- group(set)
set2 <- retcor(set1,family="s",plottype="m")
set2 <- retcor(set1, method="obiwarp",profStep=0.1,plottype="deviation")
set3 <- group(set2, bw=5,mzwid=0.015,minfrac=0.25,minsamp=1)
set4 <- fillPeaks(set3)
peaklist<-peakTable(set4, filebase="blanks")
write.table(peaklist,file="blanks.tsv",sep="\t",quote=FALSE,row.names=FALSE)
png("mzRT.png")
plot(peaklist$rt,peaklist$mz,pch=1,main="Algae",ylab="m/z",xlab="Rt[min]")
dev.off()

## Deisotoping the XCMS result ##

library(CAMERA)
file <- system.file('resultsmzdata/MM14.mzdata', package = "CAMERA")
xs   <- xcmsSet(file, method="centWave",ppm=30, peakwidth=c(5,10))
an   <- xsAnnotate(xs)
an   <- groupFWHM(an)
an   <- findAdducts(an, polarity="positive")
plotEICs(an, pspec=2, maxlabel=5)
plotPsSpectrum(an, pspec=2, maxlabel=5)
