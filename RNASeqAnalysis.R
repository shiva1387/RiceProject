# Author: Shiv
# Version : 25-March-2015
# Protocol followed: Anders_2013_Count-based differential expression analysis

##### Scripts

setwd('/data/metagenome_data/Desktop/ProfPrakash/data')
samples<-read.table("ExperimentalDesignRice.txt",sep="\t",header=T,check.names=FALSE)

gf= "Oryza_sativa.IRGSP-1.0.26.gtf"
bowind= "OryzaS_IRGSP_1_0_26"
cmd= with (samples, paste("/data/metagenome_data/Desktop/ProfPrakash/software/tophat-2.0.14.Linux_x86_64/tophat2 -G", gf, "-p 12 -o", LibraryName, bowind, fastq1, fastq2))

#[1] /data/metagenome_data/Desktop/ProfPrakash/software/tophat-2.0.14.Linux_x86_64/tophat2 -G /data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf -p 12 -o ECE1 /data/metagenome_data/Desktop/ProfPrakash/data/OryzaS_IRGSP_1_0_26 /data/metagenome_data/Desktop/ProfPrakash/data/ECE1_CGATGT_L002_R1_001.fastq /data/metagenome_data/Desktop/ProfPrakash/data/ECE1_CGATGT_L002_R2_001.fastq

#/data/metagenome_data/Desktop/ProfPrakash/software/tophat-2.0.14.Linux_x86_64/tophat2 -G /data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf -p 12 -o ECE2 /data/metagenome_data/Desktop/ProfPrakash/data/OryzaS_IRGSP_1_0_26 /data/metagenome_data/Desktop/ProfPrakash/data/ECE2_TGACCA_L002_R1_001.fastq /data/metagenome_data/Desktop/ProfPrakash/data/ECE2_TGACCA_L002_R2_001.fastq

#/data/metagenome_data/Desktop/ProfPrakash/software/tophat-2.0.14.Linux_x86_64/tophat2 -G /data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf -p 12 -o MU1 /data/metagenome_data/Desktop/ProfPrakash/data/OryzaS_IRGSP_1_0_26 /data/metagenome_data/Desktop/ProfPrakash/data/Mu1_ACAGTG_L002_R1_001.fastq /data/metagenome_data/Desktop/ProfPrakash/data/Mu1_ACAGTG_L002_R2_001.fastq   

#/data/metagenome_data/Desktop/ProfPrakash/software/tophat-2.0.14.Linux_x86_64/tophat2 -G /data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf -p 12 -o MU2 /data/metagenome_data/Desktop/ProfPrakash/data/OryzaS_IRGSP_1_0_26 /data/metagenome_data/Desktop/ProfPrakash/data/Mu2_GCCAAT_L002_R1_001.fastq /data/metagenome_data/Desktop/ProfPrakash/data/Mu2_GCCAAT_L002_R2_001.fastq

#/data/metagenome_data/Desktop/ProfPrakash/software/tophat-2.0.14.Linux_x86_64/tophat2 -G /data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf -p 12 -o WT1 /data/metagenome_data/Desktop/ProfPrakash/data/OryzaS_IRGSP_1_0_26 /data/metagenome_data/Desktop/ProfPrakash/data/WT1_CAGATC_L002_R1_001.fastq /data/metagenome_data/Desktop/ProfPrakash/data/WT1_CAGATC_L002_R2_001.fastq   

#/data/metagenome_data/Desktop/ProfPrakash/software/tophat-2.0.14.Linux_x86_64/tophat2 -G /data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf -p 12 -o WT2 /data/metagenome_data/Desktop/ProfPrakash/data/OryzaS_IRGSP_1_0_26 /data/metagenome_data/Desktop/ProfPrakash/data/WT2_CTTGTA_L002_R1_001.fastq /data/metagenome_data/Desktop/ProfPrakash/data/WT2_CTTGTA_L002_R2_001.fastq

#gunzip Sample_WT2/WT2_CTTGTA_L002_R1_001.fastq.gz
#gunzip Sample_WT2/WT2_CTTGTA_L002_R2_001.fastq.gz

#export PATH=$PATH:/data/metagenome_data/Desktop/ProfPrakash/software/bowtie2-2.2.5


for(i in seq_len(nrow(samples))) {

 lib = samples$LibraryName[i]

 ob = file.path(lib, "accepted_hits.bam")

 # sort by name, convert to SAM for htseq-count

 cat(paste0("samtools sort -n ",ob," ",lib,"_sn"),"\n")

 cat(paste0("samtools view -o ",lib,"_sn.sam ",lib,"_sn.bam"),"\n")

 # sort by position and index for IGV

 cat(paste0("samtools sort ",ob," ",lib,"_s"),"\n")

 cat(paste0("samtools index ",lib,"_s.bam"),"\n\n")

}

../software/samtools-1.2/samtools sort -n ECE1/accepted_hits.bam ECE1_sn 
../software/samtools-1.2/samtools view -o ECE1_sn.sam ECE1_sn.bam 
../software/samtools-1.2/samtools sort ECE1/accepted_hits.bam ECE1_s 
../software/samtools-1.2/samtools index ECE1_s.bam 

../software/samtools-1.2/samtools sort -n ECE2/accepted_hits.bam ECE2_sn 
../software/samtools-1.2/samtools view -o ECE2_sn.sam ECE2_sn.bam 
../software/samtools-1.2/samtools sort ECE2/accepted_hits.bam ECE2_s 
../software/samtools-1.2/samtools index ECE2_s.bam 

../software/samtools-1.2/samtools sort -n MU1/accepted_hits.bam MU1_sn 
../software/samtools-1.2/samtools view -o MU1_sn.sam MU1_sn.bam 
../software/samtools-1.2/samtools sort MU1/accepted_hits.bam MU1_s 
../software/samtools-1.2/samtools index MU1_s.bam 

../software/samtools-1.2/samtools sort -n MU2/accepted_hits.bam MU2_sn 
../software/samtools-1.2/samtools view -o MU2_sn.sam MU2_sn.bam 
../software/samtools-1.2/samtools sort MU2/accepted_hits.bam MU2_s 
../software/samtools-1.2/samtools index MU2_s.bam 

../software/samtools-1.2/samtools sort -n WT1/accepted_hits.bam WT1_sn 
../software/samtools-1.2/samtools view -o WT1_sn.sam WT1_sn.bam 
../software/samtools-1.2/samtools sort WT1/accepted_hits.bam WT1_s 
../software/samtools-1.2/samtools index WT1_s.bam 

../software/samtools-1.2/samtools sort -n WT2/accepted_hits.bam WT2_sn 
../software/samtools-1.2/samtools view -o WT2_sn.sam WT2_sn.bam 
../software/samtools-1.2/samtools sort WT2/accepted_hits.bam WT2_s 
../software/samtools-1.2/samtools index WT2_s.bam 


samples$countf = paste(samples$LibraryName, "count", sep=".")
gf = "/data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf"
cmd = paste0("htseq-count -s no -a 10 ", samples$LibraryName, "_sn.sam ", gf," > ", samples$countf)
cmd

python htseq-count -s no -a 10 /data/metagenome_data/Desktop/ProfPrakash/data/ECE1_sn.sam /data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf > ECE1.count
python htseq-count  -s no -a 10 /data/metagenome_data/Desktop/ProfPrakash/data/ECE2_sn.sam /data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf > ECE2.count
python htseq-count -s no -a 10 /data/metagenome_data/Desktop/ProfPrakash/data/MU1_sn.sam /data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf > MU1.count  
python htseq-count  -s no -a 10 /data/metagenome_data/Desktop/ProfPrakash/data/MU2_sn.sam /data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf > MU2.count  
python htseq-count  -s no -a 10 /data/metagenome_data/Desktop/ProfPrakash/data/WT1_sn.sam /data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf > WT1.count  
python htseq-count  -s no -a 10 /data/metagenome_data/Desktop/ProfPrakash/data/WT2_sn.sam /data/metagenome_data/Desktop/ProfPrakash/data/Oryza_sativa.IRGSP-1.0.26.gtf > WT2.count

library("edgeR")
counts = readDGE(samples$countf)$counts
counts.full = readDGE(samples$countf)$counts
noint = rownames(counts) %in%  c("no_feature","ambiguous","too_low_aQual", "not_aligned","alignment_not_unique") 
cpms = cpm(counts)
keep = rowSums(cpms >1) >=2 & !noint #2 is the number of replicates
counts = counts[keep,]

#Gene of interest

# > counts.full[63097,]
# ECE1 ECE2  MU1  MU2  WT1  WT2 
# 21   32    0    1    9    7 
# > rownames(counts.full)[63097]
# [1] "OS03G0140400"

#OS03G0140400 is almost 3-fold up in over-expression line

colnames(counts) = samples$shortname

head( counts[,order(samples$condition)], 5 )

d = DGEList(counts=counts, group=samples$condition)

d = calcNormFactors(d)

d = estimateCommonDisp(d)

d = estimateTagwiseDisp(d)

##Effective library sizes
d$samples$lib.size * d$samples$norm.factors

#EC.1 34690449
#EC.2 35551848 
#MU.1 31805211
#MU.2 36929702
#WT.1 34932760
#WT.2 35036621

pdf("meanVar.pdf")
plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
dev.off()
pdf("dispersion.pdf")
plotBCV(d)
dev.off()


de_WTEC = exactTest(d, pair=c("WT","EC"))
de_WTMU = exactTest(d, pair=c("WT","MU"))
de_ECMU = exactTest(d, pair=c("EC","MU"))

summary( decideTestsDGE( de_WTEC , p.value = 0.05 ) )

tt_WTEC = topTags(de_WTEC, n=nrow(d))
head(tt_WTEC$table)
tt_WTMU = topTags(de_WTMU, n=nrow(d))
head(tt_WTMU$table)
tt_ECMU = topTags(de_ECMU, n=nrow(d))
head(tt_ECMU$table)

nc = cpm(d, normalized.lib.sizes=TRUE)

rn_WTEC = rownames(tt_WTEC$table)
rn_WTMU = rownames(tt_WTMU$table)
rn_ECMU = rownames(tt_ECMU$table)

#head(nc[rn_WTEC,order(samples$condition)],5)

#depth-adjusted reads per million
write.table(nc[rn_WTEC,order(samples$condition)],"WTEC_differentialGeneTable.txt",sep="\t",quote=FALSE)
write.table(nc[rn_WTMU,order(samples$condition)],"WTMU_differentialGeneTable.txt",sep="\t",quote=FALSE)

pdf("SMEAR_WTEC.pdf")
deg_WTEC = rn_WTEC[tt_WTEC$table$FDR < .05]
plotSmear(d, de.tags=deg_WTEC)
dev.off()

pdf("SMEAR_WTMU.pdf")
deg_WTMU = rn_WTEC[tt_WTMU$table$FDR < .05]
plotSmear(d, de.tags=deg_WTMU)
dev.off()

pdf("SMEAR_ECMU.pdf")
deg_ECMU = rn_WTEC[tt_ECMU$table$FDR < .05]
plotSmear(d, de.tags=deg_ECMU)
dev.off()


write.csv(tt_WTEC$table, file="toptags_tt_WTEC_edgeR.csv")
write.csv(tt_WTMU$table, file="toptags_tt_WTMU_edgeR.csv")
write.csv(tt_ECMU$table, file="toptags_tt_ECMU_edgeR.csv")

