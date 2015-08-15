########### PREAMBLE ###########################################################################
## SICB Annual Meeting 2013
## January 3-7, 2013 
## San Francisco, CA
## Hilton San Francisco Union Square
##
## Workshop
## Genomics for Non-Model organisms: Custom Microarray Development and Analysis
##
## Bioinformatician: Dr. Matt Settles
##  University of Idaho - Genomics Resources Core <msettles@uidaho.edu>
##  http://webpages.uidaho.edu/msettles/
##
## Special Thanks to Tom Poorten, 
##   PhD student UC Berkeley, for his help in putting on this workshop
################################################################################################

## Required preliminaries in order to run
## I recommend Rstudio with R
## Download and install R from www.r-project.org
## Download and install Rstudio from http://www.rstudio.com/
## Need to install packages and dependancies for R
## To install packages, in R
##   source("http://bioconductor.org/biocLite.R")
##   biocLite(c("qvalue","oligo","RColorBrewer","limma","pdInfoBiulder)
##

## Set your primary directory, need to create all needed folders prior to running
## This Directory should have the following folders in it
##  Raw_Data - contains the xys files (or Affy:CEL, Agilent:gpr files)
##  Design_Files - contains the ndf and ngd files (or Agilent:gal files)
##  Data - is originally empty and is where intermidaite data is saved during analysis
##  Figures - is originally empty and where figures are saved during analysis
##  Tables - is originally empty and where result tables are saved during analysis
base <- "../SICB-example" ## assumes you start R in the SICB-example folder
setwd(base)

experimentName <- "zebrafish_nimblegen"

xysPath <- file.path(getwd(),"Raw_Data")
designPath <- file.path( getwd(),"Design_Files")
rawPath <- file.path(getwd(),"Data")
figurePath <- file.path(getwd(), "Figures")
tablePath <- file.path(getwd(), "Tables")

library(qvalue)
library(oligo)
library(RColorBrewer)
library(limma)
library(pdInfoBuilder)
library(Biostrings)
library(plotrix)
source("Venn.R") ## this file should be in the base directory, can ignore warnings
source("functions.R") ## this file should be in the base directory

#################################################################################
#################################################################################
## Nimblegen arrays require a custom package in order to preprocess
##  Creating and generating this custom package only needs to be done 
##  once for each array. ONLY WORKS ON MAC OR LINUX
##
## Build pdInfoBuilder packages, only need to do once
library(pdInfoBuilder)
ndf <- list.files(designPath, pattern = ".ndf$",full.names = TRUE)
xys <- list.files(xysPath, pattern = ".xys",full.names = TRUE)[1]
seed <- new("NgsExpressionPDInfoPkgSeed",
     ndfFile = ndf, xysFile = xys,
     author = "Matt Settles",
     email = "msettles@uidaho.edu",
     biocViews = "AnnotationData",
     genomebuild = "UCSC_Zv7",
     organism = "Zebrafish", species = "Danio Rerio",
     url = "http://webpages.uidaho.edu/msettles")
makePdInfoPackage(seed, destDir = ".",unlink=TRUE)
################################
## Then build and install into R
#####################################
system("R CMD build pd.090818.zv7.expr")
system("R CMD INSTALL pd.090818.zv7.expr_0.0.1.tar.gz")
#################################################################################
#################################################################################

## Generate an Annotation table using the ndf and ngd files
ndf <- list.files(designPath, pattern = ".ndf$",full.names = TRUE)
ndf_data <- read.table(ndf,sep="\t",header=T,as.is=T)
ngd <- list.files(designPath, pattern = ".ngd$",full.names = TRUE)
ngd_data <- read.table(ngd,sep="\t",header=T,comment.char="",quote="",as.is=T)
nms <- scan(ngd,what=character(0),n=2,sep="\t")
colnames(ngd_data) <- c("SEQ_ID","Desc")
annotation <- t(as.data.frame(strsplit(ngd_data$Desc,split="|",fixed=T)))
rownames(annotation) <- ngd_data$SEQ_ID
colnames(annotation) <- strsplit(nms[2],split="|",fixed=T)[[1]]
write.table(annotation,file=file.path(designPath,"annotation.txt"),sep="\t",row.names=T,col.names=T,quote=F)

## Start of Analysis
#################################################################################
### Load data, check quality, 
#################################################################################
# targets file is the experiment description file
targetsFile <- "targets.txt"

# phenodata, create phenotypic data object
phenodata <- read.AnnotatedDataFrame(targetsFile,widget=FALSE,header=TRUE, sep="\t",row.names="Name")
varMetadata(phenodata)$channel = factor(unique(as.character(phenodata$Dye)),levels=c("_ALL_"))   

# load data
RawData <- read.xysfiles(file.path(xysPath,phenodata$FileName),
    phenoData=phenodata,sampleNames=sampleNames(phenodata),
    notes=experimentName, verbose=TRUE, checkType=TRUE)   

## setup color scheme for treatment factors
fcol <- brewer.pal(4,"Set1")[as.numeric(as.factor(paste(RawData$Treatment)))]

### Quality Assurance of Raw Data
#Generate some pseudo images of the chips:
for(i in 1:nrow(pData(RawData))){
  png(file.path(figurePath, paste("QA_PseudoImage_", rownames(pData(RawData))[i], ".png", sep="")), width=6, height=4, units="in", pointsize=8, bg="white", res=300)
  image(RawData[,i])  #image won't take a main argument for some reason.. probably using the one in "graphics" instead of "oligo" ?
  #main=paste(phenodata$Sample_Label[i], phenodata$Description[i], sep=" : ")
  dev.off()
}

## Histogram
png(file.path(figurePath,"QA_rawhist.png"),width = 6, height = 5, 
    units = "in", pointsize = 8, bg = "white",res=300)   
hist(RawData,col=fcol, lwd=2, lty=1, main="Raw probe-level expression")
legend("topright", fill=fcol, legend=rownames(pData(RawData)))
dev.off()

## Boxplots
png(file.path(figurePath,"QA_rawboxplot.png"),width = 5, height = 5, 
    units = "in", pointsize = 8, bg = "white",res=300)   
par(mar=c(7,4,4,2) + .1)
boxplot(RawData, col=fcol, names=rownames(pData(RawData)), main="Raw probe-level expression", xaxt="n")
axis(1, at=1:12, labels=F)
text(1:12, par("usr")[3] - 0.25, srt=45, adj = 1, labels=rownames(pData(RawData)), xpd=TRUE)  
dev.off()

## Hierarchical clustering plot
data.ALL <- exprs(RawData)
clust.cor=hclust(dist(t(data.ALL)),method="complete")

png(file=file.path(figurePath,"QA_rawcluster.png"),point=8,width=10,height=7,units="in",res=300)
par(mar = c(2, 4,6,2))
plclust(clust.cor,main="Cluster Using Raw probe-level expression",sub=paste(nrow(data.ALL), "genes"),xlab="",ylab="Euclidean Distance")
dev.off()

## MDS plot
d <- dist(t(data.ALL),method="euc") # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# plot solution
png(file=file.path(figurePath,"QA_rawMDS.png"),point=8,width=10,height=7,units="in",res=300)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS Using Raw probe-level expression", type="n", xlim=range(x)*1.2, ylim=range(y)*1.2)
text(x, y, labels = sampleNames(RawData),col= fcol, cex=1, font=2)
dev.off()

## MA plots
png(file=file.path(figurePath,"QA_rawMAPlots.png"), point=8, width=24, height=24, units="in", res=300)
MAplot(RawData, pairs=T)
dev.off()

####################################################################
#Explore background: Nimblegen Only
bgSequences = bgSequence(RawData)
pmSequences = pmSequence(RawData)
bgIntensities = bg(RawData)
pmIntensities = pm(RawData)

counts.bg = alphabetFrequency(bgSequences)[,c('A','C','G','T')]
gcFreq.bg = rowSums(counts.bg[,c('C','G')])/rowSums(counts.bg)
avg.bg = apply(bgIntensities, 1, mean)

counts.pm = alphabetFrequency(pmSequences)[,c('A','C','G','T')]
gcFreq.pm = rowSums(counts.pm[,c('C','G')])/rowSums(counts.pm)
avg.pm = apply(pmIntensities, 1, mean)

png(file.path(figurePath,"QA_background_Intensities_by_GC.png"),width = 12, height = 5,units='in', bg='white', res=300)
par(mfrow=c(1,3))
boxplot(log2(avg.bg)~gcFreq.bg,main="Background intensities by GC")
boxplot(log2(avg.pm)~gcFreq.pm,main="Perfect Match intensities by GC")
tmp.bg = aggregate(log2(avg.bg)~gcFreq.bg,data.frame(avg.bg, gcFreq.bg), mean)
tmp.pm = aggregate(log2(avg.pm)~gcFreq.pm,data.frame(avg.pm, gcFreq.pm), mean)
plot(tmp.pm, type="o", ylab="mean log2 intensity",col=c("red"),xlab="GC content", main="Background/PM intensity by GC content")
lines(tmp.bg, type="o", ylab="mean log2 intensity",col=c("black"),xlab="GC content", main="Background/PM intensity by GC content")
legend("topleft",fill=c("black","red"),legend=c("Background","Perfect Match"),cex=0.75)
dev.off()

png(file.path(figurePath,"QA_probelevel_Intensities.png"),width = 5, height = 5,units='in', bg='white', res=300)
par(mfrow=c(1,1))
plot(density(log2(pm(RawData))), main="log2 probe level intensities", ylim=c(0,1.9), lwd=2)
points(density(log2(bg(RawData))), col="red", type="l", lwd=2)
legend("topleft", fill=c("black","red"), legend=c("PM", "BG"))
dev.off()

save(RawData,file=file.path(rawPath,"rawData.RData"))
#################################################################################
### Normalize and check quality
#################################################################################
## Years of experience and some research, suggests RMA is on average the best and 
##  safest preprocessing routine to use on the typical experiment
eset <- rma(RawData) ## If your going to remove samples for QA reason do it here, or in the targets file

## Histogram
png(file.path(figurePath,"RMA_hist.png"),width = 6, height = 5, 
    units = "in", pointsize = 8, bg = "white",res=300)   
hist(eset,col=fcol, lwd=2, lty=1, main="RMA normalized probeset expression values")
legend("topright", fill=fcol, legend=rownames(pData(eset)))
dev.off()

## Boxplots
png(file.path(figurePath,"RMA_boxplot.png"),width = 5, height = 5, 
    units = "in", pointsize = 8, bg = "white",res=300)   
par(mar=c(7,4,4,2) + .1)
boxplot(eset, col=fcol, names=phenodata$Sample_Label, main="RMA normalized probeset expression values")
axis(1, at=1:12, labels=F)
text(1:12, par("usr")[3] - 0.25, srt=45, adj = 1, labels=rownames(pData(RawData)), xpd=TRUE)  
dev.off()

## Hierarchical Clustering
data.ALL <- exprs(eset)
clust.cor=hclust(dist(t(data.ALL)),method="complete")
png(file=file.path(figurePath,"RMA_cluster.png"),point=8,width=10,height=7,units="in",res=300)
plclust(clust.cor,main="RMA normalized probeset expression values",sub=paste(nrow(data.ALL), "genes"),xlab="",ylab="Euclidean Distance")
dev.off()

## MDS plots
d <- dist(t(data.ALL),method="euc") # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# plot solution
png(file=file.path(figurePath,"RMA_mds.png"),point=8,width=10,height=7,units="in",res=300)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS\nRMA normalized probeset expression values", type="n", xlim=range(x)*1.2, ylim=range(y)*1.2)
text(x, y, labels = sampleNames(RawData),col= fcol, cex=1, font=2)
#legend(x="topright", legend=sampleNames(RawData), cex=1, text.col=plotcol)
dev.off()

save(eset,file=file.path(rawPath,"normalizedData.RData"))
#################################################################################
### Gene Filtering ###
#################################################################################

# Filter 2 interQuartile Range
iqrs <- apply(exprs(eset),1,function(x)IQR(x))

iqrsFilter.rma <- as.logical(ifelse(iqrs >= 0.5,1,0))
summary(iqrsFilter.rma)

### recored filter stats for future (aka publication) reference
#   Mode   FALSE    TRUE    NA's 
#logical   27891   10598       0 

pdf(file=file.path(figurePath,"IQR_rangefilter.pdf"),bg="white",point=8,width=7,height=7)
plot(quantile(iqrs,probs = seq(0, 1, 0.01)))
abline(h=0.5,col="red")
abline(h=1.0,col="blue")
dev.off()

#clustering using  correlation distance, complete linkage
png(file.path(figurePath, "IQR_filtercluster.png"), width=5, height=5, units='in', bg='white', res=300)
data.ALL <- exprs(eset[iqrsFilter.rma,])
clust.cor=hclust(dist(t(data.ALL)),method="complete")
plclust(clust.cor,main="Cluster Using RMA data post filter",sub=paste(nrow(data.ALL), "genes"),xlab="",ylab="Euclidean Distance")
dev.off()

#################################################################################
### Differential Expression 
#################################################################################

### Completely Randomized Design
CRD <- as.factor(paste(eset$Treatment, sep="."))
design <- model.matrix( ~ -1 + CRD)
colnames(design) <- c("isl2nTf1bMOs","isl2MO","Tf1bMO","UIC")
rownames(design) <- sampleNames(eset)

contrast.matrix <- makeContrasts(
				 "2visl" = isl2nTf1bMOs - isl2MO,
				 "2vTf1" = isl2nTf1bMOs - Tf1bMO,
				 "2vUIC" = isl2nTf1bMOs - UIC,
         "other" = Tf1bMO - UIC,
				  levels=design)

gfit <- lmFit(exprs(eset)[iqrsFilter.rma,],design)
gfit <- contrasts.fit(gfit,contrast.matrix)
gfit <- eBayes(gfit)

gdt <- decideTests2(gfit,adjust.method="BH", p.value=0.05) # at 5% Raw
summary(gdt)

png(file.path(figurePath,"DE_venn_BH.png"),width = 7, height = 5.5, 
    units = "in", pointsize = 8, bg = "white",res=300)
layout(mat=matrix(c(1,2)), heights=c(10,3))
par(mar=c(1,1,1,1))
vennDiagram(gdt)
plot(NA, xlim=c(0,10), ylim=c(0,10), axes=F, xlab="", ylab="", main="Differentially Expressed Genes (BH)")
addtable2plot(x=1,5, table=summary(gdt), display.rownames=T, cex=1.1)
dev.off()

#### Output results tables for adjustments (none, BH, and qvalue)
gdt.none <- decideTests2(gfit,adjust.method="none", p.value=0.05) # at 5% Raw
summary(gdt.none)
write.table(t(c("raw")),file=file.path(tablePath,"testresults.txt"),sep="\t",row.names=F,col.names=F,quote=F)
write.table(summary(gdt.none),file=file.path(tablePath,"testresults.txt"),sep="\t",quote=F,append=T,col.names=NA)

gdt.bh <- decideTests2(gfit,adjust.method="BH", p.value=0.05) # at 5% FDR
summary(gdt.bh)
write.table(t(c("BH")),file=file.path(tablePath,"testresults.txt"),sep="\t",row.names=F,col.names=F,quote=F,append=T)
write.table(summary(gdt.bh),file=file.path(tablePath,"testresults.txt"),sep="\t",quote=F,append=T,col.names=NA)

gdt.qv <- decideTests2(gfit,adjust.method="qvalue", p.value=0.05) # at 5% FDR
summary(gdt.qv)
write.table(t(c("qvalue")),file=file.path(tablePath,"testresults.txt"),sep="\t",row.names=F,col.names=F,quote=F,append=T)
write.table(summary(gdt.qv),file=file.path(tablePath,"testresults.txt"),sep="\t",quote=F,append=T,col.names=NA)

## Final Differential Expression Table
write.table(exprs(eset)[,order(eset$Treatment,eset$Replicate)],file=file.path(tablePath,"expressionData.txt"),sep="\t",col.names=NA)

annotation <- read.table(file.path(designPath,"annotation.txt"),sep="\t",header=T,row.names=1,as.is=T,quote="",comment.char="")
write.fit2(gfit,annotation,exprs(eset)[iqrsFilter.rma,order(eset$Treatment,eset$Replicate)],file=file.path(tablePath,"statisticaltest.txt"), digits = 3, quote=T, adjust = c("none","BH","qvalue"), 
    method = "separate",sep = "\t")

############################################################################
## END ANALYSIS
############################################################################


################################################################################################
## Record R session info
################################################################################################
sessionInfo()
# R version 2.15.2 (2012-10-26)
# Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] pd.090818.zv7.expr_0.0.1 plotrix_3.4-5            Biostrings_2.26.2        IRanges_1.16.4           pdInfoBuilder_1.22.0    
# [6] affxparser_1.30.0        RSQLite_0.11.2           DBI_0.2-5                limma_3.14.3             RColorBrewer_1.0-5      
# [11] oligo_1.22.0             Biobase_2.18.0           oligoClasses_1.20.0      BiocGenerics_0.4.0       qvalue_1.32.0           
# 
# loaded via a namespace (and not attached):
# [1] affyio_1.26.0         BiocInstaller_1.8.3   bit_1.1-9             codetools_0.2-8       ff_2.2-10            
# [6] foreach_1.4.0         GenomicRanges_1.10.5  iterators_1.0.6       KernSmooth_2.23-8     parallel_2.15.2      
# [11] preprocessCore_1.20.0 splines_2.15.2        stats4_2.15.2         tcltk_2.15.2          tools_2.15.2         
# [16] zlibbioc_1.4.0  
