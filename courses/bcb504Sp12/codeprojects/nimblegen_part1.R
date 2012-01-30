######## PREAMBLE ############
## Investigator: Trish Hartzell
## Bioinformatician: Matt Settles
## Date Created: Oct 16, 2011
## Last Modified: Oct 16,2011
##############################


#################################################################################
### Expects the folders and files to be present
#################################################################################
#UNIX
# Data - folder to dump binary result files
# Figures - folder to dump figures
# Tables - folder to dump tables
# XYS_Files - contains xys (raw) data files 
# 110728_Mxan_MS_til_exp_HX12 - contains design files
#
# targets.txt - experiment description file
#############

library(oligo)
library(pdInfoBuilder)
library(RColorBrewer)
library(limma)
library(qvalue)

sessionInfo()
#R Under development (unstable) (2011-10-10 r57206)
#Platform: x86_64-unknown-linux-gnu (64-bit)

#locale:
# [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C               LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8     LC_MONETARY=en_CA.UTF-8   
# [6] LC_MESSAGES=en_CA.UTF-8    LC_PAPER=C                 LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#[1] stats     graphics  grDevices datasets  utils     methods   base     

#other attached packages:
# [1] qvalue_1.27.0         limma_3.9.21          RColorBrewer_1.0-5    pdInfoBuilder_1.17.0  affxparser_1.25.1     RSQLite_0.10.0       
# [7] DBI_0.2-5             oligo_1.17.3          preprocessCore_1.15.0 oligoClasses_1.15.56  Biobase_2.13.10      

#loaded via a namespace (and not attached):
#[1] affyio_1.21.2      Biostrings_2.21.11 bit_1.1-7          ff_2.2-3           IRanges_1.11.30    splines_2.15.0     tcltk_2.15.0      
#[8] zlibbioc_0.1.8    

## Some of my R code files, should be in the main directory for now
source("Venn.R")
source("functions.R")

base <- getwd() ## expects R was started in the primary experiment directory

arraydesignname <- "110728_Mxan_MS_til_exp_HX12"
experimentName <- "Hartzell M_Xan"

rawPath <- file.path(getwd(),"Data")
figurePath <- file.path(getwd(), "Figures")
tablePath <- file.path(getwd(), "Tables")
xysPath <- file.path(getwd(),"XYS_Files")
designPath <- file.path( getwd(),arraydesignname)



##### Build pdInfoBuilder packages, only need to do once
#
#ndf <- list.files(designPath, pattern = ".ndf$",full.names = TRUE)
#xys <- list.files(xysPath, pattern = ".xys",full.names = TRUE)[1]
#seed <- new("NgsExpressionPDInfoPkgSeed",
#     ndfFile = ndf, xysFile = xys,
#     author = "Matt Settles",
#     email = "msettles@uidaho.edu",
#     biocViews = "AnnotationData",
#     genomebuild = "M_Xan",
#     organism = "Myxococcus Xanthus", species = "Bacteria",
#     url = "http://bioinfo-mite.ibest.uidaho.edu")
#makePdInfoPackage(seed, destDir = ".",unlink=TRUE)

####################################
## Need to be able to complile R packages!!
## Then build and install into R, I work in the R development version Rdev
## for normal R replace Rdev with R
#####################################
#packageName <- dir(pattern="^pd.")
#system(paste("Rdev CMD build ", packageName,sep=""))
#system(paste("Rdev CMD INSTALL ", packageName,"_0.0.1.tar.gz",sep=""))
#################################################################################
### Load data, check quality, 
#################################################################################

targetsFile <- "targets.txt"

# phenodata
phenodata <- read.AnnotatedDataFrame(targetsFile,widget=FALSE,header=TRUE, sep="\t",row.names="Sample_ID")
varMetadata(phenodata)$channel = factor(c("Cy3"),levels=c("_ALL_"))   

phenodata$exp <- substring(sampleNames(phenodata),3)

RawData <- read.xysfiles(file.path(xysPath,phenodata$FileName),
    phenoData=phenodata,sampleNames=sampleNames(phenodata),
    notes=experimentName, verbose=TRUE, checkType=TRUE)   


plotcol <- c(brewer.pal(9,"Set1"),brewer.pal(4,"Set2"))[as.numeric(as.factor(paste(RawData$exp)))]
### Quality Assurance
#plot histogram and boxplot
png(file.path(figurePath,"rawhist.png"),width = 5, height = 5, 
    units = "in", pointsize = 8, bg = "white",res=300) 	
hist(RawData,col=plotcol)
dev.off()

#png(file.path(figurePath,"rawboxplot.png"),width = 5, height = 5, 
#    units = "in", pointsize = 8, bg = "white",res=300) 	
#boxplot(RawData)
#dev.off()

data.ALL <- exprs(RawData)
clust.cor=hclust(dist(t(data.ALL)),method="complete")
plclust(clust.cor,main="",sub="",xlab="",ylab="Euclidean Distance")

#pdf(file=file.path(figurePath,"rawcluster.pdf"),point=8,width=10,height=7)
png(file=file.path(figurePath,"rawcluster.png"),point=8,width=10,height=7,units="in",res=300)
plclust(clust.cor,main="Cluster Using RawData",sub=paste(nrow(data.ALL), "genes"),xlab="",ylab="Euclidean Distance")
dev.off()


d <- dist(t(data.ALL),method="euc") # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim

# plot solution
pdf(file=file.path(figurePath,"mdsRawData.pdf"),width=10,height=7,pointsize=12)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
  main="Metric MDS", type="n")
text(x, y, labels = sampleNames(RawData),col= plotcol)
dev.off()

#################################################################################
### Normalize, check quality, filter 
#################################################################################

eset <- rma(RawData)

png(file.path(figurePath,paste(paste(TISSUE,collapse="."),"rmaboxplot.png",sep="-")),width = 5, height = 5, 
    units = "in", pointsize = 8, bg = "white",res=300) 	
boxplot(eset)
dev.off()

data.ALL <- exprs(eset)
clust.cor=hclust(dist(t(data.ALL)),method="complete")
plclust(clust.cor,main="",sub="",xlab="",ylab="Euclidean Distance")

pdf(file=file.path(figurePath,"rmacluster.pdf"),point=8,width=10,height=7)
plclust(clust.cor,main="Cluster Using RMA normalized data",sub=paste(nrow(data.ALL), "genes"),xlab="",ylab="Euclidean Distance")
dev.off()

d <- dist(t(data.ALL),method="euc") # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim

# plot solution
pdf(file=file.path(figurePath,"mdsRMA.pdf"),width=10,height=7,pointsize=12)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
  main="Metric MDS", type="n")
text(x, y, labels = sampleNames(RawData),col= plotcol)
dev.off()

save(eset,file=file.path(rawPath,"normalizedData.RData"))

ndf <- list.files(designPath, pattern = ".ndf$",full.names = TRUE)
ndf <- read.table(ndf,sep="\t",header=T,as.is=T)

##########################################################
### Gene Filtering ###
##########################################################

# Filter 2 interQuartile Range
#RMA
iqrs <- apply(exprs(eset),1,function(x)IQR(x))

iqrsFilter.rma <- as.logical(ifelse(iqrs >= 0.5,1,0))
summary(iqrsFilter.rma)

#   Mode   FALSE    TRUE    NA's 
#logical    3316    3626       0 


pdf(file=file.path(figurePath,"IQRrangefilter.pdf"),point=8,width=7,height=7)
plot(quantile(iqrs,probs = seq(0, 1, 0.01)))
abline(h=0.5,col="red")
abline(h=1.0,col="blue")
dev.off()

#clustering using  correlation distance, complete linkage
pdf(file=file.path(figurePath,"rmaclusterfilter.pdf"),point=8,width=10,height=7)
data.ALL <- exprs(eset[iqrsFilter.rma,])
clust.cor=hclust(dist(t(data.ALL)),method="complete")
plclust(clust.cor,main="Cluster Using RMA data post filter",sub=paste(nrow(data.ALL), "genes"),xlab="",ylab="Euclidean Distance")
dev.off()

plclust(clust.cor,main="Cluster Using RMA data post filter",sub=paste(nrow(data.ALL), "genes"),xlab="",ylab="Euclidean Distance")

#################################################################################
### Differential Expression 
#################################################################################
####### ANALYSIS #########
eset <- eset[,which(eset$exp != "h2o")]

Slides <- as.factor(eset$Slide)
CRD <- as.factor(eset$exp)

design <- model.matrix( ~ Slides + CRD)
colnames(design) <- c("Intercept", "S508695", "S511836", "delta.DKX",                 
"delta.mas", "delta.masPLUSpSF25", "DK.1622.tan.development", "DK.1622.tan.vegetative",    
"DK.1622.yellow.development", "DK.1622.yellow.vegetative", "MxH.2582.delta.228", "MxH.2583.504.disruption",   
"pkn.14.disruption","tbdr.disruption") 
rownames(design) <- sampleNames(eset)

contrast.matrix <- makeContrasts(
				 "asgB_mutantvY" =  - DK.1622.yellow.vegetative,
				 "DK.1622.tan.vegetativevY" = DK.1622.tan.vegetative - DK.1622.yellow.vegetative,
				 "delta.masvY" = delta.mas - DK.1622.yellow.vegetative,
				 "delta.DKXvY" = delta.DKX - DK.1622.yellow.vegetative,
				 "delta.masPLUSpSF25vY" = delta.masPLUSpSF25 - DK.1622.yellow.vegetative,
				 "DK.1622.tan.developmentvY" = DK.1622.tan.development - DK.1622.yellow.vegetative,
				 "DK.1622.yellow.developmentvY" = DK.1622.yellow.development - DK.1622.yellow.vegetative,
				 "MxH.2582.delta.228vY" = MxH.2582.delta.228 - DK.1622.yellow.vegetative,
				 "MxH.2583.504.disruptionvY" = MxH.2583.504.disruption - DK.1622.yellow.vegetative,
				 "pkn.14.disruptionvY" = pkn.14.disruption - DK.1622.yellow.vegetative,
				 "tbdr.disruptionvY" = delta.mas - DK.1622.yellow.vegetative,
				 levels=design)

gfit <- lmFit(exprs(eset)[iqrsFilter.rma,],design)
gfit <- contrasts.fit(gfit,contrast.matrix)
gfit <- eBayes(gfit)

gdt <- decideTests2(gfit,adjust.method="BH", p.value=0.05) # at 5% Raw
summary(gdt)

pdf("VennDiagrams.pdf",pointsize=8)
vennDiagram(gdt[,c(1,6,7)],names=c("asg","Tan.devel","Yellow.devel"))
vennDiagram(gdt[,c(2,8,4)],names=c("Tan.veg"," MxH.2582.delta.228","delta.DKX"))
vennDiagram(gdt[,c(8,10,11,4)],names=c("MxH.2582.delta.228","pkn.14.disruption","tbdr.disruption","delta.DKX"))
dev.off()

write.table(t(c("raw")),file=file.path(tablePath,"testresults.txt"),sep="\t",row.names=F,col.names=F,quote=F)
gdt.none <- decideTests2(gfit,adjust.method="none", p.value=0.05) # at 5% Raw
summary(gdt.none)
write.table(summary(gdt.none),file=file.path(tablePath,"testresults.txt"),sep="\t",quote=F,append=T,col.names=NA)

write.table(t(c("BH")),file=file.path(tablePath,"testresults.txt"),sep="\t",row.names=F,col.names=F,quote=F,append=T)
gdt.bh <- decideTests2(gfit,adjust.method="BH", p.value=0.05) # at 5% FDR
summary(gdt.bh)
write.table(summary(gdt.bh),file=file.path(tablePath,"testresults.txt"),sep="\t",quote=F,append=T,col.names=NA)

#write.table(t(c(TISSUE,"\tqvalue")),file=file.path(tablePath,"testresults.txt"),sep="\t",row.names=F,col.names=F,quote=F,append=T)
#library(qvalue)
#gdt.qv <- decideTests2(gfit,adjust.method="qvalue", p.value=0.05) # at 5% FDR
#summary(gdt.qv)
#write.table(summary(gdt.qv),file=file.path(tablePath,"testresults.txt"),sep="\t",quote=F,append=T,col.names=NA)

annotation <- read.table("110728_Mxan_MS_til_exp_HX12/annotation.txt",sep="\t",header=T)
anno2 <- read.table("MXAN_GMXORF_ORF_MX_MY_MATT.txt",sep="\t",as.is=T,fill=T,na.strings="-")
annotation <- data.frame(annotation,anno2[match(sub("_","",annotation$ID),anno2$V1),])
rownames(annotation) <- annotation$ID

write.table(data.frame(exprs(eset),annotation=annotation[match(featureNames(eset),annotation$ID),]),
  file=file.path(tablePath,"expressionData.txt"),sep="\t",col.names=TRUE,row.names=TRUE)

expr <- exprs(eset)
write.fit2(gfit,annotation,expr,file=file.path(tablePath,"statisticaltest_BH.txt"), digits = 3, quote=T, adjust = c("none","BH"), method = "separate",sep = "\t")

###### MORE CONTRASTS

contrast.matrix <- makeContrasts(
				 "Yv-V-Yd" = DK.1622.yellow.vegetative - DK.1622.yellow.development,
				 "Tv-V-Td" = DK.1622.tan.vegetative - DK.1622.tan.development,
				 "Yv-V-Tv" = DK.1622.yellow.vegetative - DK.1622.tan.vegetative,
				 "Yd-V-Td" = DK.1622.yellow.development - DK.1622.tan.development,
				 "Yv-Yd_V_Tv-Td" = (DK.1622.yellow.vegetative - DK.1622.yellow.development)-(DK.1622.tan.vegetative - DK.1622.tan.development),
				 "asgB_mutantvYv" =  - DK.1622.yellow.vegetative,
				 "asgB_mutantvTv" =  - DK.1622.tan.vegetative,
				 "mas-V-masPSF25" = delta.mas-delta.masPLUSpSF25,				
				 levels=design)

gfit <- lmFit(exprs(eset)[iqrsFilter.rma,],design)
gfit <- contrasts.fit(gfit,contrast.matrix)
gfit <- eBayes(gfit)

gdt <- decideTests2(gfit,adjust.method="BH", p.value=0.05) # at 5% Raw
summary(gdt)


write.table(t(c("raw")),file=file.path(tablePath,"testresults2.txt"),sep="\t",row.names=F,col.names=F,quote=F)
gdt.none <- decideTests2(gfit,adjust.method="none", p.value=0.05) # at 5% Raw
summary(gdt.none)
write.table(summary(gdt.none),file=file.path(tablePath,"testresults2.txt"),sep="\t",quote=F,append=T,col.names=NA)

write.table(t(c("BH")),file=file.path(tablePath,"testresults2.txt"),sep="\t",row.names=F,col.names=F,quote=F,append=T)
gdt.bh <- decideTests2(gfit,adjust.method="BH", p.value=0.05) # at 5% FDR
summary(gdt.bh)
write.table(summary(gdt.bh),file=file.path(tablePath,"testresults2.txt"),sep="\t",quote=F,append=T,col.names=NA)

#write.table(t(c(TISSUE,"\tqvalue")),file=file.path(tablePath,"testresults.txt"),sep="\t",row.names=F,col.names=F,quote=F,append=T)
#library(qvalue)
#gdt.qv <- decideTests2(gfit,adjust.method="qvalue", p.value=0.05) # at 5% FDR
#summary(gdt.qv)
#write.table(summary(gdt.qv),file=file.path(tablePath,"testresults.txt"),sep="\t",quote=F,append=T,col.names=NA)

expr <- exprs(eset)
write.fit2(gfit,annotation,expr,file=file.path(tablePath,"statisticaltest2_BH.txt"), digits = 3, quote=T, adjust = c("none","BH"), method = "separate",sep = "\t")



