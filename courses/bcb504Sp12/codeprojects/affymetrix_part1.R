# Robison/Benner Zebrafish Selenium Experiment

library(affy)

experimentName = "Selenium"

base <- getwd()
codePath <- "http://bioinfo-mite.ibest.uidaho.edu/Rcode/"
rawPath <- file.path(base,"Data")
figurePath <- file.path(base, "Figures")
tablePath <- file.path(base, "Tables")
celPath <- file.path(base,"CELfiles")

###############################################################
### LOAD RAW DATA ###
###############################################################
print("Load Raw Data")

source(file.path(codePath,"loadAffy.R"))
RawData <- loadAffy( targetsFile = "targets.txt", miameFile = "miame.txt", celPath = celPath, notes = experimentName)
RawData <- RawData[,order(paste(as.character(RawData$Sex),as.character(RawData$Treatment),sep="."))]


#############################################################
### Quality Assessment ###
#############################################################

################
## AffyQCReport
################
print("Quality Assessment")

library(affyQCReport)
try(qaInfo <- affyQAReport( RawData, overwrite = TRUE,repName=experimentName))
	
### AffyRNADeg Table ##
rnaDeg = AffyRNAdeg(RawData)
qaInfo$rnaDeg = rnaDeg
if (exists("qaInfo"))
	save(qaInfo, file=file.path(rawPath, "QAinfo.RData") )

################
## AffyPLM
################
library(affyPLM)
if ( !exists("RawData") )
	load( file=file.path(rawPath,"RawDataABatch.RData" ) )
if ( !exists("qaInfo")){ ## can use fitPLM object from qaInfo 
	load(file=file.path(rawPath, "QAinfo.RData") )
	pset <- qaInfo$affyPLM
} else {
    pset <- fitPLM(RawData, background = FALSE, normalize = FALSE)
}

# weight images
for( i in 0:floor((length(RawData)-1)/8) ){
	png(filename=file.path(base,"affyQA",paste("weight",i,".png",sep="")),width=1024,height=2048, pointsize=10)
	par(mfrow=c(4,2))
	for ( j in 1:8 ){
		if ((i*8+j) <= length(RawData)) image(pset, which = i*8 + j, type="weights")
	}
	dev.off()
}
# residuals
for( i in 0:floor((length(RawData)-1)/8) ){
	png(filename=file.path(base,"affyQA",paste("resid",i,".png",sep="")),width=1024,height=2048, pointsize=10)
	par(mfrow=c(4,2))
	for ( j in 1:8 ){
		if ((i*8+j) <= length(RawData) ) image(pset, which = i*8 + j, type="resids")
	}
	dev.off()
}

############
#### Affy Express QC Report ####
############
library(AffyExpress)
tmp <- getwd()
dir.create(file.path(base,"affyQA","AffyExpressQA"),recursive=TRUE)
setwd(file.path(base,"affyQA","AffyExpressQA"))
AffyQA(parameters=c("Sex","Treatment"),raw=RawData)
setwd(tmp)

### REMOVE 107
#RawData <- RawData[,-which(sampleNames(RawData) == "107")]

############################################################
### Normalization Routines ##
############################################################=
print("Normalization Routines")
## RMA
eSet.rma <- rma(RawData)
save(eSet.rma,file=file.path(rawPath,"rma.RData"))
## PLIER
library(plier)
eSet.plier <- justPlier(RawData,normalize=TRUE) # plier log2 by default no variance stabilizing method
exprs(eSet.plier) <- log2((2^exprs(eSet.plier)+16))
save(eSet.plier,file=file.path(rawPath,"plier.RData"))
## GCRMA
library(gcrma)
eSet.gcrma <- gcrma(RawData)
save(eSet.gcrma,file=file.path(rawPath,"gcrma.RData"))
## MAS5
library(simpleaffy)
eSet.mas <- mas5(RawData,normalize=FALSE)
tgt <- median(apply(exprs(eSet.mas),2,mean,trim=0.02))
sfs <- tgt/apply(exprs(eSet.mas),2,mean,trim=0.02)
eSet.mas <- affy.scalevalue.exprSet(eSet.mas, sc = tgt, analysis = "absolute")
preproc(eSet.mas)$sfs = sfs
preproc(eSet.mas)$tgt = tgt
exprs(eSet.mas)  <- log2(exprs(eSet.mas)+16)
save(eSet.mas,file=file.path(rawPath,"mas.RData"))

## MBEI
eSet.mbei <- expresso(RawData,normalize.method="invariantset",bg.correct=FALSE,
                      pmcorrect.method="pmonly",summary.method="liwong")
exprs(eSet.mbei) <- log2(exprs(eSet.mbei))#
save(eSet.mbei,file=file.path(rawPath,"mbei.RData"))  

#Histograms
pdf(file=file.path(figurePath,"normhistograms.pdf"),width=7,height=10,pointsize=8)
par(mfrow=c(3,1))
plotDensity(log2(pm(RawData)),main="Raw Data")
plotDensity(exprs(eSet.mas),main="MAS5 Lab")
plotDensity(exprs(eSet.rma),main="RMA Lab")
plotDensity(exprs(eSet.gcrma),main="GCRMA Lab")
plotDensity(exprs(eSet.plier),main="PLIER Lab")
plotDensity(exprs(eSet.mbei),main="MBEI Lab")
dev.off()

#boxplots comparing tissue profiles
pdf(file=file.path(figurePath,"normboxplots.pdf"),width=6,height=10,pointsize=8)
par(mfrow = c(5,1))
boxplot(exprs(eSet.mas),main="MAS5 Preprocessing",col=c("red"))
boxplot(exprs(eSet.plier),main="PLIER Preprocessing",col=c("red"))
boxplot(exprs(eSet.rma),main="RMA Preprocessing",col=c("green"))
boxplot(exprs(eSet.gcrma),main="GCRMA Preprocessing",col=c("green"))
boxplot(exprs(eSet.mbei),main="MBEI Preprocessing",col=c("purple"))
dev.off()

# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
library(RColorBrewer)
pdf(file=file.path(figurePath,"normMDS.pdf"),width=7,height=7,pointsize=10)

#RawData
d <- dist(t(log2(pm(RawData))),method="euc") # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results
# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS Raw Data", type="n")
text(x, y, labels = paste(RawData$Sex,RawData$Treatment,RawData$Tank,RawData$ID,sep="."), cex=.7,col= brewer.pal(11,"Paired")[as.numeric(as.factor(paste(RawData$Sex,RawData$Treatment,sep=".")))])

d <- dist(t(exprs(eSet.mas)),method="euc") # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
#fit # view results
# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS MAS5 Data", type="n")
text(x, y, labels = paste(eSet.mas$Sex,eSet.mas$Treatment,eSet.mas$Tank,eSet.mas$ID,sep="."), cex=.7,col= brewer.pal(11,"Paired")[as.numeric(as.factor(paste(eSet.mas$Sex,eSet.mas$Treatment,sep=".")))])

d <- dist(t(exprs(eSet.rma)),method="euc") # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
#fit # view results
# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS RMA Data", type="n")
text(x, y, labels = paste(eSet.rma$Sex,eSet.rma$Treatment,eSet.rma$Tank,eSet.rma$ID,sep="."), cex=.7,col= brewer.pal(11,"Paired")[as.numeric(as.factor(paste(eSet.rma$Sex,eSet.rma$Treatment,sep=".")))])

d <- dist(t(exprs(eSet.gcrma)),method="euc") # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
#fit # view results
# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS GCRMA Data", type="n")
text(x, y, labels = paste(eSet.gcrma$Sex,eSet.gcrma$Treatment,eSet.gcrma$Tank,eSet.gcrma$ID,sep="."), cex=.7,col= brewer.pal(11,"Paired")[as.numeric(as.factor(paste(eSet.gcrma$Sex,eSet.gcrma$Treatment,sep=".")))])

d <- dist(t(exprs(eSet.plier)),method="euc") # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
#fit # view results
# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS PLIER Data", type="n")
text(x, y, labels = paste(eSet.plier$Sex,eSet.plier$Treatment,eSet.plier$Tank,eSet.plier$ID,sep="."), cex=.7,col= brewer.pal(11,"Paired")[as.numeric(as.factor(paste(eSet.plier$Sex,eSet.plier$Treatment,sep=".")))])


d <- dist(t(exprs(eSet.mbei)),method="euc") # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
#fit # view results
# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS MBEI Data", type="n")
text(x, y, labels = paste(eSet.plier$Sex,eSet.plier$Treatment,eSet.plier$Tank,eSet.plier$ID,sep="."), cex=.7,col= brewer.pal(11,"Paired")[as.numeric(as.factor(paste(eSet.mbei$Sex,eSet.mbei$Treatment,sep=".")))])

dev.off()




############################################################
### Get PMA Calls ##
############################################################
print("Get PMA Calls")
if ( !exists("RawData") )
	load( file=file.path(rawPath,"RawDataABatch.RData" ) )

source(file.path(codePath, "Calls.R"))
pma <- pmaCalls(RawData,factor(paste(RawData$Treatment,RawData$Sex,sep=".")),percent=1.0)

summary(pma$result)
#   Mode   FALSE    TRUE    NA's 
#logical    3435   12182       0

allA <- featureNames(RawData)[rowSums(pma$grp) == 0]
length(allA)
# 1367
allP <- featureNames(RawData)[rowSums(pma$grp) == ncol(pma$grp)]
length(allP)
# 10633
save(pma,allA, allP,file=file.path(rawPath,"PMAcalls.RData"))


## Draw Venn Diagrams of each pma
source(file.path(codePath,"Venn.R"))
## too many tissues for venn
pdf(file=file.path(figurePath,"pmaVenns.pdf"))
vennDiagram(floor(pma$grp))
dev.off()
write.table(vennCounts(floor(pma$grp)),file=file.path(tablePath,"PMAVenn.txt"),sep="\t")

#mosaic plot
y <- t(apply(pma$pma,2,table))
pdf(file.path(figurePath,"mosaic.pdf"),point=8,width=10,height=7)
mosaicplot(y,
      main= "PMA call vs Array", ylab = "PMA Call", xlab = "Array")
cex <- 0.9
mtext(y[,1],side = 1,line = -1,at= cumsum(rowSums(y)/sum(y+1-cex)),cex=cex)
mtext(y[,2],side = 1,line = -11,at= cumsum(rowSums(y)/sum(y+1-cex)),cex=cex)
mtext(y[,3],side = 1,line = -20,at= cumsum(rowSums(y)/sum(y+1-cex)),cex=cex)
barplot(t(y))
dev.off()

## Cluster by PMA calls, 
clust.data <- matrix(sapply(pma$pma,switch,"A"=-1,"M"=0,"P"=1),ncol=dim(pma$pma)[2])
colnames(clust.data) <- colnames(pma$pma)
rownames(clust.data) <- rownames(pma$pma)
clust.cor=hclust(dist(t(clust.data),method="manhattan"),method="complete")
pdf(file=file.path(figurePath,"PMAcluster.pdf"),point=8,width=10,height=7)
plclust(clust.cor,main="Cluster PMA Calls",sub="",xlab="",ylab="Manhattan Distance (complete linkage)")
dev.off()

# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
library(RColorBrewer)
d <- dist(t(clust.data),method="manhattan") # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution
pdf(file=file.path(figurePath,"mdsPMA.pdf"),width=10,height=7,pointsize=8)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
  main="Metric MDS", type="n")
text(x, y, labels = paste(RawData$Sex,RawData$Treatment,sep="."), cex=.7,col= brewer.pal(11,"Paired")[as.numeric(as.factor(paste(RawData$Sex,RawData$Treatment,sep=".")))])
dev.off()


##########################################################
### Gene Filtering ###
##########################################################
print("Gene Filtering")
## Filter on PMA Calls
if ( !exists("eSet.rma") ){
	load(file.path(rawPath,"rma.RData"))
	load(file.path(rawPath,"mas.RData"))
	load(file.path(rawPath,"gcrma.RData"))
	load(file.path(rawPath,"plier.RData"))
	load(file.path(rawPath,"mbei.RData"))
}
if ( !exists("pma") )
	load(file.path(rawPath,"PMAcalls.RData"))
source(file.path(codePath, "Calls.R"))

pma <- mas5calls(RawData,alpha1=0.04,alpha2=0.06,TRUE)
## Use allA probes as negControls
allA <-  apply(exprs(pma),1,function(x) all(x=="A"))
summary(allA)
#   Mode   FALSE    TRUE    NA's 
#logical   14462    1155       0 

allP <- apply(exprs(pma),1,function(x) all(x=="P"))
summary(allP)
# Mode   FALSE    TRUE    NA's 
#logical    4984   10633       0 

class <- as.factor(paste(RawData$Sex,RawData$Treatment,sep="."))


plotFilter <- 
function(eset, allA, allP,type="quantile",filterVal = "99.5%",log=FALSE, mainlabel){
	par(mfrow=c(2,1))
	if (log) eset <- log2(eset)
	#plot quantiles
    plot(seq(0,100,0.1),quantile(eset[allA,],probs = seq(0, 1, 0.001)), type="p",col="red",ylim=c(0,16),main=mainlabel,xlab="quantile",ylab="log2 expression")
    points(seq(0,100,0.1),quantile(eset[allP,],probs = seq(0, 1, 0.001)), type="p",col="blue")
    points(seq(0,100,0.1),quantile(eset,probs = seq(0, 1, 0.001)), type="p",col="black")
    abline(h=quantile(eset[allA,],probs = seq(0, 1, 0.001))[filterVal],col="red")
    abline(h=quantile(eset[allP,],probs = seq(0, 1, 0.001))[paste(100-as.numeric(substring(filterVal,1,nchar(filterVal)-1)),"%",sep="")],col="blue")
   #plot histogram
    res <- hist(eset,breaks=100,col="black",xlab="log2 expression",main="Histogram")
    hist(eset[allP,],breaks=res$breaks,col="blue",add=TRUE)
    hist(eset[allA,],breaks=res$breaks,col="red",add=TRUE)
	if (type == "quantile")
		qf <- mean(quantile(eset[allA,],probs = seq(0, 1, 0.001))[filterVal])
    if (type == "sd")
    	qf <- mean(mean(eset[allA,]) + 3*sd(eset[allA,]))
    if (type == "mad")
    	qf <- median(eset[allA,]) + 3*mad(eset[allA,],high=TRUE)
    if(log) qf <- 2^qf
    qf
    
}


# Filter on intensity value
pdf(file=file.path(figurePath,"filterPlots.pdf"),point=8,width=10,height=7)

#RMA
qf <- plotFilter(exprs(eSet.rma),allA,allP, type="quantile", filterVal = "99.5%", mainlabel="RMA Filter plots")
qf
#     99% 
# 6.5036
calls.rma <- exprCalls(eSet.rma,Avg=qf,class=class,percent=1.0)
summary(calls.rma$result)
#   Mode   FALSE    TRUE    NA's 
#logical    6628    8989       0 

#MAS5
qf <- plotFilter(exprs(eSet.mas),allA,allP,type="quantile",mainlabel="MAS5 Filter plots")
qf
#     99% 
# 5.4981
calls.mas <- exprCalls(eSet.mas,Avg=qf,class=class,percent=1.0)
summary(calls.mas$result)
#   Mode   FALSE    TRUE    NA's 
#logical    4727   10890       0 

# GCRMA
qf <- plotFilter(exprs(eSet.gcrma),allA,allP,type="quantile",mainlabel="GCRMA Filter plots")
qf
#     99% 
#  2.1719
calls.gcrma <- exprCalls(eSet.gcrma,Avg=qf,class=class,percent=1.0)
summary(calls.gcrma$result)
#   Mode   FALSE    TRUE    NA's 
#logical    4540   11077       0 

# PLIER
qf <- plotFilter(exprs(eSet.plier),allA,allP,type="quantile",mainlabel="PLIER Filter plots")
qf
#     99% 
#  7.2324
calls.plier <- exprCalls(eSet.plier,Avg=qf,class=class,percent=1.0)
summary(calls.plier$result)
#   Mode   FALSE    TRUE    NA's 
#logical    6932    8685       0
# MBEI
qf <- plotFilter(exprs(eSet.mbei),allA,allP,type="quantile",mainlabel="MBEI Filter plots")
qf
#     99% 
#  7.9443
calls.mbei <- exprCalls(eSet.mbei,Avg=qf,class=class,percent=1.0)
summary(calls.mbei$result)
#   Mode   FALSE    TRUE    NA's 
#logical    7326    8291       0 


dev.off()

filter <-  calls.mas$result & calls.rma$result & calls.gcrma$result & calls.plier$result & calls.mbei$result
summary(filter)
#   Mode   FALSE    TRUE    NA's 
#logical    7695    7922       0 

source(file.path(codePath,"Venn.R"))
pdf(file=file.path(figurePath,"vennFilter.pdf"))
vennDiagram(cbind(calls.mas$result,calls.plier$result,calls.rma$result,calls.gcrma$result),names=c("MAS5","PLIER","RMA","GCRMA"))
dev.off()

# Filter 2 interQuartile Range

#GCRMA
iqrs <- apply(exprs(eSet.gcrma),1,function(x) IQR(x))
plot(quantile(iqrs,probs = seq(0, 1, 0.01)))
iqrsFilter.gcrma <- as.logical(ifelse(iqrs >= 0.5,1,0))
summary(iqrsFilter.gcrma)
# Mode   FALSE    TRUE    NA's 
#logical   13795    1822       0 

####################################
### FOR MAIA LATER
#iqrs <- apply(exprs(eSet.gcrma),1,function(x) diff(range(x)))
#plot(quantile(iqrs,probs = seq(0, 1, 0.01)))
#iqrsFilter.gcrma <- as.logical(ifelse(iqrs >= 1.0,1,0))
#summary(iqrsFilter.gcrma)
#####################################

summary(iqrsFilter.gcrma & calls.gcrma$result)
#   Mode   FALSE    TRUE    NA's 
#logical   13853    1764       0 

#clustering using  correlation distance, complete linkage
pdf(file=file.path(figurePath,"GCRMAcluster.pdf"),point=8,width=10,height=7)
## GCRMA
data.ALL <- exprs(eSet.gcrma[iqrsFilter.gcrma & calls.gcrma$result,])
clust.cor=hclust(dist(t(data.ALL)),method="complete")
plclust(clust.cor,main="Cluster Using GCRMA data post filter",sub=paste(nrow(data.ALL), "genes"),xlab="",ylab="Euclidean Distance")

dev.off()

save(filter,iqrsFilter.gcrma,calls.gcrma,file=file.path(rawPath,"filters.RData"))

###########################################################################################################
print("Differential Expression")
library(limma)
# Global significant 
if (!exists("eSet.gcrma")){
    load(file.path(rawPath,"gcrma.RData"))
    load(file=file.path(rawPath,"filters.RData"))
}
design <- model.matrix( ~ as.numeric(as.character(eSet.gcrma$Weight)) + eSet.gcrma$Sex + eSet.gcrma$Treatment + eSet.gcrma$Sex:eSet.gcrma$Treatment)
rownames(design) <- sampleNames(eSet.gcrma)
colnames(design) <- c("Intercept","Weight","Male", "Selenium","Interaction")

contrast.matrix <- makeContrasts( Sex= Male, Selenium=Selenium, Interaction=Interaction, levels=design)

gfit.gcrma <- lmFit(eSet.gcrma[iqrsFilter.gcrma & calls.gcrma$result,],design)
gfit.gcrma <- contrasts.fit(gfit.gcrma,contrast.matrix)
gfit.gcrma <- eBayes(gfit.gcrma)
## no multiple testing correction
gdt.gcrma <- decideTests(gfit.gcrma,adjust.method="BH", p.value=0.4)
summary(gdt.gcrma) 
#Sex Selenium Interaction
#-1  231        3           0
#0  1243     1760        1764
#1   290        1           0

source(file.path(codePath, "sigTable.R"))
library("zebrafish.db")
write.fit(gfit.gcrma, results=gdt.gcrma, file=file.path(tablePath,"limma.gcrma.fit.txt"), digits=3, adjust="BH", sep="\t")
write.exprs(eSet.gcrma[calls.gcrma$result & iqrsFilter.gcrma,],file=file.path(tablePath,"limma.gcrma.expr.txt"))
write.exprs(eSet.gcrma,file=file.path(tablePath,"limma.gcrma.exprNOFILTER.txt"))

sigdif <- rownames(gdt.gcrma)[which(rowSums(abs(gdt.gcrma))>0)]
affyTable(sigdif,gfit.gcrma,eSet.gcrma,
						filename=file.path(tablePath,"limma.gcrma.annoaffy"),
						tests=c("BH"),anncol=aaf.handler(),annarray="zebrafish.db",
						html=TRUE,excel=TRUE)
						

### CRD  desgin
design <- model.matrix( ~ -1 + as.numeric(as.character(eSet.gcrma$Weight)) + as.factor(paste(eSet.gcrma$Sex,eSet.gcrma$Treatment,sep=":")))
rownames(design) <- sampleNames(eSet.gcrma)
colnames(design) <- c("Weight","FC","FS", "MC","MS")

contrast.matrix <- makeContrasts( FS = (FS-FC), MS = MS-MC,Sex= (MC + MS)/2 - (FC + FS)/2, Selenium=(FS + MS)/2 - (MC + FC)/2,  levels=design)

gfit.gcrma <- lmFit(eSet.gcrma[iqrsFilter.gcrma & calls.gcrma$result,],design)
gfit.gcrma <- contrasts.fit(gfit.gcrma,contrast.matrix)
gfit.gcrma <- eBayes(gfit.gcrma)
## no multiple testing correction
gdt.gcrma <- decideTests(gfit.gcrma,adjust.method="BH", p.value=0.4)
summary(gdt.gcrma) 
#      FS   MS  Sex Selenium
# -1    3    0  401        1
# 0  1760 1763  798     1762
# 1     1    1  565        1

source(file.path(codePath, "sigTable.R"))
library("zebrafish.db")
write.fit(gfit.gcrma, results=gdt.gcrma, file=file.path(tablePath,"limma.gcrma.CRD.fit.txt"), digits=3, adjust="BH", sep="\t")


sigdif <- rownames(gdt.gcrma)[which(rowSums(abs(gdt.gcrma))>0)]
htmlout <- topTable(gfit.gcrma,number=200)$ID
affyTable(sigdif,gfit.gcrma,eSet.gcrma,
						filename=file.path(tablePath,"limma.gcrma.CRD.annoaffy"),
						tests=c("BH"),anncol=aaf.handler(),annarray="zebrafish.db",
						html=TRUE,excel=TRUE)
					
