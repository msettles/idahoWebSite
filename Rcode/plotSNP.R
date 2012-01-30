
pNames <- sfpProbeNames(RawData)

bg.method <- "none"
normalize.method <- "none"
summary <- "medianpolish"
pmcorrect <- "pmonly"

# compute expression estimates
try(norm.object <- bg.correct(RawData,method = bg.method))
if ( normalize.method != "none")
	try(norm.object <- normalize(norm.object, method=normalize.method))
	try(eset.object <- rma(norm.object,subset=probenames, verbose = FALSE, destructive = TRUE, normalize = FALSE, background = FALSE))

# Normalized Object
norm.object <- log(pm(norm.object,pNames),2) # log base2 for fold change in SAM
pcounts <- table(probeNames(abatch,pNames))
# now record and remove the main effect,  save residuals
res.mat <- norm.object -(apply(exprs(eset.object)[pNames,],2,function(x) rep.int(x,times=pcounts[pNames])))

K <- sfpProbePairsCount(RawData)
S <- length(RawData)
P <- length(pNames)

n.pm <- matrix( t(norm.object), nrow = S * K, ncol = P )
colnames(n.pm) <- pNames

r.pm <- matrix( t(res.mat), nrow = S * K, ncol = P )
colnames(r.pm) <- pNames

genoSFP <- read.table("Tables/SFPresultsSMALL2.txt",sep="\t",header=TRUE)
sfpProbeSetName <- t(matrix(unlist(strsplit(rownames(genoSFP),split = "_at")),nr=2)) # split off _at probe+number
sfpProbeSetName[,1] <- paste(sfpProbeSetName[,1],"_at",sep="")
snp <-  tapply(sfpProbeSetName[,2], sfpProbeSetName[,1], function(x) c(x))
probes <- names(snp)

ctype <- as.numeric(RawData$Genotype)
ltype <- as.numeric(RawData$Tissue)
xpos <- 1:K

i <- 1
pdf(file="Figures/snpPlots2.pdf",width=10,height=5,pointsize=8)
for ( i in 1:length(probes) ) {

par(mfrow=c(1,2))
matplot(xpos, t(matrix(n.pm[,probes[i]],nr = S)),
         type = "l", col = ctype, lty = ltype,
         xlab = NA, ylab = "log intensity", main = "")
axis(1,1:K,labels=1:K)
legend(x="topright",levels(RawData$Genotype),col=palette()[unique(RawData$Genotype)],pch=3)

abline(v=snp[[i]])
abline(h=unlist(lapply(split(exprs(eset.object)[probes[i],],ctype),mean)),col=c(1,2))

matplot(xpos, t(matrix(r.pm[,probes[i]],nr = S)),
         type = "l", col = ctype, lty = ltype,
         xlab = NA, ylab = "log intensity", main = "")
axis(1,1:K,labels=1:K)
#legend(x="topright",levels(RawData$Genotype),col=palette()[unique(RawData$Genotype)],pch=3)

abline(v=snp[[i]])

mtext(paste(probes[i], "probe:",paste(sort(snp[[i]]),collapse=", "),sep=" "),side=1,outer=TRUE,padj=-5)

#abline(h=unlist(lapply(split(exprs(eset.object)[probes[i],],ctype),mean)),col=c(1,2))
}
dev.off()

