#############################################################
## R Script to Plot Feature Maps onto Genes or Chromosomes ##
#############################################################
# Author: Thomas Girke
# Usage:
# source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/featureMap.txt")  

# Generate Sample Data
startx <- abs(rnorm(400, 0, sd=2000))
endx <- startx+abs(rnorm(400, 0, sd=150))
cat("\nAssign gene or chromosome length to 'genelength', e.g.: \n genelength <- max(endx)\n")
genelength <- max(endx)
strand <- sign(rnorm(400, 0, sd=2000))
mapDF <- data.frame(startx, endx, strand)
cat("\nGenerate input data frame of this format: \n")
print(mapDF[1:4,])
cat("...\n")

# Function to transform DF to proper format
transformDF <- function(plotDF=mapDF, genelength=genelength, method="simple") {
	# Simple and fast arrangement of fragments
	if(method=="simple") {
	# Identify overlaps to plot them with offset
		plotDF <- plotDF[order(-plotDF$strand, plotDF$startx),]
		plotDF <- data.frame(plotDF, OL=sign(c(1,c(plotDF$startx[-1],1)-plotDF$endx)[-length(plotDF[,1])]))
		plotDF[plotDF[,3]==-1 & plotDF[,4]==-1,][1,4] <- 1
		group <- c(which(plotDF[,3]==1 & plotDF[,4]==1), which(plotDF[,3]==-1 & plotDF[,4]==1))
		count <- c(group[-1], length(plotDF[,1]))-c(group[-length(group)],group[length(group)]-1)
		group <- rep(1:length(group), count)
		plotDF <- data.frame(plotDF, group)
		offset <- NULL
		for(i in unique(group)) {
			offset <- c(offset, seq(0.1, length(group[group==i]), by=0.1)[1:length(group[group==i])])
		}
		plotDF <- data.frame(plotDF, offset)
		plotDF[,6] <- plotDF[,3] * plotDF[,6] 
		plotDF <- plotDF[,c(1,6,2,6)]
		names(plotDF)[c(2,4)] <- c("starty", "endy")
	}
	# Arrange fragments more space efficiently in a line-wise manner
	if(method=="arrange") {
		collectDF <- data.frame(startx=0, endx=0, strand=0, count=0)
		cat("\n Please wait. This step can take ~10 sec per 500 fragments!\n")
		mapDF <- plotDF
		for(i in rev(sort(unique(mapDF$strand)))) {
			tempDF <- mapDF[mapDF$strand==i,]
			count <- i
			myiteration <- 0
			tempDF <- tempDF[order(-tempDF$strand, tempDF$startx),]
			while(length(tempDF[,1])>0) {
				myiteration <- myiteration + 1
				if(count==-1 & myiteration==1) {
					previousEnd <- 0
					} else {
					previousEnd <- collectDF$endx[length(collectDF[,1])]
				}
				Condnextline <- which(tempDF$startx > previousEnd)[1]
				if(is.na(Condnextline)) { 
					count <- count + i
					mynext <- 1 
					} else {
					mynext <- Condnextline 
				}
				collectDF <- rbind(collectDF, cbind(tempDF[mynext,], count))
				tempDF <- tempDF[-c(mynext),]
				tempDF <- tempDF[order(-tempDF$strand, tempDF$startx),]
#				cat("\ndim tempDF: ", dim(tempDF), "\n") # debugging
#				cat("\ndim collectDF: ", dim(collectDF), "\n") # debugging
			}
		}
		plotDF <- collectDF[-1,]
		plotDF <- plotDF[, c(1,4,2,4)]
		names(plotDF)[c(2,4)] <- c("starty","endy")
		plotDF[,2] <- plotDF[,2]/10
		plotDF[,4] <- plotDF[,4]/10
	}
	# Scale feature position to x-axis (1:11)
	plotDF[,1] <- (plotDF$startx*10/genelength) + 1
	plotDF[,3] <- (plotDF$endx*10/genelength) + 1
	plotDF
}
cat("\nUse the function transformDF() like this to transform the input mapping coordinates for plotting: \n plotDF <- transformDF(plotDF=mapDF, genelength, method=\"arrange\")\n", "\tThe argument method=\"arrange\" arranges the fragments more space efficiently than method=\"simple\".\n")
plotDF <- transformDF(plotDF=mapDF, genelength=genelength, method="arrange")
cat("\nThe result should have this format:\n")
print(plotDF[1:4,])
cat("...\n")

# Plot mapping frame 'plotDF'
Mymain <- "My Feature Map"
Mysub <- paste("Gene length", as.integer(genelength), "bp")
featurePlot <- function(plotDF=plotDF, genelength=11, geneline=6, genecol="green", featureline=2, featurecol="red", showxy=F) {
	# This first step creates an empty plotting device. To show x-y axis, leave out the arguments: xaxt="n", yaxt="n", bty="n"
	if(showxy==TRUE) { 
		plot(x=1, y=1, type="h", col="white", xlim=c(0,12), ylim=c(-5,5), xlab="", ylab="", main=Mymain, sub=Mysub)
	} else {
		plot(x=1, y=1, type="h", col="white", xlim=c(0,12), ylim=c(-5,5), xlab="", ylab="", main=Mymain, sub=Mysub, xaxt="n", yaxt="n", bty="n")
	}	
			
	segments(1.0, 0.0, genelength, 0.0, col=genecol, lwd=geneline)
	# arrows(1.0, 0.0, 11, 0.0, col="green", lwd=6)
	segments(plotDF$startx, plotDF$starty, plotDF$endx, plotDF$endy, col=featurecol, lwd=featureline)
	text(0.5, 0, "5'"); text(11.5, 0, "3'")
}
cat("\nThis command plots the mapping frame plotDF: \n featurePlot(plotDF=plotDF, geneline=6, genecol=\"green\", featureline=2, featurecol=sample(c(2,4,8), length(plotDF[,1]), replace=T), showxy=F)\n")
featurePlot(plotDF=plotDF, geneline=6, genecol="green", featureline=2, featurecol=sample(c(2,4,8), length(plotDF[,1]), replace=T))
cat("\tThe setting 'showxy=T' will include the xy-axes in the plot\n")
cat("\nThe main title and subtitle can be changed by assigning them to 'Mymain' and 'Mysub'. The defaults are:\n Mymain <- \"My Feature Map\"; Mysub <- paste(\"Gene length\", as.integer(genelength), \"bp\")\n")
