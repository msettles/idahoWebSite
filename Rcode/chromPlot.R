## cPlot.R
plotPoly <- function(chromNum, locs1,locs2, xPoints, chromLens, fg,
                      scale = c("relative","max"),glen=0.4)
{
	## TEST
#	chromNum <- cName
#	locs1 <- -cmeth$start
#	locs2 <- -cmeth$stop 
#	chromLens <- lens
#	fg <- "green" 
	##
	stopifnot(length(locs1) == length(locs2))
    ## Get the scaling factor
    scale <- match.arg(scale)

    scaledX <- cScale(xPoints, chromLens, scale, chromNum)
	
    nlocs <- length(locs1)
    
	
    ## APply the scaling factor to the x positions
    cNum <- match(chromNum, names(chromLens))
    locs1 <- locs1*scaledX
    locs2 <- locs2*scaledX
    mids <- rowMeans(cbind(locs1,locs2))
    if (length(mids) == 0) {
        if (scale == "relative")
            return()
    }
    else {
    	polys <- c()
        ## Determine the direction of the Y plot (+ or -)
        ypos <- rep(cNum, nlocs)
        ytop <- ifelse(locs1>0, ypos+glen, ypos-glen)
        for ( i in 1:length(mids) ){
        	polys<- rbind(polys,cbind(abs(mids[i]),ypos),cbind(abs(mids[i])-1000000*scaledX, ytop),cbind(abs(mids[i])+1000000*scaledX, ytop),cbind(1,NA))}
		#mids, ypos
		#mids-1000, ypos
		#mids+1000, ypos
		
        ## Plot
        polygon(polys, col=fg,border=fg)

        ## Drawn last to ensure that that the lines are actually displayed
    }
    if (scale == "max") {
        lines(c(1,xPoints-1),c(cNum,cNum),col="blue")
    }
    else {
        lines(c(1,max(abs(locs1[!is.na(locs1)]))),c(cNum,cNum),col="blue")
    }
}

# cPlot.R
plotPoly2 <- function(chromNum, locs1,locs2, xPoints,xstart, chromLens, fg,
                      scale = c("relative","max"),glen=0.4)
{
	## TEST
#	chromNum <- cName
#	locs1 <- -cmeth$start
#	locs2 <- -cmeth$stop 
#	chromLens <- lens
#	fg <- "green" 
	##
	stopifnot(length(locs1) == length(locs2))
    ## Get the scaling factor
    scale <- match.arg(scale)

    scaledX <- xPoints/chromLens#cScale(xPoints, chromLens, scale, chromNum)
	
    nlocs <- length(locs1)
    
	
    ## APply the scaling factor to the x positions
    cNum <- match(chromNum, names(chromLens))
    locs1 <- (locs1-xstart)*scaledX
    locs2 <- (locs2-xstart)*scaledX
    mids <- rowMeans(cbind(locs1,locs2))
    if (length(mids) == 0) {
        if (scale == "relative")
            return()
    }
    else {
    	polys <- c()
        ## Determine the direction of the Y plot (+ or -)
        ypos <- rep(1, nlocs)
        ytop <- ypos-(glen/2)
        for ( i in 1:length(mids) ){
        	polys<- rbind(polys,cbind(ypos,abs(mids[i])),cbind(ytop,abs(mids[i])-500000*scaledX),cbind(ytop,abs(mids[i])+500000*scaledX),cbind(1,NA))}
		#mids, ypos
		#mids-1000, ypos
		#mids+1000, ypos
		
        ## Plot
        polygon(polys, col=fg,border=fg)

        ## Drawn last to ensure that that the lines are actually displayed
    }
    if (scale == "max") {
        lines(c(1,1),c(1,xPoints-1),col="blue")
    }
    else {
        lines(c(1,1),c(1,max(abs(locs1[!is.na(locs1)]))),col="blue")
    }
}

plotCluster2 <- function(chromNum, locs1,locs2, Pval, xPoints,xstart, chromLens, fg,
                      scale = c("relative","max"),glen=0.4)
{

	stopifnot(length(locs1) == length(locs2))
    ## Get the scaling factor
    scale <- match.arg(scale)

    scaledX <- xPoints/chromLens
	
    nlocs <- length(locs1)
    ord <- order(clusters$Pval,decreasing=TRUE)
	locs1 <- locs1[ord]
	locs2 <- locs2[ord]
	Pval <- Pval[ord]    

    ## APply the scaling factor to the x positions
    #cNum <- match(chromNum, names(chromLens))
    locs1 <- (locs1-xstart)*scaledX
    locs2 <- (locs2-xstart)*scaledX
    if (length(locs1) == 0) {
        if (scale == "relative")
            return()
    }
    else {
        ## Determine the direction of the Y plot (+ or -)
        ypos <- rep(1, nlocs)
        ytop <- ypos-glen

        ## Plot
        fg <- rep("lightblue4",time=length(locs1))
        fg[which(Pval < 0.01)] <- "blue4"
        rect( ypos,abs(locs1),  ytop,abs(locs2), col="snow1",border=fg)

        ## Drawn last to ensure that that the lines are actually displayed
    }
    if (scale == "max") {
        lines(c(1,xPoints-1),c(cNum,cNum),col="blue")
    }
    else {
        lines(c(1,1),c(1,xPoints),col="blue")
    }
}

plotCluster <- function(chromNum, locs1,locs2, Pval, xPoints, chromLens, fg,
                      scale = c("relative","max"),glen=0.4)
{

	stopifnot(length(locs1) == length(locs2))
    ## Get the scaling factor
    scale <- match.arg(scale)

    scaledX <- cScale(xPoints, chromLens, scale, chromNum)
    nlocs <- length(locs1)
    ord <- order(Pval,decreasing=TRUE)
	locs1 <- locs1[ord]
	locs2 <- locs2[ord]
	Pval <- Pval[ord]    

    ## APply the scaling factor to the x positions
    cNum <- match(chromNum, names(chromLens))
    locs1 <- locs1*scaledX
    locs2 <- locs2*scaledX
    if (length(locs1) == 0) {
        if (scale == "relative")
            return()
    }
    else {
        ## Determine the direction of the Y plot (+ or -)
        ypos <- rep(cNum, nlocs)
        ytop <- ifelse(locs1>0, ypos+glen, ypos-glen)

        ## Plot
        fg <- rep("lightblue4",time=length(locs1))
        fg[which(Pval < 0.01)] <- "blue4"
        rect(abs(locs1), ypos, abs(locs2), ytop, col=fg,border=fg)

        ## Drawn last to ensure that that the lines are actually displayed
    }
    if (scale == "max") {
        lines(c(1,xPoints-1),c(cNum,cNum),col="blue")
    }
    else {
        lines(c(1,max(abs(locs1[!is.na(locs1)]))),c(cNum,cNum),col="blue")
    }
}
## cPlot.R
plotData2 <- function(chromNum, locs, xPoints, chromLens,chromStart, fg,
                      scale = c("relative","max"),glen=0.4)
{
    ## Get the scaling factor
    scale <- match.arg(scale)

    scaledX <- xPoints/chromLens

    nlocs <- length(locs)

    ## APply the scaling factor to the x positions
    cNum <- match(chromNum, names(chromLens))
    locs <- (locs-chromStart)*scaledX
    if (length(locs) == 0) {
        if (scale == "relative")
            return()
    }
    else {
        ## Determine the direction of the Y plot (+ or -)
        ypos <- rep(1, nlocs)
        ytop <- ifelse(locs>0, ypos+glen, ypos-glen)

        ## Plot
        segments(ypos,abs(locs), ytop, abs(locs), col=fg)

        ## Drawn last to ensure that that the lines are actually displayed
    }
    if (scale == "max") {
        lines(c(cNum,cNum),c(1,xPoints-1),col="blue")
    }
    else {
        lines(c(1,1),c(1,xPoints),col="blue")
    }
}

## cPlot.R
.plotData <- function(chromNum, locs, xPoints, chromLens, fg,
                      scale = c("relative","max"),glen=0.4)
{
    ## Get the scaling factor
    scale <- match.arg(scale)

    scaledX <- cScale(xPoints, chromLens, scale, chromNum)

    nlocs <- length(locs)

    ## APply the scaling factor to the x positions
    cNum <- match(chromNum, names(chromLens))
    locs <- locs*scaledX
    if (length(locs) == 0) {
        if (scale == "relative")
            return()
    }
    else {
        ## Determine the direction of the Y plot (+ or -)
        ypos <- rep(cNum, nlocs)
        ytop <- ifelse(locs>0, ypos+glen, ypos-glen)

        ## Plot
        segments(abs(locs), ypos, abs(locs), ytop, col=fg)

        ## Drawn last to ensure that that the lines are actually displayed
    }
    if (scale == "max") {
        lines(c(1,xPoints-1),c(cNum,cNum),col="blue")
    }
    else {
        lines(c(1,max(abs(locs[!is.na(locs)]))),c(cNum,cNum),col="blue")
    }
}

cColor <- function(probes, color, plotChroms,
                   scale=c("relative","max"), glen=0.4) {
    ## Passed a vector of probe names, a color and an instance of a
    ## chromLocation class.  Will recolor the specific probes in the
    ## cPlot created plot to match the specified color.  Scale should
    ## be the same as the scale from cPlot
    scale <- match.arg(scale)
    xPoints <- 1000

    gc <- unlist(mget(probes,env=probesToChrom(plotChroms), ifnotfound=NA))
    gchr <- split(names(gc),gc)

    gchr[["NA"]] <- NULL

    ## Look up the locations of these probes in each chromosome,
    ## plotting any results.
    locList <- chromLocs(plotChroms)
    lens <- chromLengths(plotChroms)
    names(lens) <- chromNames(plotChroms)

    sel <- (names(gchr)%in% names(locList))
    chrToDo <- names(gchr)[sel]
    if(!all(sel))
      warning(sprintf("No locations available for probes on chromosome%s %s",
                      c("", "s")[1+(sum(!sel)>1)],
                      names(gchr)[!sel]))
      
    for (cName in chrToDo) {
        locs <- locList[[cName]][gchr[[cName]]]
        locs <- as.numeric(locs[!is.na(locs)])
        if (length(locs) > 0) {
            .plotData(cName, locs, xPoints, lens,
                      color, scale, glen)
        }
    }
}


cPlot <- function(plotChroms, useChroms=chromNames(plotChroms),
                  scale=c("relative", "max"), fg="white",
                  bg="lightgrey", glen=0.4, xlab="", ylab="Chromosome",
                  main = organism(plotChroms)) {
    ## Passed an instance of a chromLocation class, and the number of
    ## points to represent on the X axis, will utilize that data
    ## to plot a set of genes on their proper chromosome locations.
    scale <- match.arg(scale)

    xPoints <- 1000

    chromNames <- chromNames(plotChroms)
    labs <- chromNames[chromNames %in% useChroms]

    lens <- chromLengths(plotChroms)
    whichLabs <- chromNames %in% labs
    lens <- lens[whichLabs]
    names(lens) <- chromNames[whichLabs]

    ## Build the initial plot structure
    op <- par(bg=bg)
    plot(c(1, xPoints), c(1-glen,length(labs)+glen), type="n", xlab=xlab,
         ylab=ylab, axes=FALSE, las=2, main=main)
    par(op)

    axis(2, c(1:length(labs)), labs, las=2)

    chromLocs <- chromLocs(plotChroms)
    byChroms <- chromLocs[labs]

    for (cName in labs) {
        .plotData(cName,byChroms[[cName]], xPoints,
                  lens, fg, scale,glen);
    }
}



