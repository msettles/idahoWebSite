###############################################################################
###
### Generic MAplot - mixture of code from maplot from affy and limma
### Matt Settles
### Bioinformatics Core
### Washington State University, Pullman, WA
### 
### Created July 7, 2008
### 
###############################################################################
##############
### things to add
### marker name on plot for significant markers
##############
## Function paramters
## object - can be one of MAList, RGList, R, A
## G      - if R,G values: G vector
## M      - if M,A values: M vector
## array  - which array in the matrix
## status - character vector giving status of genes
##  values: character vector giving values of 'status' to be highlighted
##          on the plot. Defaults to unique values of 'status'. Ignored
##          if there is no 'status' vector. 
## subset - which genes to show, if NA and length > 10000 random sample is taken
## show.statistics - show median and sd on plot
## span   -
## family.loess - how to compute loess curve
## cex    - point size
## plot.method - type of plot
## add.loess - add the loess curve to the plot
## lwd    - line width
## lty    - line type
## loess.col - color of loess line
## hline  - vector of where to draw horizontal lines
## legend - whether or not to draw a legend for status
## ...    - extra plot parameters

## USAGE
## tt <- read.csv("chrtest.csv",header=TRUE,row.names=1)
## tt[,3] <- 10^(-tt[,3])

## png("wgplot.png",units="in",width=8,height=5,res=300)
## par(las="2",cex=0.6,pch=21,bg="white")
## wgplot(tt,cutoffs = c(1,3, 5, 7, 9),color=palette()[2:5],labels=as.character(1:4))
## title("Whole Genome Associaton Plot of Significance for Chromosomes 1 to 4")
## dev.off()
##

# G = NULL;M = NULL; array=1; col="black";subset = NA;show.statistics = TRUE; span = 2/3; family.loess = "gaussian"; 
# cex=0.3; pch=16; plot.method = c("normal", "smoothScatter", "add"); add.loess = TRUE; lwd = 1; lty = 1; loess.col = "red";hline=NULL; legend=FALSE
# object
# values <- NULL
# status <- NULL

"maplot" <- 
function (object,G = NULL,M = NULL, array=1, status, values, col,subset = NA , 
    show.statistics = TRUE, span = 2/3, family.loess = "gaussian", 
    cex=0.3, pch=16, plot.method = c("normal", "smoothScatter", "add"), 
    add.loess = TRUE, lwd = 1, lty = 1, loess.col = "red",hline=NULL,
	legend=FALSE, ...) 
{
	if (inherits(object, "MAList")){  ## limma MAList
		A <- object$A[,array]
		M <- object$M[,array]
		if (missing(status))  status <- object$genes$Status
	} else if (inherits(object, "RGList")){  ## limma RGList
		A <- as.numeric((log2(object$R[,array]) + log2(object$G[,array]))/2)
		M <- as.numeric((log2(object$R[,array]) - log2(object$G[,array])))
		if (missing(status))  status <- object$genes$Status
	} else if (!is.null(G)){ ## R/G vectors
		A <- (log2(object) + log2(G))/2
		M <- (log2(object) - log2(G))
		if (missing(status))  status <- NULL
	} else if (!is.null(M)){ ## M/A vectors
		A <- object
		if (missing(status))  status <- NULL
	} else {
		stop("Unknown object type")
	}
    plot.method <- match.arg(plot.method)
    fn.call <- list(...)
    if (plot.method == "smoothScatter") {
        require("geneplotter")
        plotmethod <- "smoothScatter"
    } else if (plot.method == "add") {
        plotmethod <- "add"
    } else {
        plotmethod <- "normal"
    }

    if (is.na(subset))	
		subset = sample(1:length(M), min(c(50000, length(M))))

    if (!is.element("ylim", names(fn.call))) {
        yloc <- max(M,na.rm=TRUE)
    } else {
        yloc <- max(fn.call$ylim)
    }
    if (!is.element("xlim", names(fn.call))) {
        xloc <- max(A,na.rm=TRUE)
    } else {
        xloc <- max(fn.call$xlim)
    }

    if (is.null(status) || all(is.na(status))) {

        if (plotmethod == "smoothScatter") {
            smoothScatter(A[subset], M[subset],cex = cex, ...)
        } else if (plotmethod == "add") {
            points(A[subset], M[subset], cex = cex, ...)
        } else {
            plot(A[subset], M[subset], cex = cex,xlab="A",ylab="M",...)
        }
    } else {
        plot(A[subset], M[subset],type="n", cex = cex,xlab="A",ylab="M",...)
 
        if (missing(values)) {
            if (is.null(attr(status, "values"))) 
                values <- names(sort(table(status), decreasing = TRUE))
            else values <- attr(status, "values")
        }
        sel <- !(status %in% values)
        nonhi <- any(sel)
        if (nonhi)
		        points(A[sel], M[sel], cex = cex, ...)
        nvalues <- length(values)
        if (missing(col)) {
            if (is.null(attr(status, "col"))) {
                col <- nonhi + 1:nvalues
				col <- c("black",sample(colors(),length(col)-1))
            }
            else col <- attr(status, "col")
        }
        col <- rep(col, length = nvalues)
        for (i in 1:nvalues) {
            sel <- status == values[i]
            points(A[sel], M[sel],pch=pch,cex=cex,col=col[i], ...)
        }
        if (legend) {
			legend("bottom",horiz=TRUE, legend = values, 
                pch = pch, , col = col, cex = 0.5)
        }
    }

    abline(0, 0, col = "blue")

    if ( add.loess ) {
        aux <- loess(M[subset] ~ A[subset], degree = 1, span = span, 
                         family = family.loess)$fitted
        o <- order(A[subset])
        A <- A[subset][o]
        M <- aux[o]
        o <- which(!duplicated(A))
        lines(approx(A[o], M[o]), col = loess.col, lwd = lwd, 
            lty = lty)
    }
    if ( show.statistics ) {
        sigma <- IQR(M,na.rm=TRUE)
        mean <- median(M,na.rm=TRUE)
        txt <- format(sigma, digits = 3)
        txt2 <- format(mean, digits = 3)
        text(xloc-1, yloc-1, paste(paste("Median:", txt2), paste("IQR:", 
            txt), sep = "\n"), adj = c(1, 1),cex=0.5)
    }
	if ( !is.null(hline) )
		abline(h = hline) 
    invisible()
}


"volcanoplot" <- 
function (M ,P, status, values, col,subset = NA , 
    span = 2/3, show.statistics=FALSE,
    cex, pch, plot.method = c("normal", "add"), 
	hline=NULL, vline=NULL,
	legend=FALSE, ...) 
{
	M <- as.vector(M)
	P <- as.vector(-log10(P))
	if(missing(status)) status <- NULL

    plot.method <- match.arg(plot.method)
    fn.call <- list(...)
    if (plot.method == "add") {
        plotmethod <- "add"
    }
    else {
        plotmethod <- "normal"
    }

	if (is.na(subset))	
		subset = sample(1:length(M), min(c(10000, length(M))))

    if (!is.element("ylim", names(fn.call))) {
        yloc <- max(P)
    } else {
        yloc <- max(fn.call$ylim)
    }
    if (!is.element("xlim", names(fn.call))) {
        xloc <- max(M)
    } else {
        xloc <- max(fn.call$xlim)
    }
    if (missing(pch)) 
        pch = 16
    if (missing(cex)) 
        cex = 0.3
    if (is.null(status) || all(is.na(status))) {

		if (plotmethod == "add") {
            points(M[subset], P[subset], cex = cex, ...)
        } else {
            plot(M[subset], P[subset], cex = cex,xlab="M",ylab="-log10(p-value)",...)
        }
	} else {
        plot(M[subset], P[subset],type="n", cex = cex,xlab="M",ylab="-log10(p-value)",...)
 
        if (missing(values)) {
            if (is.null(attr(status, "values"))) 
                values <- names(sort(table(status), decreasing = TRUE))
            else values <- attr(status, "values")
        }
        sel <- !(status %in% values)
        nonhi <- any(sel)
        if (nonhi)
		        points(M[sel], P[sel], cex = cex, ...)
        nvalues <- length(values)
        if (missing(col)) {
            if (is.null(attr(status, "col"))) {
                col <- nonhi + 1:nvalues
				col <- c("black",sample(colors(),length(col)-1))
            }
            else col <- attr(status, "col")
        }
        col <- rep(col, length = nvalues)
        for (i in 1:nvalues) {
            sel <- status == values[i]
            points(M[sel], P[sel],pch=pch,cex=cex,col=col[i], ...)
        }
        if (legend) {
			legend("topright", legend = values, 
                pch = pch, , col = col, cex = 0.5)
        }
    }

	abline(v= vline, col = "blue")

	if ( !is.null(hline) ){
            abline(h = -log10(hline),col="blue")
            if (show.statistics)
                text(xloc, -log10(hline)+0.1, paste("p value:", hline), adj = c(1, .1))
	} 
    invisible()
}


