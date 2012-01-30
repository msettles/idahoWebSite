
###########modification to gplots barplot2 for anova results

### debug
#	height = ests
#	groups = groups
#	width = 1 
#	space = NULL 
#	names.ind = NULL
#	legend.text = NULL
#	beside = FALSE
#	horiz = FALSE 
#	density = NULL 
#   angle = 45
#	col = NULL 
#	prcol = NULL 
#	border = par("fg") 
#   main = NULL
#	sub = NULL 
#	xlab = NULL 
#	ylab = NULL 
#	xlim = NULL 
#   ylim = NULL 
#	xpd = TRUE 
#	log = ""
#	axes = TRUE 
#	axisnames = TRUE 
#   cex.axis = par("cex.axis") 
#	cex.names = par("cex.axis") 
#    inside = TRUE
#	plot = TRUE
#	axis.lty = 0 
#	offset = 0
#	plot.ci = FALSE 
#   ci.l = NULL
#	ci.u = NULL 
#	ci.color = "black" 
#	ci.lty = "solid"
#   ci.lwd = 1
#	plot.grid = FALSE 
#	grid.inc = NULL
#	grid.lty = "dotted" 
#   grid.lwd = 1
#	grid.col = "black" 
#	add = FALSE
#	panel.first = NULL 
#   panel.last = NULL
####


"barplot.anova" <-
function (
	height,
	groups= NULL,
	width = 1, 
	space = NULL, 
	names.ind = NULL,
	names.grp = NULL,
	plot.pv=FALSE,
	pvalues = NULL,
	pv.color="red",
	legend.text = NULL,
	beside = FALSE, 
	horiz = FALSE, 
	density = NULL, 
    angle = 45, 
	col = NULL, 
	prcol = NULL, 
	border = par("fg"), 
    main = NULL, 
	sub = NULL, 
	xlab = NULL, 
	ylab = NULL, 
	xlim = NULL, 
    ylim = NULL, 
	xpd = TRUE, 
	log = "", 
	axes = TRUE, 
	axisnames = TRUE, 
    cex.axis = par("cex.axis"), 
	cex.names = par("cex.axis"), 
    inside = TRUE, 
	plot = TRUE, 
	axis.lty = 0, 
	offset = 0, 
	plot.ci = FALSE, 
    ci.l = NULL, 
	ci.u = NULL, 
	ci.color = "black", 
	ci.lty = "solid", 
    ci.lwd = 1, 
	plot.grid = FALSE, 
	grid.inc = NULL, 
	grid.lty = "dotted", 
    grid.lwd = 1, 
	grid.col = "black", 
	add = FALSE, 
	panel.first = NULL, 
    panel.last = NULL, ...) 
{
	## check args
    if (!missing(inside)) 
        .NotYetUsed("inside", error = FALSE)
    if (missing(space))
        space <- if ((is.matrix(height) && beside) || length(unique(groups)) > 1) c(0, 1) else 0.2
    space <- space * mean(width)
    if (plot && axisnames && missing(names.ind)) 
        names.ind <- if (is.matrix(height)) colnames(height) else names(height)

### need to mess with this
    if (is.vector(height) || (is.array(height) && (length(dim(height)) == 1))) {
        height <- cbind(height)
        beside <- TRUE
        if (is.null(col)) 
            col <- "grey"
    } else if (is.matrix(height)) {
        if (is.null(col)) 
            col <- heat.colors(nrow(height))
    } else stop(paste(sQuote("height"), "must be a vector or a matrix"))
    if (is.logical(legend.text)) 
        legend.text <- if (legend.text && is.matrix(height)) rownames(height)
###
	### log
    logx <- FALSE
    logy <- FALSE
    if (log != "") {
        if (any(grep("x", log))) 
            logx <- TRUE
        if (any(grep("y", log))) 
            logy <- TRUE
    }
    if ((logx || logy) && !is.null(density)) 
        stop("Cannot use shading lines in bars when log scale is used")
	### log
    NR <- nrow(height)
    NC <- ncol(height)

    if (beside) {
        if (length(space) == 2) 
			if(NC > 1){
	            space <- rep.int(c(space[2], rep.int(space[1], NR - 1)), NC)
			} else if (!missing(groups))
				space <- unlist(sapply(table(groups),function(x) c(space[2], rep.int(space[1],x - 1)),simplify=TRUE))
        width <- rep(width, length.out = NR)
    } else width <- rep(width, length.out = NC)

    offset <- rep(as.vector(offset), length.out = length(width))
    delta <- width/2
    w.r <- cumsum(space + width)
    w.m <- w.r - delta
    w.l <- w.m - delta
    if (!beside && (NR > 1) && plot.ci) 
        plot.ci = FALSE
    if (plot && plot.ci) {
        if ((missing(ci.l)) || (missing(ci.u))) 
            stop("confidence interval values are missing")
        if (is.vector(ci.l) || (is.array(ci.l) && (length(dim(ci.l)) == 1))) 
            ci.l <- cbind(ci.l)
        else if (!is.matrix(ci.l)) 
            stop(paste(sQuote("ci.l"), "must be a vector or a matrix"))
        if (is.vector(ci.u) || (is.array(ci.u) && (length(dim(ci.u)) == 1))) 
            ci.u <- cbind(ci.u)
        else if (!is.matrix(ci.u)) 
            stop(paste(sQuote("ci.u"), "must be a vector or a matrix"))
        if (any(dim(height) != dim(ci.u))) 
            stop(paste(sQuote("height"), "and", sQuote("ci.u"), 
                "must have the same dimensions."))
        else if (any(dim(height) != dim(ci.l))) 
            stop(paste(sQuote("height"), "and", sQuote("ci.l"), 
                "must have the same dimensions."))
    }
    if ((logx && horiz) || (logy && !horiz)) {
        height.na <- sum(is.na(height))
        if (height.na > 0) {
            warning(sprintf("%.0f values == NA in 'height' omitted from logarithmic plot", 
                height.na), domain = NA)
        }
        height.lte0 <- sum(height <= 0, na.rm = TRUE)
        if (height.lte0 > 0) {
            warning(sprintf("%0.f values <=0 in 'height' omitted from logarithmic plot", 
                height.lte0), domain = NA)
            if (beside) 
                height[height <= 0] <- NA
        }
        if (plot.ci && (min(ci.l) <= 0)) 
            stop("log scale error: at least one lower c.i. value <= 0")
        if (logx && !is.null(xlim) && (xlim[1] <= 0)) 
            stop("log scale error: 'xlim[1]' <= 0")
        if (logy && !is.null(ylim) && (ylim[1] <= 0)) 
            stop("'log scale error: 'ylim[1]' <= 0")
        if (plot.ci) {
            rectbase <- c(height[is.finite(height)], ci.l)
            rectbase <- min(0.9 * rectbase[rectbase > 0])
        }
        else {
            rectbase <- height[is.finite(height)]
            rectbase <- min(0.9 * rectbase[rectbase > 0])
        }
        if (logy && !is.null(ylim) && !horiz) 
            rectbase <- ylim[1]
        else if (logx && !is.null(xlim) && horiz) 
            rectbase <- xlim[1]
        if (!beside) 
            height <- rbind(rectbase, apply(height, 2, cumsum))
        lim <- if (plot.ci) 
            c(height, ci.l, ci.u)
        else height
        rangeadj <- c(0.9 * lim + offset, lim + offset)
        rangeadj <- rangeadj[rangeadj > 0]
    } else {
        rectbase <- 0
        if (!beside) 
            height <- rbind(rectbase, apply(height, 2, cumsum))
        lim <- if (plot.ci) c(height, ci.l, ci.u) else height
		rangeadj <- c(-0.01 * lim + offset, lim + offset)
    }
    if (horiz) {
        if (missing(xlim)) 
            xlim <- range(rangeadj, na.rm = TRUE)
		if(plot.pv) xlim[2] <- xlim[2]+xlim[2]*0.06
        if (missing(ylim)) 
            ylim <- c(min(w.l), max(w.r))
    } else {
        if (missing(xlim)) 
            xlim <- c(min(w.l), max(w.r))
        if (missing(ylim)) 
            ylim <- range(rangeadj, na.rm = TRUE)
		if(plot.pv) ylim[2] <- ylim[2]+ylim[2]*0.06
    }
    if (beside) 
        w.m <- matrix(w.m, nc = NC)
    if (plot) {
        opar <- if (horiz) par(xaxs = "i", xpd = xpd) else par(yaxs = "i", xpd = xpd)
        on.exit(par(opar))
        if (!add) {
            plot.new()
            plot.window(xlim, ylim, log = log, ...)
        }
        panel.first
        usr <- par("usr")
        if (logx) {
            usr[1] <- 10^usr[1]
            usr[2] <- 10^usr[2]
        }
        if (logy) {
            usr[3] <- 10^usr[3]
            usr[4] <- 10^usr[4]
        }
        if (!missing(prcol)) 
            rect(usr[1], usr[3], usr[2], usr[4], col = prcol)
        if (plot.grid) {
            par(xpd = FALSE)
            if (is.null(grid.inc)) {
                if (horiz) {
                  grid <- axTicks(1)
                  abline(v = grid, lty = grid.lty, lwd = grid.lwd, 
                    col = grid.col)
                } else {
                  grid <- axTicks(2)
                  abline(h = grid, lty = grid.lty, lwd = grid.lwd, 
                    col = grid.col)
                }
            } else {
                if (horiz) {
                  grid <- pretty(xlim, n = grid.inc)
                  abline(v = grid, lty = grid.lty, lwd = grid.lwd, 
                    col = grid.col)
                } else {
                  grid <- pretty(ylim, n = grid.inc)
                  abline(h = grid, lty = grid.lty, lwd = grid.lwd, 
                    col = grid.col)
                }
            }
            par(xpd = xpd)
        }
        xyrect <- function(x1, y1, x2, y2, horizontal = TRUE, 
            ...) {
            if (horizontal) { 
                rect(x1, y1, x2, y2, ...)
            } else rect(y1, x1, y2, x2, ...)
        }
        if (beside) { 
            xyrect(rectbase + offset, w.l, c(height) + offset, 
                w.r, horizontal = horiz, angle = angle, density = density, 
                col = col, border = border)
        } else {
            for (i in 1:NC) xyrect(height[1:NR, i] + offset[i], 
                w.l[i], height[-1, i] + offset[i], w.r[i], horizontal = horiz, 
                angle = angle, density = density, col = col, 
                border = border)
        }
        panel.last
        if (plot.ci) {
            ci.width = width/4
            if (horiz) {
                segments(ci.l, w.m, ci.u, w.m, col = ci.color, 
                  lty = ci.lty, lwd = ci.lwd)
                segments(ci.l, w.m - ci.width, ci.l, w.m + ci.width, 
                  col = ci.color, lty = ci.lty, lwd = ci.lwd)
                segments(ci.u, w.m - ci.width, ci.u, w.m + ci.width, 
                  col = ci.color, lty = ci.lty, lwd = ci.lwd)
            } else {
                segments(w.m, ci.l, w.m, ci.u, col = ci.color, 
                  lty = ci.lty, lwd = ci.lwd)
                segments(w.m - ci.width, ci.l, w.m + ci.width, 
                  ci.l, col = ci.color, lty = ci.lty, lwd = ci.lwd)
                segments(w.m - ci.width, ci.u, w.m + ci.width, 
                  ci.u, col = ci.color, lty = ci.lty, lwd = ci.lwd)
            }
        }
        if (axisnames && (!is.null(names.ind) || !is.null(names.grp))) {
			if (!is.null(names.ind)){	         
			   at.l <- if (length(names.ind) != length(w.m)) {
	                		if (length(names.ind) == NC) 
	                  			colMeans(w.m)
							else if (!missing(groups))
								tapply(w.m,groups,mean)
	                		else stop("incorrect number of names")
	            		} else w.m
	            axis(if (horiz) 2 else 1, 
					at = at.l, labels = names.ind, lty = axis.lty, line=-1,
	                cex.axis = cex.axis, ...)
			}
			if (!is.null(names.grp)){
		        at.l <- if (length(names.grp) != length(w.m)) {
		            		if (length(names.grp) == NC) 
		              			colMeans(w.m)
							else if (!missing(groups))
								tapply(w.m,groups,mean)
		            		else stop("incorrect number of names")
		        		} else w.m
				mtext(side = 1, at = at.l, line = 2, 
				      text = names.grp, col = "black",cex=cex.names,...)
			}
        }
        if (!is.null(legend.text)) {
            legend.col <- rep(col, length = length(legend.text))
            if ((horiz & beside) || (!horiz & !beside)) {
                legend.text <- rev(legend.text)
                legend.col <- rev(legend.col)
                density <- rev(density)
                angle <- rev(angle)
            }
            if (logx) 
                legx <- usr[2] - ((usr[2] - usr[1])/10)
            else legx <- usr[2] - xinch(0.1)
            if (logy) 
                legy <- usr[4] - ((usr[4] - usr[3])/10)
            else legy <- usr[4] - yinch(0.1)
            legend(legx, legy, legend = legend.text, angle = angle, 
                density = density, fill = legend.col, xjust = 1, 
                yjust = 1)
        }
        title(main = main, sub = sub, xlab = xlab, ylab = ylab, 
            ...)
        if (axes) {
            if (plot.grid){ 
                axis(ifelse(horiz,1,2), at = grid, cex.axis = cex.axis, ...)
            } else axis(ifelse(horiz,1,2), cex.axis = cex.axis, ...)
        }
		if (plot.pv) {
			if(!missing(pvalues)){
				Signif <- symnum(pvalues, corr = FALSE, na = FALSE, 
    	               cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
    	               symbols = c("***", "**", "*", ".", " "))
				Signif[1] <- " "
				text(w.m, ci.u+ci.u*0.03, labels=Signif,col = pv.color,cex=1.0)
			} else 
				error("Missing pvalues")
		}
		box()
        invisible(w.m)
    } else w.m
}


