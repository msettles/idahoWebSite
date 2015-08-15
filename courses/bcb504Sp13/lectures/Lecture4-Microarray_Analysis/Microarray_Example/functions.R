
### test
#fit <- gfit
#adjust <- c("none","BH","qvalue") 
#results.cutoff = 0.05
#results.lfc=0
###
write.fit2 <- 
function (fit,anno=NULL,expr=NULL,file,  results.cutoff = 0.05,results.lfc=0, digits = 3, adjust = "none", 
    method = "separate", sep = "\t", ...) 
{
    if (!is(fit, "MArrayLM")) 
        stop("fit should be an MArrayLM object")
    if (is.null(fit$t) || is.null(fit$p.value)) 
        fit <- eBayes(fit)
    method <- match.arg(method, c("separate", "global"))
    p.value <- as.matrix(fit$p.value)
    lad <- length(adjust)
    p.value.adj <- list()
    results.adj <- list()
    for (adjust.method in adjust){
	    if (adjust.method == "none") {
		p.value.adj[[adjust.method]] <- p.value
		results.adj[[adjust.method]] <- unclass(decideTests2(fit,method=method,adjust.method=adjust.method,p.value=results.cutoff,lfc=results.lfc))
	    } else if (adjust.method == "qvalue") {
		require("qvalue")
		p.value.adj[[adjust.method]] <- p.value
		if (method == "separate") 
		    for (j in 1:ncol(p.value)) p.value.adj[[adjust.method]][, j] <- qvalue(p.value[,j])$qvalues
		if (method == "global") 
		    p.value.adj[[adjust.method]] <- qvalue(p.value)$qvalues
		results.adj[[adjust.method]] <- unclass(decideTests2(fit,method=method,adjust.method=adjust.method,p.value=results.cutoff,lfc=results.lfc))
	    } else {
		p.value.adj[[adjust.method]] <- p.value
		if (method == "separate") 
		    for (j in 1:ncol(p.value)) p.value.adj[[adjust.method]][, j] <- p.adjust(p.value[,j], method = adjust.method)
		if (method == "global") 
		    p.value.adj[[adjust.method]] <- p.adjust(p.value, method = adjust.method)
		results.adj[[adjust.method]] <- unclass(decideTests2(fit,method=method,adjust.method=adjust.method,p.value=results.cutoff,lfc=results.lfc))
	    }
    }
   rn <- function(x, digits = digits) if (is.null(x)) 
        NULL
    else {
        if (is.matrix(x) && ncol(x) == 1) 
            x <- x[, 1]
        round(x, digits = digits)
    }
    tab <- list()
    tab$A             <- rn(fit$Amean, digits = digits - 1)
    tab$Coef          <- rn(fit$coef, digits = digits)
    tab$p.value.adj   <- rn(data.frame(p.value.adj), digits = digits + 6)
    tab$F.p.value     <- rn(fit$F.p.value, digits= digits + 6)
    tab$results <- data.frame(results.adj)
    tab <- data.frame(tab, check.names = FALSE)
    if (!is.null(expr)) tab <- data.frame(tab,expr[match(fit$genes[[1]],rownames(expr)),])
    tab$Genes <- fit$genes[[1]]
    if (!is.null(anno)) tab <- data.frame(tab,anno[match(fit$genes[[1]],rownames(anno)),])     
    write.table(tab, file = file,  row.names = FALSE, sep = sep, ...)
}

decideTests2 <- 
function (object, method = "separate", adjust.method = "BH", 
    p.value = 0.05, lfc = 0) 
{
    if (!is(object, "MArrayLM")) 
        stop("Need MArrayLM object")
    if (is.null(object$t)) 
        object <- eBayes(object)
    method <- match.arg(method, c("separate", "global", "hierarchical", 
        "nestedF"))
    adjust.method <- match.arg(adjust.method, c("none", "bonferroni", 
        "holm", "BH", "fdr", "BY", "qvalue"))
    if (adjust.method == "fdr") 
        adjust.method <- "BH"
    switch(method, separate = {
        p <- as.matrix(object$p.value)
        tstat <- as.matrix(object$t)
	if (adjust.method == "qvalue"){	
		for (j in 1:ncol(p)) {
		    o <- !is.na(p[, j])
		    p[o, j] <- qvalue(p[o, j])$qvalues
		}
	} else{
		for (j in 1:ncol(p)) {
		    o <- !is.na(p[, j])
		    p[o, j] <- p.adjust(p[o, j], method = adjust.method)
		}
	}
        results <- new("TestResults", sign(tstat) * (p < p.value))
    }, global = {
        p <- as.matrix(object$p.value)
        tstat <- as.matrix(object$t)
        o <- !is.na(p)
	if (adjust.method == "qvalue")
		p[o] <- qvalue(p[o])$pvalues
	else
	        p[o] <- p.adjust(p[o], method = adjust.method)
        results <- new("TestResults", sign(tstat) * (p < p.value))
    })
    if (lfc > 0) {
        if (is.null(object$coefficients)) 
            warning("lfc ignored because coefficients not found")
        else results@.Data <- results@.Data * (abs(object$coefficients) > 
            lfc)
    }
    results
}


