
###############################################################
## Produce heatmap based on sample to sample correlation
## 
###############################################################
## last modified Feb 22, 2010
#
#
##
###############################################################

##############################################################################
## PLOT Correlations
## NEEDS WORK
"plotCor" <-
function (data, probenames = NULL , which= c("pm","mm","both"), labels = NULL,
         log=TRUE, reorder=TRUE, lines=NULL,  method = c("pearson", "kendall", "spearman"),
         output=c("screen","jpeg","png","pdf") ,file = "corimage", text = NULL, ...) {

    require("RColorBrewer")
    # determine data type, currently ok for AffyBatch, ExpressionSet and Matrix
    if ( class(data) == "AffyBatch"){
	which <- match.arg(which)
        if (which=="pm")
            abatch <- pm(data, probenames)
        if (which=="mm")
            abatch <- mm(data, probenames)
        if (which=="both")
            abatch <- exprs(data,probenames)
        if(is.null(labels)) labels <- sampleNames(data)
        if(log) abatch <- log2(abatch)
    }else if ( class(data) == "ExpressionSet" ){
        abatch <- exprs(data)
        if (!is.null(probenames)) abatch <- abatch[probenames,]
        if(is.null(labels)) labels <- sampleNames(data)
    }else if ( is.matrix(as.matrix(data)) ){
        abatch <- data
        if (is.numeric(probenames)) abatch <- abatch[probenames,]
        if (is.null(labels)) labels <- colnames(data)
        if(log) abatch <- log2(abatch)
    } else stop("data is of unknown type")

    # compute correlations
    method <- match.arg(method)
    abatch <- cor(abatch,method=method)
    if(reorder){
        order <- hclust(as.dist(1-abatch))$order
        labels <- labels[order]
        abatch <- abatch[order,order]
    }

    output <- match.arg(output)
    switch(output,
	screen = NULL,
	jpeg   = jpeg(file = paste(file,output,sep=".")),
	png    = png(file = paste(file,output,sep=".")),
	pdf    = pdf(file = paste(file,output,sep="."),width = 10, height = 6,pointsize=8))

    par( mai=c(1.5, 1.2, 0.15, 0.3) )
    layout(cbind(1,2), widths=c(4, 1))

    imcol=colorRampPalette(brewer.pal(9, "RdPu"))(256)
    rg=range(abatch, na.rm=TRUE)
    image(1:ncol(abatch), 1:ncol(abatch), abatch, xlab="", ylab="", zlim=rg,
        main="", col=imcol, axes=FALSE)
    axis(1, at=1:ncol(abatch), labels=labels, las=3)
    axis(2, at=1:ncol(abatch), labels=labels, las=2)
    if(!is.null(lines)){
        abline(h=lines)
        abline(v=lines)
    }
    mtext (text, 3)
    image(1, seq(rg[1], rg[2], length=length(imcol)), rbind(seq_along(imcol)),
        xaxt="n", xlab="", ylab="", col=imcol)
    text(1, 0, "Array Correlations", xpd=NA)
    if (output != "screen") dev.off()
}
