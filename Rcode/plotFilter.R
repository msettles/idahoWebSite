"plotFilter" <- 
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
