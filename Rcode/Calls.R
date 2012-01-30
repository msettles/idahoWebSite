##############################################################################
#### R script to
####  Calls given groups, these function compute calls that can be used
####	 as filters based upon group membership
####		
#### Written By: Matt Settles
####				Postdoctoral Research Associate
####				Washington State University
####
##############################################################################
####
#### Change Log: 
####	Feb 8, 2008: 
####		formalized code
####		modified pmaCalls, added alpha1 adn alpha2 as parameters
####
##############################################################################
####
####	Usage:
####	source("http://bioinfo-mite.crb.wsu.edu/Rcode/Calls.R")
####
##############################################################################

############################################################
##
##	Function pmaCalls
##
## 	Given an AffyBatch object, gets PMA calls and computes
##	 percentage present calls based
##	on class membership, 
##
## 	returns a list with 4 elements
##	pma - pma calls
##  p.value - computed p values
##  grp - percent present per class
##	alpha1 <- p-value < alpha1 Present 
##	alpha2 <- p-value > alpha2 Absent , otherwise Marginal
##	Note: default in affy mas5calls is (0.04,0.06)
## 		  default in gcos mas5calls    16-20      11 probe pairs/probeset
##								Alpha1  0.04    0.05 
##								Alpha2  0.06   0.065
##
##  result - list of 
##		pma - actual pma calls
##		p.value - pma p value
##		grp - group value
##		result - vector of TRUE/FALSE for filtering
#############################################################	

"plotFilter" <- 
function(eset, allA, allP,type="quantile",filterVal = "99.5%",mainlabel=NULL){
	par(mfrow=c(2,1))
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
		qf <- quantile(eset[allA,],probs = seq(0, 1, 0.001))[filterVal]
    if (type == "sd")
    	qf <- 3*sd(eset[allA,])
    qf
}

"pmaCalls" <- 
function(abatch, class=NULL, percent=1.0,marginalAbsent=TRUE,verbose=TRUE, ...) 
{
    require("affy")
	pmaCalls <- mas5calls(abatch,verbose=verbose, ...)
	calls <- list()
	calls$pma <- assayData(pmaCalls)$exprs
	calls$p.value <- assayData(pmaCalls)$se.exprs

	if ( !is.null(class) ) {
		class <- as.factor(class)
		members <- levels(class)
		result <- rep(FALSE,dim(calls$pma)[1])
		all <- list()
		if (marginalAbsent)
		    disCalls <- matrix(sapply(calls$pma,switch,"A"=0,"M"=0,"P"=1),ncol=dim(calls$pma)[2])
	    else
		    disCalls <- matrix(sapply(calls$pma,switch,"A"=0,"M"=0.5,"P"=1),ncol=dim(calls$pma)[2])
		for (i in members) {
			if(verbose)  print(paste("Checking member", i))
        	sums <- rowSums(as.data.frame(disCalls[, class %in% i]))
        	ok <- sums/ifelse(length(which(class %in% i)) > 0,
        				length(which(class %in% i)),1) 
			if(verbose) print(summary(as.factor(ok)))        	
			all <- cbind(all,ok)
        	result <- result | (ok >= percent)
		}
		names(result) <- rownames(calls$pma)
		calls$result <- result
	
		calls$grp <- apply(all,2,as.numeric)
		colnames(calls$grp) <- members
		rownames(calls$grp) <- rownames(calls$pma)
	}
		
	if(verbose) {print("Summary of results"); print(summary(calls$result))}
	calls
}

############################################################
##
##	Function exprCalls
##
## 	Given an ExprSet object computes Average expression
##	on class membership, 
##
## 	returns a list with 2 elements
##  grp - Average expr per class
##  result - vector of TRUE/FALSE for filtering if at least 
##  	one class average >= Avg
#############################################################	
"exprCalls" <- 
function(eset, class = NULL, percent=1.0,Avg = log2(100), verbose=TRUE) 
{
	if(class(eset) == "ExpressionSet")
		eset <- exprs(eset)
	if(is.matrix(eset)) 
		eset <- as.data.frame(eset)
	
	if (!is.data.frame(eset)) 
        stop("Argument eset is not an ExpressionSet, data.frame or matrix")
		
	calls <- list()
	
	if ( !is.null(class) ) {
		class <- as.factor(class)
		members <- levels(class)
		result <- rep(FALSE,dim(eset)[1])
		all <- list()
		disCalls <- ifelse(eset >= Avg, 1, 0)
	
		for (i in members) {
			tmp <- which( class %in% i )
			if (length(tmp) > 0){
     		    if(verbose) print(paste("Checking member", i))
        	    sums <- rowSums(as.data.frame(disCalls[, class %in% i]))
        		ok <- sums/ifelse(length(which(class %in% i)) > 0,
        							length(which(class %in% i)),1) 
				if(verbose) print(summary(as.factor(ok)))        	
				all <- cbind(all,ok)
        		result <- result | (ok >= percent)
        	}
		}
		names(result) <- rownames(eset)
		calls$result <- result
		calls$pma <- disCalls
		calls$grp <- apply(all,2,as.numeric)
		colnames(calls$grp) <- members
		rownames(calls$grp) <- rownames(eset)
	}
	if(verbose) summary(calls$result)
	calls
}

"GemCalls" <- 
function(eset, cl=NULL) 
{
	uni.cl <- sort(unique(cl2))
    if (length(uni.cl) != 2) stop ("need two classes\n")
	g1 <- which(cl2 == uni.cl[1])
	N1 <- length(cl2[g1])
	N2 <- length(cl2[-g1])
	
	exprSet <- scaled.eset
	#exprSet <- exprs(eset)

   	MEAN1 <- apply(exprSet[,g1],1,mean)
   	MEDIAN1 <- apply(exprSet[,g1],1,median)
	  VAR1 <- apply(exprSet[,g1],1,var)
   	MEAN2 <- apply(exprSet[,-g1],1,mean)
   	MEDIAN2 <- apply(exprSet[,-g1],1,median)
 	  VAR2 <- apply(exprSet[,-g1],1,var)
	
	  ratio <- 2^MEAN1/2^MEAN2
	  RANGE1 <- t(apply(2^exprSet[,g1],1,range))
	  RANGE2 <- t(apply(2^exprSet[,-g1],1,range))
	  overlap <- pmax(RANGE1[,1]/RANGE2[,2],RANGE2[,1]/RANGE1[,2])
   gem <- (ratio >= 2 | ratio <= 0.5) & overlap > 1
	
   num <- MEAN2 - MEAN1
   denom <- sqrt(((N1*VAR1)+(N2*VAR2))/(N1+N2))
	 var0.genes <- which(round(denom, 10) == 0)
	
	  d.stat <- num/denom
	  r.stat <- d.stat/sqrt((d.stat*d.stat)+4)
		
}

