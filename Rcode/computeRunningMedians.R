computeRM <- 
function (xSet, probeAnno, modColumn = "Cy5", allChr = c(1:19, 
    "X", "Y"), winHalfSize = 400, min.probes = 5, quant = 0.5, 
    combineReplicates = FALSE, checkUnique = TRUE, uniqueCodes = c(0), 
    verbose = TRUE) 
{
    stopifnot(inherits(xSet, "ExpressionSet"), all(is.character(allChr)), 
        is.numeric(quant), (quant >= 0) & (quant <= 1), length(quant) == 
            1)
    if (combineReplicates) 
        grouping <- factor(pData(xSet)[[modColumn]])
    else grouping <- factor(sampleNames(xSet))
    newExprs <- matrix(NA, nrow = nrow(exprs(xSet)), ncol = nlevels(grouping))
    rownames(newExprs) <- featureNames(xSet)
    for (chr in allChr) {
        if (verbose) 
            cat("\nChromosome", chr, "...\n")
        chrsta <- get(paste(chr, "start", sep = "."), env = probeAnno)
        chrend <- get(paste(chr, "end", sep = "."), env = probeAnno)
        chrmid <- round((chrsta + chrend)/2)
        chridx <- get(paste(chr, "index", sep = "."), env = probeAnno)
        if (checkUnique) {
            chruni <- get(paste(chr, "unique", sep = "."), env = probeAnno)
            stopifnot(length(chruni) == length(chridx))
            chridx <- chridx[chruni %in% uniqueCodes]
            chrmid <- chrmid[chruni %in% uniqueCodes]
        }
		ind <- which(chridx %in% intersect(chridx,featureNames(xSet)))        
		for (i in 1:nlevels(grouping)) {
            modSamples <- which(grouping == levels(grouping)[i])
            if (verbose) 
                cat(sampleNames(xSet)[modSamples], "... ")
            combined.dat <- as.vector(t(exprs(xSet)[chridx[ind], modSamples, 
                drop = FALSE]))
            combined.pos <- rep(chrmid[ind], each = length(modSamples))
            slidingRes <- sliding.quantile(positions = combined.pos, 
                scores = combined.dat, half.width = winHalfSize, 
                prob = quant, return.counts = TRUE)
            slidingRes <- slidingRes[seq(1, nrow(slidingRes) + 
                1 - length(modSamples), by = length(modSamples)), 
                , drop = FALSE]
            chrrm <- slidingRes[, "quantile"]
            slidingRes[, "count"] <- slidingRes[, "count"]/length(modSamples)
            areBelow <- slidingRes[, "count"] < min.probes
            if (any(areBelow)) 
                chrrm[areBelow] <- NA
            stopifnot(length(chrrm) == length(chrmid[ind]))
            newExprs[chridx[ind], i] <- chrrm
        }
    }
    sample.labels <- vector("character", nlevels(grouping))
    for (i in 1:nlevels(grouping)) sample.labels[i] <- as.character(levels(grouping)[i])
    newPD <- new("AnnotatedDataFrame", data = data.frame(label = sample.labels, 
        row.names = sample.labels), varMetadata = data.frame(varLabel = c("label"), 
        row.names = c("label")))
    newEset <- new("ExpressionSet", exprs = newExprs, phenoData = newPD)
    featureNames(newEset) <- featureNames(xSet)
    sampleNames(newEset) <- sample.labels
    return(newEset)
}

asExprSetRG <- 
function (from) 
{
    stopifnot(inherits(from, "RGList"), !is.null(from$targets), 
        !is.null(from$genes))
    from$R <- as.matrix(from$R)
    from$G <- as.matrix(from$G)
    stopifnot(nrow(from$targets) == ncol(as.matrix(from$R)))
    if (is.null(rownames(from$targets))) 
        rownames(from$targets) <- as.character(from$targets[[1]])
    colnames(from$R) <- paste(from$targets$Cy5,from$targets$Array,sep="")
    colnames(from$G) <- paste(from$targets$Cy3,from$targets$Array,sep="")
    newtargetG <- from$targets[,-(grep("Cy5",colnames(from$targets)))]
    names(newtargetG) <-  c("SlideNumber","FileName","Array","Pair","Cell.line","Trt")
    newtargetG$RG <- "G"
    newtargetR <- from$targets[,-(grep("Cy3",colnames(from$targets)))]
    names(newtargetR) <-  c("SlideNumber","FileName","Array","Pair","Cell.line","Trt")
    newtargetR$RG <- "R"
    newtargets <- rbind(newtargetG,newtargetR)
    rownames(newtargets) <- c(colnames(from$G),colnames(from$R))
    myPD <- new("AnnotatedDataFrame", data = newtargets, varMetadata = data.frame(varLabel = colnames(newtargets), 
        row.names = colnames(newtargets)))
    myEset <- new("ExpressionSet", exprs = cbind(from$G,from$R), phenoData = myPD)
    if (!is.null(from$genes$PROBE_ID)) {
        featureNames(myEset) <- make.names(from$genes$PROBE_ID, 
            unique = TRUE)
    }
    else {
        if (!is.null(from$genes$ID)) 
            featureNames(myEset) <- make.names(from$genes$ID, 
                unique = TRUE)
    }
    featureData(myEset) <- as(from$genes, "AnnotatedDataFrame")
    return(myEset)
}


