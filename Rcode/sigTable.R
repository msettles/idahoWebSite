##########
##  affyTable - Creates HTML and TXT Annotation Files
##
##  Written By: Matt Settles
##  Created 3/20/2007
##
##	probeids - vector of probe ids (genelist)
##  fitobj - a limma results object
##  intobj - ExprSet object
##  filename - filename to output to (do not include .html or .csv)
##  tests - a vector of multiple testing adjustments methods default=c("BH") (see p.adjust)
##  anncol - annotation columns to use default = aaf.handler()
##  annarray - annotation package to use
##  html - output html (TRUE/FALSE)
##  txt - output tab delim txt file (TRUE/FALSE)
## 
###################################################
##  2/26/08 - removed requirement for annoation packages
###########


require(limma)
require(annaffy)

probeids=sigdif
fitobj=gfit.rma
intobj=eSet.rma
filename=file.path(tablePath,"limma.rma.CRD.annoaffy")
tests=c("BH")
anncol=aaf.handler()
annarray="yeast2.db"
html=TRUE
excel=TRUE

affyTable <- function(probeids=featureNames(intobj),fitobj,intobj,filename="annoaffy",
						tests=c("BH"),anncol=aaf.handler(),annarray=annotation(intobj),
						html=TRUE,excel=TRUE){

	# Convert p-values to adjusted p-values
	adjp.BH <- sapply(as.data.frame(fitobj$p.value),function(x){p.adjust(as.matrix(x),method=tests)})
	adjp.BH <- as.data.frame(adjp.BH)
	colnames(adjp.BH) <- paste("adj.p.value",names(adjp.BH),sep=" ")
	rownames(adjp.BH) <- rownames(fitobj)
	# significance
	testtable <- aafTableSig(adjp.BH[probeids,,drop=FALSE])
	# Get Affy annotations if available
        anntable <-try(aafTableAnn(as.character(probeids) , annarray, anncol))
        if( class(anntable) == "try-error") 
	    anntable <- aafTableFrame(as.data.frame(probeids),probeids=probeids,colnames="ProbeIds")

	# M and As
	fitobj <- fitobj[probeids,]
	MAs <- aafTableFrame(cbind(fitobj$Amean,fitobj$coefficients),
				colnames=c("A",paste("M",colnames(fitobj$coefficients),sep=" ")),signed=TRUE)

	# Expression values
	exprtable <- aafTableInt(intobj[probeids,])

	table <- merge(anntable,testtable)
	table <- merge(table,MAs)
	table <- merge(table,exprtable)
	if(html) saveHTML(table,paste(filename,"html",sep="."),title=NULL)
	if(excel) saveText(table,paste(filename,"txt",sep="."))
}

aafTableAnn <-
function (probeids, chip, colnames = aaf.handler(chip = chip), 
    widget = FALSE) 
{
    colnames <- intersect(colnames, aaf.handler(chip = chip))
    if (widget) 
        colnames <- selectorWidget(aaf.handler(), colnames, ordernsel = TRUE, 
            title = "Select Annotation Data Columns")
    table <- vector("list", length(colnames))
    for (i in 1:length(colnames)) {
        table[[i]] = aaf.handler(probeids, chip, colnames[i])
    }
    names(table) = colnames
    return(new("aafTable", probeids = probeids, table = table))
}

aaf.handler <-
function (probeids, chip, name) 
{
    deps <- list(Probe = character(0), Symbol = "SYMBOL", Alias="ALIAS", GeneName = "GENENAME", Description = "DESCRIPTION", 
        Chromosome = "CHR", `Chromosome Location` = "CHRLOC", 
        GenBank = "ACCNUM", Gene = "ENTREZID", Cytoband = c("MAP", "ENTREZID"), Ensembl = "ENSEMBL", PubMed = "PMID", 
        `Gene Ontology` = "GO", Pathway = c("PATH", "ENZYME"))
    if (!missing(chip)) {
        require(paste(chip,".db",sep=""), character.only = TRUE) || stop(paste("Couldn't load data package", chip))
        use <- rep(TRUE, length(deps))
        pkgSyms <- ls(paste("package:", chip,".db", sep = ""))
        prefix <- annaffy:::annpkg_prefix(chip)
        for (i in seq(along = deps)[-1]) {
            if (any(!(paste(prefix, deps[[i]], sep = "") %in% 
                pkgSyms))) 
                use[i] <- FALSE
        }
        deps <- deps[use]
    }
    if (missing(probeids)) 
        return(names(deps))
    else switch(name, 
                    Probe = aafProbe(probeids), 
                    Symbol = aafSymbol(probeids, chip), 
                    Description = aafDescription(probeids, chip), 
                    Chromosome = aafChromosome(probeids, chip), 
                    `Chromosome Location` = aafChromLoc(probeids, chip), 
                    GenBank = aafGenBank(probeids, chip), 
                    Gene = aafLocusLink(probeids, chip), 
                    LocusLink = aafLocusLink(probeids, chip), 
                    Cytoband = aafCytoband(probeids, chip), 
                    UniGene = aafUniGene(probeids, chip), 
                    PubMed = aafPubMed(probeids, chip), 
                    `Gene Ontology` = aafGO(probeids, chip), 
                    Pathway = aafPathway(probeids, chip))
}

aafTableSig <- function(frame, colnames = names(frame), 
                          probeids = row.names(frame), pvalue = 0.05) {

    options(pvalue.sig = pvalue)
    len <- dim(frame)[1]
    
    if (sum(duplicated(colnames)))
        stop("All column names must be unique")
    
    for (col in colnames)
        if (!nchar(col))
            stop("Blank column names not allowed")
    
    table <- vector("list", dim(frame)[2])
    for (col in 1:dim(frame)[2]) {
        table[[col]] <- new("aafList", as.list(frame[,col]))
        if (pvalue)
            for (row in 1:len)
                class(table[[col]][[row]]) <- "aafSig"
    }
    names(table) <- colnames
    
    return(new("aafTable", probeids = probeids, table = table))
}

setClass("aafSig", "numeric", prototype = numeric(0))

setMethod("getTD", "aafSig", function(object) {
    sig <- getOption("pvalue.sig")
    if (object <= sig)
        class <- "aafSigValue"
    else
        class <- "aafSignedZero"
    
    return(paste("<td class=\"", class, "\">", getHTML(object), "</td>", sep = ""))       
})

setMethod("getCSS", "aafSig", function(object) {
    
    return(c("td.aafSigValue { background-color: #ff9 }"))
})

##########
## NOT USED ANYMORE
#table.rawp2adjp <-
#function (rawp, proc = c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY")) 
#{
#	m <- length(rawp)
#   n <- length(proc)
#    index <- order(rawp)
#    spval <- rawp[index]
#    adjp <- matrix(0, m, n + 1)
#    dimnames(adjp) <- list(NULL, c("rawp", proc))
#    adjp[, 1] <- spval
#    if (is.element("Bonferroni", proc)) {
#        tmp <- m * spval
#        tmp[tmp > 1] <- 1
#        adjp[, "Bonferroni"] <- tmp
#    }
#   if (is.element("Holm", proc)) {
#        tmp <- spval
#        tmp[1] <- min(m * spval[1], 1)
#        for (i in 2:m) tmp[i] <- max(tmp[i - 1], min((m - i + 
#            1) * spval[i], 1))
#        adjp[, "Holm"] <- tmp
#    }
#    if (is.element("Hochberg", proc)) {
#        tmp <- spval
#        for (i in (m - 1):1) {
#            tmp[i] <- min(tmp[i + 1], min((m - i + 1) * spval[i], 
#                1, na.rm = TRUE), na.rm = TRUE)
#            if (is.na(spval[i])) 
#                tmp[i] <- NA
#        }
#        adjp[, "Hochberg"] <- tmp
#    }
#    if (is.element("SidakSS", proc)) 
#        adjp[, "SidakSS"] <- 1 - (1 - spval)^m
#    if (is.element("SidakSD", proc)) {
#        tmp <- spval
#        tmp[1] <- 1 - (1 - spval[1])^m
#        for (i in 2:m) tmp[i] <- max(tmp[i - 1], 1 - (1 - spval[i])^(m - 
#            i + 1))
#        adjp[, "SidakSD"] <- tmp
#    }
#    if (is.element("BH", proc)) {
#        tmp <- spval
#        for (i in (m - 1):1) {
#            tmp[i] <- min(tmp[i + 1], min((m/i) * spval[i], 1, 
#                na.rm = TRUE), na.rm = TRUE)
#            if (is.na(spval[i])) 
#                tmp[i] <- NA
#        }
#        adjp[, "BH"] <- tmp
#    }
#    if (is.element("BY", proc)) {
#        tmp <- spval
#        a <- sum(1/(1:m))
#        tmp[m] <- min(a * spval[m], 1)
#        for (i in (m - 1):1) {
#            tmp[i] <- min(tmp[i + 1], min((m * a/i) * spval[i], 
#                1, na.rm = TRUE), na.rm = TRUE)
#            if (is.na(spval[i])) 
#                tmp[i] <- NA
#        }
#        adjp[, "BY"] <- tmp
#    }
#    adjp.frame <- adjp[order(index),-1]
#    adjp.frame
#}

