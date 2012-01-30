############################################################
##
##	Function loadNimblegen
##
##	Creates an nbsExprsSet object from tab delimited files
##  targets.txt, description.txt and miame.txt
##
##	Returns an nbsExpersSet object
#############################################################	
#debug

"loadNimblegen" <-
function ( targetsFile = NULL, metaFile = NULL, miameFile = NULL, dataPath = choose.dir(), notes = "" ) {
	if (!require("oligo")) 
		stop("Need to install the library: oligo")
	if ( is.null(targetsFile) ) targetsFile = "targets.txt"
	if ( is.null(metaFile) ) metaFile = "description.txt"
	if ( is.null(miameFile) ) miameFile = "miame.txt"
	# Read in Targets File
	targets <- read.table(targetsFile, header = TRUE, as.is = TRUE, sep = "\t", quote = "\"", fill = TRUE,colClasses=c("factor"))
	if ( "Name" %in% colnames(targets) )
		rownames(targets) <- targets$Name
	# metadata
	metadata <- read.table(metaFile, header=FALSE, row.names=1,  as.is = TRUE, sep = "\t", quote="\"",fill=TRUE,colClasses= c("character"))
	colnames(metadata) <- c("labelDescription")	
	# phenodata
	phenodata <- new("AnnotatedDataFrame",data=targets, varMetadata=metadata)
	# miamedata
	miamedata <- read.MIAME(filename=miameFile,widget=FALSE)
	
	if( length(grep("pair",phenodata$FileName) ) > 0 )
		NblBatch <- read.pairfiles(filenames=targets$FileName,
						pairfile.path=dataPath,
                    	pkgname=NULL,
                    	sampleNames=targets$Name,
                    	phenoData=phenodata,
                    	featureData=NULL,
                    	experimentData=miamedata,
                    	notes=notes,
                    	verbose = FALSE)
	else if( length(grep("xys",phenodata$FileName) ) > 0 )
		NblBatch <- read.xysfiles(filenames=targets$FileName,
						xysfile.path=dataPath,
                    	pkgname=NULL,
                    	sampleNames=targets$Name,
                    	phenoData=phenodata,
                    	featureData=NULL,
                    	experimentData=miamedata,
                    	notes=notes,
                    	verbose = FALSE)
	else stop("Unrecognized file type must be pair or xys files")
	
	NblBatch
}

"read.xysfiles3" <- function (..., filenames, ndf, pos, phenoData, notes, normalize=T,verbose = TRUE, sampleNames, checkType = TRUE) 
{
	require("oligo")
	readxysHeader <-
	function (filename) 
	scan(filename, nlines = 1, quiet = TRUE, what = character(0))

	getFilenames <-
	function (filenames, ...) 
	{
	    if (!missing(filenames)) {
		filenames <- c(filenames, unlist(list(...)))
	    } else {
		filenames <- unlist(list(...))
	    }
	    filenames
	}

    filenames <- getFilenames(filenames = filenames, ...)
	checkValidFilenames <-
	function (filenames) 
	{
	    stopifnot(is.character(filenames))
	    dirs <- file.info(filenames)[["isdir"]]
	    if (any(is.na(dirs))) {
		msg <- paste("These do not exist:", paste("\t", filenames[is.na(dirs)], 
		    collapse = "\n"), sep = "\n")
		stop(msg, call. = FALSE)
	    }
	    if (any(dirs)) {
		msg <- paste("These are directories:", paste("\t", filenames[dirs], 
		    collapse = "\n"), sep = "\n")
		stop(msg, call. = FALSE)
	    }
	    readable <- file.access(filenames, 4) == 0
	    if (any(!readable)) {
		msg <- paste("These are not readable:", paste("\t", filenames[!readable], 
		    collapse = "\n"), sep = "\n")
		stop(msg, call. = FALSE)
	    }
	    TRUE
	}

    checkValidFilenames(filenames)
    if (checkType)
	checkChipTypes <-
	function (filenames, verbose = TRUE, manufacturer, useAffyio) 
	{
	    if (missing(manufacturer)) 
		stop("'checkChipTypes' needs 'manufacturer'")
	    if (manufacturer == "affymetrix") {
		chips <- sapply(filenames, getCelChipType, useAffyio)
		ok <- length(unique(chips)) == 1
		if (!ok & verbose) 
		    message("All the CEL files must be of the same type.")
	    }
	    else if (manufacturer == "nimblegen") {
		designnamelist <- NULL
		for (xysfile in filenames) {
		    firstline <- readxysHeader(xysfile)
		    designname <- unlist(strsplit(firstline[grep("designname", 
		        firstline, fixed = TRUE, useBytes = TRUE)], "="))[2]
		    designnamelist <- rbind(designnamelist, designname)
		}
		ok <- length(unique(designnamelist)) == 1
		if (!ok & verbose) 
		    message("All the XYS files must be of the same type.")
	    }
	    else {
		stop("'manufacturer' ", manufacturer, " unknown")
	    }
	    ok
	}
 
        stopifnot(checkChipTypes(filenames, verbose, "nimblegen"))
    if (!missing(sampleNames)) 
        stopifnot(length(sampleNames) == length(filenames))
    firstline <- readxysHeader(filenames[1])
    designname <- unlist(strsplit(firstline[grep("designname", 
        firstline, fixed = TRUE, useBytes = TRUE)], "="))[2]
    tmp <- .Call("R_read_xys_files", filenames, verbose,PACKAGE="oligo")
    tmpExprs <- tmp[["intensities"]]
    datetime <- tmp[["date"]]
    print(datetime)
    rm(tmp)
    colnames(tmpExprs) <- Biobase::sampleNames(phenoData)
    if (normalize)
	tmpExprs <- log2(normalize.quantiles(tmpExprs))
    coordinates <- tmp$coordinates
	NDF <- read.table( ndf , sep="\t",header=TRUE,as.is=TRUE)
	POS <- read.table( pos , sep="\t",header=TRUE,as.is=TRUE)
	NDF <- NDF[match(apply(coordinates,1,paste,collapse="."),apply(NDF[,c("X","Y")],1,paste,collapse=".")),]
	NDF$CHROMOSOME <- "RANDOM"
	NDF[match(POS$PROBE_ID,NDF$PROBE_ID),"CHROMOSOME"] <- POS$CHROMOSOME
	NDF$LENGTH <- nchar(NDF$PROBE_SEQUENCE)
	NDF$COUNT <- 0
	NDF[match(POS$PROBE_ID,NDF$PROBE_ID),"COUNT"] <- POS$COUNT
	library(Biostrings)
	NDF$GC <-sapply(gregexpr("[CG]",NDF$PROBE_SEQUENCE) ,length)
	NDF <- NDF[order(NDF$CHROMOSOME,NDF$POSITION),]
	POS <- NDF

    return( new("ExpressionSet", phenoData = phenoData, featureData = new("AnnotatedDataFrame",data=POS), 
	experimentData = new("MIAME"), annotation = character(0), protocolData = phenoData[,integer(0)],
	exprs = tmpExprs) )

}

