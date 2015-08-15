############################################################
##
##	Function loadAffy
##
##	Creates an AffyBatch object from tab delimited files
##
##	Returns an AffyBatch object
#############################################################	
"loadAffyOLD" <-
function ( targetsFile = "targets.txt", metaFile = "description.txt", miameFile = "miame.txt", celPath = ".", notes = "", cdfname=NULL) {
	if (!require("affy")) 
		stop("Need to install the library: affy")
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
	# Read in CEL files
	ReadAffy(filenames=targets$FileName,
              celfile.path=celPath,
              sampleNames=rownames(targets),
              phenoData=phenodata,
              description=miamedata,
              notes=notes,cdfname)
}

"loadAffy" <-
function ( targetsFile = "targets.txt", miameFile = "miame.txt", celPath = ".", notes = "", cdfname=NULL ) {
	if (!require("affy")) 
		stop("Need to install the library: affy")
	# Read in Targets File
	phenodata <- read.AnnotatedDataFrame(targetsFile, header = TRUE, sep = "\t")
	# miamedata
	miamedata <- read.MIAME(filename=miameFile,widget=FALSE)
	# Read in CEL files
	ReadAffy(filenames=phenodata$FileName,
              celfile.path=celPath,
              phenoData=phenodata,
              description=miamedata,
              notes=notes,cdfname=cdfname)
}

"loadGPR" <-
function ( targetsFile = "targets.txt", metaFile = "description.txt", miameFile = "miame.txt", spotfile = "spottypes.txt", gprPath = ".", notes = "",pkgName ) {
	if (!require("limma")) 
		stop("Need to install the library: limma")
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
	# Read in CEL files
    targets$Label <- targets$Name
    mywtfun <- function(x) as.numeric(x$Flags > -50.5)
    rawdata <- read.maimages(files = targets,
                            source = "genepix",
                            path = gprPath,
                            ext = NULL,
                            names = NULL,
                            other.columns = c("Flags"),
                            annotation = NULL,
                            wt.fun = mywtfun,
                            verbose = TRUE,
                            sep = "\t",
                            quote = NULL,
                            columns = list(R = "F633 Mean", G = "F543 Mean", 
                                           Rb = "B633 Median", Gb = "B543 Median"))
    spottypes <- readSpotTypes(spotfile)
    rawdata$genes$Status <- controlStatus(spottypes, rawdata)
    rawdata$miameData <- miamedata
    rawdata$notes <- notes
    rawdata$pData <- phenodata
    rawdata
           
}

"loadAffyOligo" <-
function ( targetsFile = "targets.txt", metaFile = "description.txt", miameFile = "miame.txt", celPath = ".", notes = "", pkgname) {
	if (!require("oligo")) 
		stop("Need to install the library: oligo")
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
	# Read in CEL files
	read.celfiles(filenames=file.path(celPath,targets$FileName),
			phenoData=phenodata,
			experimentData=miamedata, 
			notes=notes,
			sampleNames=rownames(targets),
            pkgname=pkgname,
			verbose = TRUE)
	
}
