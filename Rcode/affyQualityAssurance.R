###############################################################
## Generate Quality Reports from Affymetrix Data 
###############################################################
## Matt Settles
## Feb 25, 2010
##
##    Includes, the affyQAReport output
##       weight and residual pseudoimages
##       and RData object of all data
##
###############################################################

#############################################################
### BEGIN Quality Assessment ###
#############################################################

################
## PRODUCE QC REPORTS
################
"ProduceQCreports" <- 
function(RawData, name= deparse(substitute(RawData)),outdir = paste("QA",deparse(substitute(RawData)),sep="_"),overwrite=TRUE,verbose=T) {
## debug
#RawData <- E_TABM_110
#name= deparse(substitute(E_TABM_110))
#outdir = paste("QA",deparse(substitute(E_TABM_110)),sep="_")
#overwrite=TRUE
## end debug
    require(affyQCReport)
    codePath <- "http://bioinfo-mite.ibest.uidaho.edu/Rcode"

    ## make sure simpleaffy qcdef for array is present
    ## extdata directly must be writable
    if (verbose) print("checking for qcdef information")
    qcname <- paste(cleancdfname(cdfName(RawData)),".qcdef",sep="")
    if (!file.exists(file.path(system.file("extdata",package="simpleaffy"),qcname))){
    ## try to get it from my site, simpleaffy/extdata must be writable
        if (file.access(system.file("extdata",package="simpleaffy"),mode=2) < 0)
        	stop(paste("Error:", system.file("extdata",package="simpleaffy"), "\nmust be writable"))
        
        qcfile <- file.path(codePath,"simpleaffy-qcdef",qcname)
        if (download.file(qcfile, file.path(system.file("extdata",package="simpleaffy"),qcname),method="auto") != 0){
            stop(paste("Error: qcdef file \"", qcname, "\" is not available\n",
                "contact Matt Settles <msettles@uidaho.edu> to request one be created"))
        }
    }
    if (verbose) print("creating affyQAReport")
    ## there is a undesired side effect to affyQAReport when given a folder for outdir
	source("http://bioinfo-mite.ibest.uidaho.edu/Rcode/affyQAReport.R")
    try(qaInfo <- affyQAReport2( RawData, outdir = outdir, overwrite = overwrite,repName="affyQAReport"))
	if (inherits(qaInfo, "try-error")) stop("affyQAReport failed")
    if (verbose) print("gathering RNA Degradation data")	
    ### AffyRNADeg Table ###
    rnaDeg = AffyRNAdeg(RawData)
    qaInfo$rnaDeg = rnaDeg
    if (exists("qaInfo"))
	    save(qaInfo, file=file.path(outdir, "QAinfo.RData") )
	if (verbose) print("copying affyQAReport.pdf to main QA folder")
	## copy the affyQAReport to main folder
    file.copy(from=file.path(outdir,"affyQAReport","affyQAReport.pdf"), overwrite=TRUE,to=file.path(outdir,"affyQAReport.pdf"))     

    ################
    ## AffyPLM
    ################
	if (verbose) print("getting affyPLM data")    

	file.remove(dir(path=outdir,pattern="png"))
    require(affyPLM)
    if ( exists("qaInfo") ){
        pset <- qaInfo$affyPLM
    } else if ( !exists("qaInfo") & file.exists(outdir, "QAinfo.RData")){ ## can use fitPLM object from qaInfo 
	    load(file=file.path(outdir, "QAinfo.RData") )
	    pset <- qaInfo$affyPLM
    } else {
        pset <- fitPLM(RawData, background = FALSE, normalize = FALSE)
    }

	if (verbose) print("creating weight pseudoimages")
    # weight images
    for( i in 0:floor((length(RawData)-1)/8) ){
	    png(filename=file.path(outdir,paste("weight",i,".png",sep="")),width=1024,height=2048, pointsize=10)
	    par(mfrow=c(4,2))
	    for ( j in 1:8 ){
		    if ((i*8+j) <= length(RawData)) image(pset, which = i*8 + j, type="weights")
	    }
	    dev.off()
    }
    if (verbose) print("creating residual pseudoimages")
    # residuals
    for( i in 0:floor((length(RawData)-1)/8) ){
	    png(filename=file.path(outdir,paste("resid",i,".png",sep="")),width=1024,height=2048, pointsize=10)
	    par(mfrow=c(4,2))
	    for ( j in 1:8 ){
		    if ((i*8+j) <= length(RawData) ) image(pset, which = i*8 + j, type="resids")
	    }
	    dev.off()
    }
	if (verbose) print("Finished Generating QA info")    

}
#############################################################
### END Quality Assessment ###
#############################################################

