"QualityAssurance" <-
function () 
{
    Try(ArraysLoaded <- get("ArraysLoaded", envir = affylmGUIenvironment))
    Try(if (ArraysLoaded == FALSE) {
        Try(tkmessageBox(title = "Quality Assurance", message = "Error: No arrays have been loaded.", 
            icon = "error", default = "ok"))
        return()
    })
    Require("affy")
    Require("affyQCReport")
    
    if (as.numeric(version$minor) < 6) # load info for wheat barley and bovine
	    source("http://bioinfo-mite.crb.wsu.edu/Rcode/modifiedQA.R")
    
    Try(RawAffyData <- get("RawAffyData", envir = affylmGUIenvironment))
    Try(qaInfo <- affyQAReport( RawAffyData, overwrite = FALSE))
    Try(outdir <- file.path(getwd(), "affyQA"))
    Try(openQAReport(outdir))
}
