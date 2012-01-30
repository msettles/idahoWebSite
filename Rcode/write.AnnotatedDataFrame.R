"write.AnnotatedDataFrame" <- 
function (
    x, 
    file, 
    path, 
    sep = "\t", 
    quote = FALSE, 
    varMetadata.char = "#", 
    ...
) 
{
    require("Biobase")
    if (!is(x,"AnnotatedDataFrame"))
        stop("object not of type AnnotatedDataFrame")

    if (!(is.character(varMetadata.char) && (identical(nchar(varMetadata.char), 1L)))) 
        stop("Invalid  'varMetadata.char'")

    if (!missing(path)) 
        file = file.path(path, file)

    vmd <- paste(varMetadata.char," ",rownames(varMetadata(x)),":", varMetadata(x)$labelDescription,sep="")
    writeLines(vmd,file)

    pData <- data.frame(rownames=rownames(pData(x)),pData(x))
    suppressWarnings(write.table(pData,file,sep=sep,quote=quote,col.names=TRUE,row.names=FALSE,append=TRUE,...))
}


### format cross_match output
parse_cm <- function(filename){
    lines <- readLines(filename)
    align <- grep("^ALIGNMENT",lines,value=T)
    align <- sub(" C ", " ",align)
#    cm_out <- read.table(filename)
    cm_out <- matrix(unlist(strsplit(align,split=" +")),ncol=13,byrow=T)
    cm_out <- data.frame(cm_out,stringsAsFactors=FALSE)
    cm_out$FC <- "F"
    cm_out$FC[grep("(",cm_out$X11,fixed=TRUE)] <- "C"
    cm_out[cm_out$FC == "C",c("X11","X12","X13")] <- cm_out[cm_out$FC == "C",c("X13","X12","X11")]
    colnames(cm_out) <- c("Alignment","score","perc_sub","perc_del","perc_ins","read_id",
        "read_start","read_end","read_remain","adapt","adapt_start","adapt_end","adapt_remain","FC")   
    cm_out$read_start <- as.numeric(cm_out$read_start)
    cm_out$read_end <- as.numeric(cm_out$read_end)
    cm_out$read_remain <- as.numeric(gsub("[()]","",cm_out$read_remain))
    cm_out$adapt_remain <- as.numeric(gsub("[()]","",cm_out$adapt_remain))
    cm_out$adapt_start <- as.numeric(cm_out$adapt_start)
    cm_out$adapt_end <- as.numeric(cm_out$adapt_end)
    cm_out$score <- as.numeric(cm_out$score)
    cm_out$perc_sub <- as.numeric(cm_out$perc_sub)
    cm_out$perc_del <- as.numeric(cm_out$perc_del)
    cm_out$perc_ins <- as.numeric(cm_out$perc_ins)

    cm_out$perc_sub <- round(cm_out$perc_sub *(cm_out$read_end - cm_out$read_start +1),digits=0)/100
    cm_out$perc_del <- round(cm_out$perc_del *(cm_out$read_end - cm_out$read_start +1),digits=0)/100
    cm_out$perc_ins <- round(cm_out$perc_ins *(cm_out$read_end - cm_out$read_start +1),digits=0)/100

    cm_out$read_len <- cm_out$read_end+cm_out$read_remain

    cm_out$err <- apply(cm_out[,c("perc_sub","perc_del","perc_ins","adapt_remain","adapt_start")],1,sum) -1
    cm_out
}


