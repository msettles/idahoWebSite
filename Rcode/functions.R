### Functions


checksystem <- function(){
 print("Checking for required system applications")
 stopifnot(all(
 	system("which sfffile") ==0,
 	system("which sffinfo") ==0,
 	system("which sff_extract") ==0,
 	system("which parseSFF_XML.pl") ==0,
 	system("which cp") ==0,
 	system("which sed") ==0,
 	system("which cross_match") ==0,
 	system("which grep") ==0,
 	system("which mv") ==0,
 	system("which wc") ==0,
 	system("which lucy") ==0,
 	system("which lucy2tab.pl") ==0))
}

## trimmed standard deviation
sd.trim <- function(x, trim=0, na.rm=FALSE, ...)
{
  if(!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(NA_real_)
  }
  if(na.rm) x <- x[!is.na(x)]
  if(!is.numeric(trim) || length(trim) != 1)
    stop("'trim' must be numeric of length one")
  n <- length(x)
  if(trim > 0 && n > 0) {
     if(is.complex(x)) stop("trimmed sd are not defined for complex data")
     if(trim >= 0.5) return(0)
     lo <- floor(n * trim) + 1
     hi <- n + 1 - lo
     x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
  }
  sd(x)
}

### format cross_match output
parse_cm <- function(filename){
    lines <- readLines(filename)
    align <- grep("^ALIGNMENT",lines,value=T)
    align <- sub(" C ", " C_",align)
#    cm_out <- read.table(filename)
    cm_out <- matrix(unlist(strsplit(align,split=" +")),ncol=13,byrow=T)
    cm_out <- data.frame(cm_out,stringsAsFactors=FALSE)
    cm_out$FC <- "F"
    cm_out$FC[grep("^C_",cm_out$X10)] <- "C"
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

    cm_out$err <- apply(cm_out[,c("perc_sub","perc_del","perc_ins","adapt_remain","adapt_start")],1,sum)
    cm_out
}


### relies on SampleNote column in the annotation
abundance.table <- function(rdp.otu, groups, anno,samples_note = "ID",rdpThres = 0.5,outfile="abundance.out"){
    numberoflevels <- (ncol(rdp.otu)-5)/3
    rdp.otu$assign <- as.character(rdp.otu$V6)
    rdp.otu$level <- as.character(rdp.otu$V7)

    rdp.otu$assign[(rdp.otu$V8 < 0.5)] <- "Unclassified"
    rdp.otu$level[(rdp.otu$V8 < 0.5)] <- "Unknown"

    for (i in seq.int(1,numberoflevels-1)){
      column = 6+3*i
      prows <- rdp.otu[, (column + 2)] >= rdpThres & !is.na(rdp.otu[, (column + 2)])
      rdp.otu$assign[prows] <- as.character(rdp.otu[,column][prows])
      rdp.otu$level[prows] <- as.character(rdp.otu[,column+1][prows])
    }
 
    rdp.otu$id <- groups$ID[match(rdp.otu$V1,groups$Acc)]
    
    rdp.otu$name <- anno[match(rdp.otu$id,anno$ID),samples_note] ### Annotation

    abundanceTable <- table(rdp.otu$assign,rdp.otu$name)
    abundanceTable <- cbind(abundanceTable,groupTotal=rowSums(abundanceTable))
    abundanceTable <- rbind(abundanceTable,ReadTotals=colSums(abundanceTable))
    abundanceTable <- cbind(Level=rdp.otu$level[match(rownames(abundanceTable),rdp.otu$assign)],abundanceTable)
    abundanceTable <- data.frame(Taxon_Name = rownames(abundanceTable),abundanceTable)
    write.table(abundanceTable,file=outfile,sep="\t",quote=F,col.names=T,row.names=F)


}


