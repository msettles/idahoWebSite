## Read in a Roche 454 quality file  which is paired to a fasta file
## function expects the two files to match each other, 
## ie the same number of reads, of the same length, in the same order
## A situation occurs when there are a large (unknown actual value) number of reads
## R is unable to hash the factor(f) below
## this function splits the reads into a more managable 500K
## this function still takes a long time to execute and should be superseded by code
## in C to read in Roche 454 Qual files
##
## usage: quals <- read.bigqual("myqualityfile.qual",DNAStringSet)
##
## parameters:
## filename: is the filename of the quality file
## reads: is a DNAStringSet object of the matching fasta file
##
## output: A BStringSet object of the corresponding qualities
##
"read.bigqual" <- 
function (filename,reads){
    nums <- scan(filename, integer(0), n = sum(width(reads)),comment.char=">")
    
    ## too big, split into 500K chunks
    y <- vector("list",length(reads))
    rcycle <- ceiling(length(reads)/500000)

    for( i in seq.int(rcycle) ){
        ireads <- seq(500000*(i-1)+1,min(length(fa),500000*i))
        f <- rep(ireads, width(reads[ireads]))
        f <- factor(f)
        storage.mode(f) <- "integer"

        tmp <- nums[seq.int(sum(width(reads[ireads])))]
        nums <- nums[-(seq.int(width(reads[ireads])))]
        y[ireads] <- .Internal(split(tmp, f))

    }

    names(y) <- names(reads)
    return(BStringSet(sapply(y, function(elt) rawToChar(as.raw(elt+33)))))
}

## Write out a Roche 454 quality data (converted to character format) back into 
## integer format
##
## usage: write.RocheQual(BStringSet)
##
## parameters:
## x: is a BStringSet of Qualities
## file: file to writ out to
##
## output: No output returned
##
## TODO: break lines into 80
"write.RocheQual" <- 
function (x, file) {
    unlink(file)
    zz <- file(file, "a")  # open an output file connection
    for (i in seq.len(length(x))){
	cat( ">",names(x)[[i]],"\n",sep="", file = zz)
	cat( as.integer(charToRaw(as.character(x[[i]])))-33, "\n",file=zz)
    }
    close(zz)
}

## Not sure how this would effect reading in the data
#    writeBString <- function(bstring)
#    {
#        if (length(bstring) == 0L)
#            return()
#        nlines <- (length(bstring) - 1L) %/% width + 1L
#        lineIdx <- seq_len(nlines)
#        start <- (lineIdx - 1L) * width + 1L
#        end <- start + width - 1L
#        if (end[nlines] > length(bstring))
#            end[nlines] <- length(bstring)
#        bigstring <- paste(
#            as.character(Views(bstring, start = start, end = end)),
#            collapse="\n"
#        )
#        cat(bigstring, "\n", file=file, sep="")
#    }

