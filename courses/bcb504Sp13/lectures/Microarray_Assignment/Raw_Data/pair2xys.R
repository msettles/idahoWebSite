
### function to convert pair files to xys files given ndf file
### Matt Settles

## dir - directory of pair files
## ndf - ndf file

"pair2xys" <- function(dir,ndf){

  read.ndf <- read.table(ndf,sep="\t",as.is=T,header=T)
  pair.files <- dir(dir,pattern=".pair",full.names=TRUE)
  for (i in pair.files){
    firstLine <- readLines(i,n=1)
    pair <- read.table(i,sep="\t",as.is=T,header=T)
    order.pair <- match(paste(read.ndf$X,read.ndf$Y),paste(pair$X,pair$Y))
    pair <- pair[order.pair,]
    pair$X <- read.ndf$X
    pair$Y <- read.ndf$Y
    xysData <- data.frame(X=pair$X,Y=pair$Y,SIGNAL=pair$PM,COUNT=0)
    newfile <- sub("pair","xys",i)
    writeLines(firstLine,newfile)
    suppressWarnings(write.table(xysData,newfile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE,append=TRUE))
    cat(paste("Processed file:",i,"\n"))
  }
}

