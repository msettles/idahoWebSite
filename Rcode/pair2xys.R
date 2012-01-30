pairfolder <- "PairData"
xysfolder <- "XYSData"
pair2xys <- function(pairfolder,xysfolder, pattern="_pair.txt"){
	files <- dir(pairfolder,pattern=pattern,full.names=F)
        for (file in files){
		firstline <- readLines(file.path(pairfolder,file),n=1)
		dd <- read.table(file.path(pairfolder,file),sep="\t",comment.char="#",header=T,stringsAsFactors=T)
		xys <- dd[,c("X","Y","PM")]
                colnames(xys) <- c("X","Y","SIGNAL")
		xys$COUNT <- NA
		xys <- xys[order(xys$Y,xys$X),]
		file <- sub(pattern,".xys",file)
		writeLines(firstline,file.path(xysfolder,file))
		suppressWarnings(TRUE)
		write.table(xys,file.path(xysfolder,file),sep="\t",row.names=F,col.names=T,quote=F,append=T)
		suppressWarnings(FALSE)
	}
}


