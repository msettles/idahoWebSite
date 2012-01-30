##  Created by Matt Settles
##             mattsettles@gmail.com
##             University of Idaho
## 
##  
## Adds affyQA support for:
##    wheatcdf, barley1cdf, bovinecdf
##
##  Edits: 
##    8-30-2007 by Sam Hunter: added bovine, notes, and modified approach for qc.probes 
##

##---------------------------------------------------------------------
##  Some additional notes on how/why this code works:
##  the .qcEnv contains three objects: ls(.qcEnv) ---->  "alpha"  "qc.probes" "spikes"
##  Typical Entries:
##    .qcEnv$alpha    --->  rhesuscdf 0.05  0.065
##    .qcEnv$qc.probes ---> yeast2cdf        AFFX-YER148w3_at AFFX-YER148wM_at AFFX-YER148w5_a
##    .qcEnv$spikes ---> celeganscdf      AFFX-r2-Ec-bioB-3_at AFFX-r2-Ec-bioC-3_at AFFX-r2-Ec-bioD-3_at AFFX-r2-P1-cre-3_at

## To add new items, you can follow these steps:
## 1: Add a new entry to newalpha
## 2: Add entries for spikes.  If you have an abatch, 
##    use featureNames(abatch)[grep("bio|crex", featureNames(abatch), ignore.case=T)]
##    to extract useful featurenames.  Create a new spikes object containing these, and put it in newspikes.
## 3: Add entries to qc.probes, same as above, but use:
##    featureNames(abatch)[grep("AFFX", featureNames(abatch), ignore.case=T)] to find proper names
##    for the qc.probes
##---------------------------------------------------------------------

if (is.na(get("alpha",envir=.qcEnv)["wheatcdf",][[1]] ) )
{
  ##------- 1: alpha
  newalpha = rbind(.qcEnv$alpha, wheatcdf=c(0.05, 0.065), barley1cdf=c(0.05,0.065),bovinecdf=c(0.05,0.065))
  assign("alpha",newalpha,envir=.qcEnv)
  rm(newalpha)

  ##------- 2: spikes
  wheatspikes = data.frame( "biob"="AFFX-BioB-3_at" , "bioc"="AFFX-BioC-3_at" ,
              "biod"="AFFX-BioDn-3_at" , "crex"="AFFX-CreX-3_at" ,row.names="wheatcdf")
  barleyspikes = data.frame( "biob"="AFFX-BioB-3_at" , "bioc"="AFFX-BioC-3_at" ,
              "biod"="AFFX-BioDn-3_at" , "crex"="AFFX-CreX-3_at" ,row.names="barley1cdf")
  bovinespikes = data.frame("biob"="AFFX-BioB-3_at", "bioc"="AFFX-BioC-3_at",
              "biod"="AFFX-BioDn-3_at", "crex"="AFFX-CreX-3_at", row.names="bovinecdf")

  newspikes = rbind(.qcEnv$spikes, wheatspikes, barleyspikes, bovinespikes)
  assign("spikes", newspikes, envir=.qcEnv)
  rm(wheatspikes, barleyspikes, bovinespikes, newspikes)

  ##------- 3: qc.probes
  ## This is sort of a hack, but deals with the problem of having multiple chip types with
  ## specialized QC probes.
  ## Define each set of probes, if additional probes need to be added, use cbind, then add the
  ## new probes to extraprobes, and add the extraprobes at the end of each probe list

  ## add extra probes here:
  #extraprobes = c(tubulin3=NA,tubulinM=NA,tubulin5=NA)
  extraprobes = list(tubulin3=NA,tubulinM=NA,tubulin5=NA)
  qc.probes <- cbind(.qcEnv$qc.probes, tubulin3=NA,tubulinM=NA,tubulin5=NA)

  ## Define new probsets, remember to include extraprobes
  wheatprobes = data.frame("actin3"="AFFX-Ta-actin-3_at", "actinM"="AFFX-Ta-actin-M_at",
    "actin5"="AFFX-Ta-actin-5_at", "gapdh3"="AFFX-Ta-gapdh-3_at", "gapdhM"="AFFX-Ta-gapdh-M_at",
    "gapdh5"="AFFX-Ta-gapdh-5_at", "rnapolII3"=NA, "rnapolIIM"=NA, "rnapolII5"=NA,
    "tatabp3"=NA, "tatabpM"=NA, "tatabp5"=NA, extraprobes,
    row.names="wheatcdf")

  barleyprobes = data.frame("actin3"="Contig1390_3_s_at", "actinM"= "Contig1390_M_at",
    "actin5"="Contig1390_5_at", "gapdh3"=NA, "gapdhM"=NA, "gapdh5"=NA,
    #"gapdh3"="Contig865_3_s_at" , "gapdhM"=NA , "gapdh5"="Contig865_5_s_at" ,
    "rnapolII3"=NA, "rnapolIIM"=NA, "rnapolII5"=NA,
    "tatabp3"=NA, "tatabpM"=NA, "tatabp5"=NA,
    "tubulin3"="Contig333_3_x_at", "tubulinM"="Contig333_M_at", "tubulin5"="Contig333_5_at",
    row.names="barley1cdf")

  bovineprobes = data.frame("actin3"="AFFX-Bt-actin-3_at", "actinM"="AFFX-Bt-actin-M_at",
    "actin5"="AFFX-Bt-actin-5_at", "gapdh3"="AFFX-Bt-gapd-3_at", "gapdhM"="AFFX-Bt-gapd-M_at",
    "gapdh5"="AFFX-Bt-gapd-5_at", "rnapolII3"=NA, "rnapolIIM"=NA, "rnapolII5"=NA,
    "tatabp3"=NA, "tatabpM"=NA, "tatabp5"=NA, extraprobes,
    row.names="bovinecdf")

##  TEMPLATEprobes = data.frame("actin3"=NA, "actinM"=NA,
##    "actin5"=NA, "gapdh3"=NA, "gapdhM"=NA,
##    "gapdh5"=NA, "rnapolII3"=NA, "rnapolIIM"=NA, "rnapolII5"=NA,
##    "tatabp3"=NA, "tatabpM"=NA, "tatabp5"=NA, extraprobes,
##    row.names="TEMPLATEcdf")

  qc.probes = rbind(qc.probes, wheatprobes, barleyprobes, bovineprobes)

  assign("qc.probes", qc.probes ,envir=.qcEnv)
  rm(extraprobes, qc.probes, wheatprobes, barleyprobes, bovineprobes)
  gc()
}


if (exists(".getRatios")) remove(.getRatios)

.getRatios <- function(x)
{
  namegrep3 <- function(stems,all)
    {sapply(stems,function(stem) {grep(paste(stem,"[_-]3.?_?.?_at$",sep=""),all,value=T)});}

  namegrepM <- function(stems,all)
    {sapply(stems,function(stem) {grep(paste(stem,"[_-]M.?_?.?_at$",sep=""),all,value=T)});}

  namegrep5 <- function(stems,all)
    {sapply(stems,function(stem) {grep(paste(stem,"[_-]5.?_?.?_at$",sep=""),all,value=T)});}

  vals <- x@qc.probes;
  unique.names <- colnames(vals)
  unique.names <- sub("[_-]5.?_?.?_at$","",unique.names,perl=T);
  unique.names <- sub("[_-]3.?_?.?_at$","",unique.names,perl=T);
  unique.names <- sub("[_-]M.?_?.?_at$","",unique.names,perl=T);
  unique.names <- unique(unique.names);
  p3 <- namegrep3(unique.names,colnames(vals))
  p5 <- namegrep5(unique.names,colnames(vals));
  pM <- namegrepM(unique.names,colnames(vals));
  res1 <- rbind(c(),(vals[,p3] - vals[,p5]))

  colnames(res1) <- paste(unique.names,".3'/5'",sep="")

  res2 <- rbind(c(),(vals[,p3] - vals[,pM]))
  colnames(res2) <- paste(unique.names,".3'/M",sep="")
  r <- cbind(res1,res2)
  return(r)

}

setMethod("ratios","QCStats",function(object) .getRatios(object))

