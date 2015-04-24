GSEfile <-"C:/Users/weidewind/Documents/Asymmetry/2015/RawDataSets/GEO/Caulo/GSE26977_family.soft.gz"
gse <- getGEO(filename = GSEfile)
gsmlist <- GSMList(gse)
Table(gsmlist[[1]])
CSVfile <-"C:/Users/weidewind/Documents/Asymmetry/2015/DataSetsProcessed/CSV/Caulo.csv"
basic <- read.csv (CSVfile,  header = TRUE, blank.lines.skip=T, stringsAsFactor=F)
basic<-try(data.frame(Locus=basic$"Locus", 
                      Strand=basic$"Strand", 
                      Evoless=basic$"Evoless",
                      Essential=basic$"Essential", 
                      CAI=as.numeric(basic$"CAI")))
basic<-basic[basic$"Strand" != 2,]
gsm <- gsmlist[[1]]
gpls <- GPLList(gse)
idsToLocus <-hash(keys=Table(gpls[[1]])$ID, values=Table(gpls[[1]])$ORF)

gpls[[gsm@header$platform_id]]
gsm@header$channel_count