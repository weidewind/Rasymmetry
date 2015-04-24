##  Anfisa, February 2015
##  Example script showing how to wield a GEO entry.



## Installing Bioconductor and biocLite package, connecting to libraries.
## Use fresh = T on a system without Bioconductor.

install <- function (fresh = F){
    if (fresh) {
        source("http://www.bioconductor.org/biocLite.R")
        biocLite("GEOquery")
    }
    library(Biobase)
    library(GEOquery)
}

##  Reading the NCBI's GEO microarray SOFT files (GDS), which can be zipped.
##  Showing gds features. 

showGDS <- function(file){
    gds <- getGEO(filename = file)
  
    ##   some GDS features  
    Meta(gds)$sample_organism
    Meta(gds)$sample_count
    colnames(Table(gds))
    gds
}


##  Use to delete "_at" at the end of rownames (to get locus names) in expression set.
chopRowNames <- function(eset, length = 3){
    v <- vapply(rownames(exprs(eset)),
                function (e) substring(e, 1,nchar(e)-length), 
                character(1))
    names(v) <- NULL
    rownames(exprs(eset)) <- v
    eset
}


##  Merges (by locus column) expression data (gds) for one sample (name) with our basic data.
##  Adds boolean column, where True means that gene expression is high (top 10%)

attachGSM <- function (basic, gsm, gds, name, locus){
    if (is.na(gsm)){
        gnew<-data.frame(locus, "Expression" = as.numeric(Table(gds)[[name]]))
    }
    else {
        
        if (!is.na(Table(gsm)$ABS_CALL)){
            gnew<-data.frame(locus, "Expression" = as.numeric(Table(gsm)$VALUE), 
                                    "Call" = as.character(Table(gsm)$ABS_CALL))
            gnew<-gnew[gnew$"Call" != "A" & gnew$"Call" != "M",] ## Delete data for Absent and Marginal calls
        }
        else if (!is.na(Table(gsm)$Detection)){
            gnew<-data.frame(locus, "Expression" = as.numeric(Table(gsm)$VALUE), 
                                    "Call" = as.character(Table(gsm)$Detection))
            gnew<-gnew[gnew$"Call" != "A" & gnew$"Call" != "M",] ## Delete data for Absent and Marginal calls
        }
        else {
            gnew<-data.frame(locus, "Expression" = as.numeric(Table(gsm)$VALUE))
        }
    }
    merged <- merge (basic, gnew, by.x = "Locus", by.y = "locus" )
    print(merged)
    sorted<<-merged[order(-merged["Expression"]),]
    exprThreshold<-sorted[round(nrow(sorted)/10),"Expression"]
    exprBoolean<-sorted[["Expression"]]>=exprThreshold
    sorted<-data.frame(sorted, exprBoolean)
}

processGSM <- function(gsm, basic){
    locus <- vapply(as.character(Table(gsm)$ID_REF),
                    function (e) { t <- substring(e, 1,nchar(e)-3)
                                   gsub("T", "T_", t, perl=TRUE) },   
                    character(1))
    data <- attachGSM(basic=basic, gsm=gsm, locus=locus)
    x <-data$"Strand"
    y <-data$"exprBoolean"
    acc <- Meta(gsm)$geo_accession
    desc <- Meta(gsm)$description
    sn <- Meta(gsm)$source_name_ch1
    print(acc)
    print(desc)
    print (sn)
    print (table(x,y))
    pval<-fisher.test(x, y,
                      or = 1, alternative = "two.sided",
                      conf.int = FALSE)
    print (pval)
    finalList <- list(as.character(acc), as.character(desc), as.character(sn), pval[["p.value"]])
    return (finalList)
}


# gse <- getGEO("GSE781",GSEMatrix=FALSE)
##  Takes path to GSE, path to basic data (in csv format), output filename
##  Prints in output Fisher p-values for Expression(Strand) (top 10%, defined in attachGSM) for all GSMs.
##  Unlike #processGDE, uses lapply instead of for-loop
##  Locus trimming and Fisher statistics are done in #processGSM

processGSE <- function(GSEfile, CSVfile, OUTfile){
    gse <- getGEO(filename = GSEfile)
    gsmlist <- GSMList(gse)

    
    basic <- read.csv (CSVfile,  header = TRUE, blank.lines.skip=T, stringsAsFactor=F)
    basic<-try(data.frame(Locus=basic$"Locus", 
                          Strand=basic$"Strand", 
                          Evoless=basic$"Evoless",
                          Essential=basic$"Essential", 
                          CAI=as.numeric(basic$"CAI")))
    if (inherits(basic, "try-error")){}
    basic<-basic[basic$"Strand" != 2,]
    
    sink(OUTfile, append=TRUE, split=TRUE)
    final <<- lapply(gsmlist, function(e) processGSM(e, basic=basic))
    names(final) <- NULL
    ids <- vapply(final, function(x) as.character(x[1]), character(1))
    desc <- vapply(final, function(x) as.character(x[2]), character(1))
    sn <- vapply(final, function(x) as.character(x[3]), character(1))
    pval <- vapply(final, function(x) as.numeric(x[4]), 0)
    final <- data.frame(ids, desc, sn, pval)
    print(final)
    sink()
  
}


processGSE('C:/Users/weidewind/Documents/Asymmetry/2015/RawDataSets/GEO/Caul/GSE52375_family.soft.gz',
           "C:/Users/weidewind/Documents/Asymmetry/2015/DataSetsProcessed/CSV/Caulo.csv",
           "Caulo_GSE52375_family")

##  Takes path to GDS, path to basic data (in csv format), output filename
##  Prints in output Fisher p-values for Expression(Strand) (top 10%, defined in attachGSM) for all GSMs.
processGDS <- function(GDSfile, CSVfile, OUTfile){
    g <-showGDS(GDSfile)
    
    ##  Convert GDS into expression set object (using base 2 logarithms)
    eset <- GDS2eSet(g, do.log2=TRUE)  
    
    myDf <- read.csv (CSVfile),  header = TRUE, blank.lines.skip=T, stringsAsFactor=F)
    df<-try(data.frame(Locus=myDf$"Locus", Strand=myDf$"Strand", Evoless=myDf$"Evoless", Essential=myDf$"Essential", CAI=as.numeric(myDf$"CAI")))
    if (inherits(df, "try-error")){
    }
    df<-df[df$"Strand" != 2,]
    sample_count <-  Meta(g)$sample_count
    locus <- vapply(as.character(Table(g)$ID_REF),
                    function (e) { t <- substring(e, 1,nchar(e)-3)
                                   gsub("T", "T_", t, perl=TRUE) },   
                    character(1))
    
    
    sink(OUTfile, append=TRUE, split=TRUE)
    final <- data.frame(txt=rep("u", sample_count), txt=rep("u", sample_count), num=rep(NA, sample_count), # as many cols as you need
                        stringsAsFactors=FALSE)          # you don't know levels yet
    i <- 1
    for (name in colnames(Table(g))){
      if (substring(name,1,3) == "GSM"){
          data <- attachGSM(df, g, name, locus)
          x <-data$"Strand"
          y <-data$"exprBoolean"
          print (name)
          print (pData(eset[,name])$growth)
          print (table(x,y))
          pval<-fisher.test(x, y,
                            or = 1, alternative = "two.sided",
                            conf.int = FALSE)
          print (pval)
          
          final[i,]<-list(as.character(pData(eset[,name])$growth), as.character(pData(eset[,name])$time),pval)
          i <- i+1
      }
    }
    
    colnames(final) <- c("growth.protocol", "time", "pval")
    final
    sink()
}


## Example usage of eset:

##> eset["1320_at","GSM14504"]
##  Expression Set (exprSet) with 
##  1 genes
##  1 samples
##  phenoData object with 4 variables and 1 cases
##  varLabels
##  : sample
##  : infection
##  : genotype/variation
##  : description

##> exprs(eset["1320_at","GSM14504"])
##  GSM14504
##  1320_at  6.70044
