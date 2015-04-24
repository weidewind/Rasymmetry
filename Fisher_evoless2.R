require(gdata)
require(vcd)

files <<- NULL
files <-list.files("C:/Users/weidewind/Documents/Asymmetry/2015/DataSetsProcessed/CSV/")

sink("Evoless_all_corr2", append=TRUE, split=TRUE)
myRows <- integer()


for(f in files){
  
  name <- substring(f, 1,nchar(f)-4)
  print (name)
  myDf <- read.csv (paste("C:/Users/weidewind/Documents/Asymmetry/2015/DataSetsProcessed/CSV/", f, sep=""),  header = TRUE, blank.lines.skip=T, stringsAsFactor=F)
  
  df<-try(data.frame(Strand=myDf$"Strand", Evoless=myDf$"Evoless", Essential=myDf$"Essential", CAI=as.numeric(myDf$"CAI")))
  if (inherits(df, "try-error")){
    myRows<-c(myRows, NA, NA, NA, NA, NA)
    next
  }
  
  df<-df[!is.na(df$"Strand"),]
  #remove(myDf)
  
  cleanSub <-subset(df, df$"Strand" != 2)
  cleanSub <-subset(cleanSub, ((cleanSub$"Essential" == "essential") | (cleanSub$"Essential" == "nonessential")))
  cleanSub <-drop.levels(cleanSub)
  
  print ("Evoless(Strand)")
  x <-cleanSub$"Strand"
  y <-cleanSub$"Evoless"
  
  print (table(x,y))
  pval<-fisher.test(x, y,
                    or = 1, alternative = "two.sided",
                    conf.int = FALSE)
  print (pval)
  p0<-pval[["p.value"]]
  
  
  print ("Evoless(Essential)")
  x <-cleanSub$"Essential"
  y <-cleanSub$"Evoless"
  
  print (table(x,y))
  pval<-fisher.test(x, y,
                    or = 1, alternative = "two.sided",
                    conf.int = FALSE)
  print (pval)
  p1<-pval[["p.value"]]
  
  
  
  
  a<-assocstats(table(x,y))
  
  myRows<-c(myRows, p0, p1, a$phi, a$cont, a$cramer)
  print (myRows)
  
}

final <- matrix(myRows,ncol=5,byrow=T)
colnames(final) <- c("p-val_evoless(strand)", "p-val", "phi", "cont", "cramer")
rownames(final) <- files
final <- as.table(final)
final

sink()
