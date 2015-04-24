require(gdata)

sheetNames<-sheetNames("C:/Users/weidewind/Documents/Asymmetry/2015/DataSetsProcessed/CleanData.xlsx")


sink("Fisher_evoless_all", append=TRUE, split=TRUE)
myRows <- integer()


for(name in sheetNames){
  myDf <- read.xls ("C:/Users/weidewind/Documents/Asymmetry/2015/DataSetsProcessed/CleanData.xlsx", sheet = name, header = TRUE, blank.lines.skip=T, stringsAsFactor=F)
  
  df<-data.frame(Strand=myDf$"Strand", Essential=myDf$"Evoless", CAI=as.numeric(myDf$"CAI"))
  df<-df[!is.na(df$"Strand"),]
  #remove(myDf)

  cleanSub <-subset(df, df$"Strand" != 2)
  cleanSub<-drop.levels(cleanSub)
  
  print ("Essentiality(Strand)")
  x <-cleanSub$"Strand"
  y <-cleanSub$"Essential"
  
  print (table(x,y))
  pval<-fisher.test(x, y,
                    or = 1, alternative = "two.sided",
                    conf.int = FALSE)
  print (pval)
  p1<-pval[["p.value"]]
  
  myRows<-c(myRows, p1)
  print (myRows)
}

final <- matrix(myRows,ncol=1,byrow=T)
colnames(final) <- c("Evoless")
rownames(final) <- sheetNames
final <- as.table(final)
final

sink()