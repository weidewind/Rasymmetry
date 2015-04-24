require(gdata)
sheetNames<-sheetNames("C:/Users/weidewind/Documents/Asymmetry/2015/DataSetsProcessed/BasicData.xlsx")
for(name in sheetNames){
print(name)
myDf <- read.xls ("C:/Users/weidewind/Documents/Asymmetry/2015/DataSetsProcessed/BasicData.xlsx", sheet = name, header = TRUE)
print ("% on leading:")
print(mean(myDf$"Strand"[myDf$"Strand" != 2]))
mySub <-subset(myDf, (myDf$"Essential" != "EXCLUDED" & myDf$"Essential" != "absent" & myDf$"Essential" != "undetermined"))
mySub <-subset(mySub, mySub$"Strand" != 2)
mySub<-drop.levels(mySub)
x <-mySub$"Essential"
y <-mySub$"Strand"
print (table(x,y))
print(fisher.test(x, y,
            or = 1, alternative = "two.sided",
            conf.int = FALSE))
}

