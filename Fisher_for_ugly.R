
require(gdata)

sink("Caje_fisher", append=TRUE, split=TRUE)
myRows <- integer()
essentials<-c("Essential_FBA", "Essential_transp_FBA_article", "Essential_transp_Stahl")
myDf <- read.xls ("C:/Users/weidewind/Documents/Asymmetry/2015/DataSetsProcessed/BasicData.xlsx", sheet = "Caje", header = TRUE, blank.lines.skip=T, stringsAsFactor=F)

for(e in essentials){

  df<-data.frame(Strand=myDf$"Strand", Essential=myDf$e, CAI=as.numeric(myDf$"CAI"))
  df<-df[!is.na(df$"Strand"),]
  #remove(myDf)
  
  print ("% on leading:")
  percent <- mean(df$"Strand"[df$"Strand" != 2 & !is.na(df$"Strand")])
  print(percent)
  cleanSub <-subset(df, df$"Essential" == "essential" | df$"Essential" == "nonessential"))
cleanSub <-subset(cleanSub, cleanSub$"Strand" != 2)
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




sortedDf<-df[order(-df$"CAI"),]
cai_threshold<-sortedDf[round(nrow(sortedDf)/10),3]
expr_f<-sortedDf$"CAI">=cai_threshold
sortedDf<-data.frame(sortedDf, expr_f)
sortedDf <-subset(sortedDf, sortedDf$"Strand" != 2)
sortedDf<-drop.levels(sortedDf)
cleanSub <-subset(sortedDf, (sortedDf$"Essential" != "EXCLUDED" & sortedDf$"Essential" != "absent" & sortedDf$"Essential" != "undetermined"  & sortedDf$"Essential" != "potential_essential"))
cleanSub<-drop.levels(cleanSub)


print("Essentiality(Strand) in HExp subgroup")
HExpsub<-subset(cleanSub, cleanSub$"expr_f" == TRUE)
x <-HExpsub$"Strand"
y <-HExpsub$"Essential"
print (table(x,y))
pval<-fisher.test(x, y,
                  or = 1, alternative = "two.sided",
                  conf.int = FALSE)
print (pval)
p7<-pval[["p.value"]]

print("Essentiality(Strand) in NHExp subgroup")
NHExpsub<-subset(cleanSub, cleanSub$"expr_f" == FALSE)
x <-NHExpsub$"Strand"
y <-NHExpsub$"Essential"
print (table(x,y))
pval<-fisher.test(x, y,
                  or = 1, alternative = "two.sided",
                  conf.int = FALSE)
print (pval)
p8<-pval[["p.value"]]


print("CAI(Essentiality)")
x <-cleanSub$"Essential"
y <-cleanSub$"expr_f"
print (table(x,y))
pval<-fisher.test(x, y,
                  or = 1, alternative = "two.sided",
                  conf.int = FALSE)

print (pval)
p2<-pval[["p.value"]]

print("CAI(Strand)")
x <-sortedDf$"Strand"
y <-sortedDf$"expr_f"
print (table(x,y))
pval<-fisher.test(x, y,
                  or = 1, alternative = "two.sided",
                  conf.int = FALSE)
print (pval)
p3<-pval[["p.value"]]


print("CAI(Strand) in essential subgroup")
Essub<-subset(cleanSub, cleanSub$"Essential" == "essential")
x <-Essub$"Strand"
y <-Essub$"expr_f"
print (table(x,y))
pval<-fisher.test(x, y,
                  or = 1, alternative = "two.sided",
                  conf.int = FALSE)
print (pval)
p4<-pval[["p.value"]]


print("CAI(Strand) in nonessential subgroup")
NEssub<-subset(cleanSub, cleanSub$"Essential" == "nonessential")
x <-NEssub$"Strand"
y <-NEssub$"expr_f"
print (table(x,y))
pval<-fisher.test(x, y,
                  or = 1, alternative = "two.sided",
                  conf.int = FALSE)
print (pval)
p5<-pval[["p.value"]]

print("CAI(Strand) in undetermined subgroup")
NDEssub<-subset(sortedDf, (sortedDf$"Essential" != "essential" & sortedDf$"Essential" != "nonessential" ))
x <-NDEssub$"Strand"
y <-NDEssub$"expr_f"
print (table(x,y))
if(length(levels(factor(x))) == 2 & length(levels(factor(y))) == 2){
  pval<-fisher.test(x, y,
                    or = 1, alternative = "two.sided",
                    conf.int = FALSE)
  print (pval)
  p6<-pval[["p.value"]]
} else p6 <- NA


myRows<-c(myRows, percent, p1,p7,p8, p2,p3,p4,p5, p6)
print (myRows)
}

final <- matrix(myRows,ncol=9,byrow=T)
colnames(final) <- c("%_on_leading", "Essentiality(Strand)", "Essentiality(Strand)_in_HExp_group", "Essentiality(Strand)_in_NHExp_group", "CAI(Essentiality)", "CAI(Strand)", "CAI(Strand)_in_essential_subgroup", "CAI(Strand)_in_nonessential_subgroup", "CAI(Strand)_in_undetermined_subgroup")
rownames(final) <- sheetNames[10]
final <- as.table(final)
final

sink()