rm(list = ls())

source("https://bioconductor.org/biocLite.R")
biocLite("piano")

pkgs <- c("yeast2.db", "limma", "affy", "plier", "affyPLM", "gtools", "plotrix")
source("http://www.bioconductor.org/biocLite.R")
biocLite(pkgs)

library(piano)

load("C:/Users/Suruchi Ahuja/Desktop/R prog/New folder/525 hw3/GeneMap.RData")
load("C:/Users/Suruchi Ahuja/Desktop/R prog/New folder/525 hw3/Array.RData")

GeneGroup <- factor(c(rep("A",3),rep("B",3), rep("A",7), rep("B",7)))

Values <- function(i) t.test(ArrayData[i,c(1:3,7:13)],ArrayData[i,c(4:6,14:20)],var.equal=T)$p.value
Pvalue <- mapply(Values, 1:1000)

adjPvalue <- p.adjust(Pvalue, method = "bonferroni")
plot(1:20,adjPvalue,col=GeneGroup)
abline(h=0.05,col="blue")

# let's look at p-values for the truly DE genes:
Pvalue[1:10]

# note: with Bonferroni correction
50*Pvalue[1:10]

# factor object for significant genes
SigFactor=factor(Pvalue <- 0.05, levels <- c(T,F))
levels(SigFactor)=c("Significant", "Non-Significant")

# creating a table
table(SigFactor,GeneGroup)

# test association of group with significance using Fisher's Exact Test
fisher.test(SigFactor,GeneGroup)



