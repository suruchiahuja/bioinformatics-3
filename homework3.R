#####################################
#
# This code is for STA525 homework3
# Group2 
# Modified: April 25,2016
#
#####################################

rm(list = ls())
load("Array.RData")
load("GeneMap.RData")
load("Valid.RData")

#----------------------------------#
# Question 1
#----------------------------------#



#----------------------------------#
# Question 2
#----------------------------------#


dim(ArrayData)
ArrayData
popDX=c(rep("A",3),rep("B",3),rep("A",7),rep("B",7))


###################
# Try Different Distance Metrics
# Euclidean
d1 = dist(t(ArrayData))
hcd1 = hclust(d1, "complete")

pdf("figures/Euclidean.pdf",width=5,height=4.5)
plot(hcd1,labels=popDX, main="Euclidean Distance Matrics")
dev.off()


clust1=cutree(hcd1,k=2)
table(cluster=factor(clust1),truth=factor(popDX))

fisher.test(x=as.integer(factor(clust1)),y=as.integer(factor(popDX)))

# Manhattan
d2 = dist(t(ArrayData), method = "manhattan")
hcd2 <- hclust(d2, "complete")
pdf("figures/Manhattan.pdf",width=5,height=4.5)
plot(hcd2,labels=popDX, main="Manhattan Distance Matrics")
dev.off()

clust2=cutree(hcd2,k=2)
table(cluster=factor(clust2),truth=factor(popDX))
fisher.test(x=as.integer(factor(clust2)),y=as.integer(factor(popDX)))


######################3
# Try different clustering methods
# Complete Method, Euclidean distance
par(mfrow=c(1,1))
d = dist(t(ArrayData))
h.complete=hclust(d)
pdf("figures/Complete.pdf",width = 5,height = 4.5)
plot(h.complete,labels=popDX,main="Complete Method")
dev.off()

clust=cutree(h.complete,k=2)
table(cluster=factor(clust),truth=factor(popDX))

fisher.test(x=factor(clust),y=factor(popDX))



# ward Method, Euclidean distance
h.ward=hclust(d,method="ward.D")
pdf("figures/Ward.pdf",width = 5,height = 4.5)
plot(h.ward,labels=popDX,main="Ward Method")
dev.off()

clust=cutree(h.ward,k=2)
table(cluster=factor(clust),truth=factor(popDX))

fisher.test(x=factor(clust),y=factor(popDX))


#######################
# Try different K values
# Euclidean distance, Complete Method

# k=2 
d1 = dist(t(ArrayData))
hcd1 = hclust(d1, "complete")
pdf("figures/K=2.pdf",width = 5,height = 4.5)
plot(hcd1,labels=popDX,main="K=2")
rect.hclust(hcd1, k=2, border="red")
dev.off()


# k=3
d1 = dist(t(ArrayData))
hcd1 = hclust(d1, "complete")
pdf("figures/K=3.pdf",width = 5,height = 4.5)
plot(hcd1,labels=popDX,main="K=2")
rect.hclust(hcd1, k=3, border="red")
dev.off()

# k=4
d1 = dist(t(ArrayData))
hcd1 = hclust(d1, "complete")
pdf("figures/K=4.pdf",width = 5,height = 4.5)
plot(hcd1,labels=popDX,main="K=2")
rect.hclust(hcd1, k=4, border="red")
dev.off()


###############################
#
# Compare the results with univariate t-test
#



#--------------------------#
#  QUESTION 3
#--------------------------#
#Cluster the features of your array data
#(a) Load the file GeneMap.RData which provides a variable 
#GeneMap that maps each feature (e.g., a known set)


head(ArrayData)
dim(ArrayData)

class(GeneMap)
#GeneMap is a vector of integers that maps each feature to a known set
#Check how many unique sets we have
length(unique(GeneMap))
#We have 50 unique sets 

#lets look at the distribution of the GeneMap
hist(GeneMap)

#(b) Perform hiearchical clustering of features.
#Do any of the clusters appear to be "enriched" sets?
mosaic <- function(data, distance.metric, linkage.method, k){
  d <- dist(data, method = distance.metric)
  hc <- hclust(d, method = linkage.method)
  clusters <- as.factor(cutree(hc, k=k))
  df <- data.frame(rownames(data), sets=GeneMap, clusters=clusters, row.names=1)
  table <- table(sets=df$sets, clusters=df$clusters)
  colors <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33')
  # mosaicplot plots contingency in array form 
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  mosaicplot(~sets + clusters,
             data=table,
             cex.axis=.9,
             main=paste0("Feature Clusters (k = ", k, ", ", distance.metric, " distance/", linkage.method, " linkage) Gene Sets"),
             col=colors[1:k])
  legendText <- c()
  for(i in 1:k){
    legendText[i] <- paste("Cluster", i, sep = " ")
  }
  legend("bottomright",
         inset=c(-0.1,0),
         legendText,
         pch=c(rep(22,k)),
         pt.bg = colors[1:k],
         col=colors[1:k])          
}

pdf("figures/hclustring.pdf")
par(mfrow=c(4,1), mar=c(2,2,2,2))
mosaic(ArrayData, "euclidean", "complete", k=2)
mosaic(ArrayData, "euclidean", "complete", k=3)
mosaic(ArrayData, "euclidean", "complete", k=4)
mosaic(ArrayData, "euclidean", "complete", k=5)
dev.off()

#A mosaic plot is a graphical display that examines the relationship
#among two or more categorical variables. 
#The mosaic plot starts as a square with a length one.
#For this example, the bars are then split vertically into the number of genesets
#The width of the vertical bars are proportional to the total number of features in the geneset
#the horizontal splits are the number the proportion of observations in each cluster (depending on k).




####   Part C
#(c) Examine the features that you identified for validation.
#Does that collection appear to be "enriched" for any GeneMap sets?

load("Valid.RData")

#we validated 10 features -- grab only non-null from dataset
validation <- ValidData[!sapply(ValidData, is.null)]
validation_df <- do.call("rbind", validation)

#now lets look at genesets
validation_set <- GeneMap[which(rownames(ArrayData) %in% rownames(validation_df))]

#these are the genesets that appear to be "enriched" for in the GeneMap sets
barplot(table(validation_set),
        ylab="Frequency",
        xlab="GeneMap Sets",
        main="Validated Features in GeneMap sets")


#-----------------------#
# Question 4
#-----------------------#
# According to our biomarker challendge result, the top differencialy expressed features are 
#   145,705,235,63,754 and 954
# We return which gene set these features belong to
GeneMap[c(145,705,235,63,754,954)]
# 1  1 18 33 27 19
# these sets are also validated according to question 3(1,3,12,18,19,27,33,50)

# So we assume that the set 1,18,33,27,19 are differentially expressed, as fearture group "G1"
#   other features are "G2"


group=c()
for (i in 1:1000) {if (GeneMap[i] %in% c(1,18,33,27,19)) group[i]="G1" else group[i]="G2"}
GeneGrp=as.factor(group)

Values = function(i) t.test(ArrayData[i,c(1:3,7:13)],ArrayData[i,c(4:6,14:20)],var.equal=T)$p.value
ObsPvalue = mapply(Values, 1:1000)


plot(1:1000,ObsPvalue,col=GeneGrp)
abline(h=0.05,col="blue")


# factor object for significant genes
SigFactor=factor(ObsPvalue<=0.05, levels <- c(T,F))
levels(SigFactor)=c("Significant", "Non-Sig")

# creating a table
table(SigFactor,GeneGrp)

# test association of group with significance using Fisher's Exact Test
fisher.test(SigFactor,GeneGrp)


# Try other p-value's cut-off
# 0.01
SigFactor1=factor(ObsPvalue<=0.01, levels <- c(T,F))
levels(SigFactor1)=c("Significant", "Non-Sig")
table(SigFactor1,GeneGrp)
fisher.test(SigFactor1,GeneGrp)

# 0.1
SigFactor10=factor(ObsPvalue<=0.1, levels <- c(T,F))
levels(SigFactor10)=c("Significant", "Non-Sig")
table(SigFactor10,GeneGrp)
fisher.test(SigFactor10,GeneGrp)




#-----------------------#
# Question 5
#-----------------------#

library(piano)

#Creating input 1: gene set collection 
genes2genesets = cbind(row.names(ArrayData), paste("GeneSet", GeneGrp))
myGeneSetCollection = loadGSC(genes2genesets)

#viewing all the genes in gene set G1
myGeneSetCollection$gsc[2]



















#--------------------------------------#
#
# Question 6
#
#--------------------------------------#

# Using the GeneMap.RData, perform a GSEA analysis of the type that we discussed in lecture.

SmplGrp=c(rep(0,3),rep(1,3),rep(0,7),rep(1,7)) # 0=Group A, 1=Group B

# We could/(should?) use t-test values instead of correlation but
# for illustrative purposes, let's use correlation

MyCor=function(i) cor(SmplGrp,ArrayData[i,])
RhoHats=mapply(MyCor,1:1000)

# aside, let's compare t-test p-vals and cor estimates
plot(ObsPvalue,RhoHats)

# let's look at corr values
plot(1:1000,RhoHats,col=GeneMap)

# let's sort the values in descending order
oDX=order(RhoHats,decreasing=T)
plot(1:1000,RhoHats[oDX],col=GeneMap[oDX],xlab="ordered gene list")

# let's make this look like a component in the GSEA plots that we have seen:
plot(1:1000,RhoHats[oDX],col=GeneMap[oDX],xlab="ordered gene list",ylab=expression(hat(rho)),pch=" ")
abline(h=0)
for(i in 1000:1) lines(c(i,i),c(0,RhoHats[oDX[i]]),col=GeneMap[oDX[i]])

ranks=1001-rank(RhoHats) # rank such that highest corr has rank 1

# put vertical lines for group A
abline(v=ranks[1:10])



















