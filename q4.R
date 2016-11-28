rm(list = ls())

source("https://bioconductor.org/biocLite.R")
biocLite("piano")

pkgs <- c("yeast2.db", "limma", "affy", "plier", "affyPLM", "gtools", "plotrix")
source("http://www.bioconductor.org/biocLite.R")
biocLite(pkgs)

library(piano)

load("C:/Users/Suruchi Ahuja/Desktop/R prog/New folder/525 hw3/GeneMap.RData")


GeneGroup <- factor(c(rep("A",10),rep("B",15),rep("C",25)))

SimDatMVN2=function(p.alt1=10, # the number of differentially expressed features in set 1
                    p.alt2=10, # the number of differentially expressed features in set 2       
                    p.null=80, # the number of non-differentially expressed features
                    n=20, # the number of samples in each of the treatment and control groups
                    rho.alt1=0.7, # correlation of the alt1 hyp variables
                    rho.alt2=0.7, # correlation of the alt2 hyp variables
                    rho.null=0.0, # correlation of the null hyp variables
                    delta1=2,# the mean of the p.alt1 features in the "treatment" group 
                    delta2=1 # the mean of the p.alt1 features in the "treatment" group 
                    
)
  {
  
  
  p.alt=p.alt1+p.alt2
  p=p.alt+p.null
  Sigma=array(rep(0,p^2),dim=c(p,p))
  Sigma[1:p.alt1,1:p.alt1]=rho.alt1
  Sigma[(p.alt1+(1:p.alt2)),(p.alt1+(1:p.alt2))]=rho.alt2
  Sigma[(p.alt+(1:p.null)),(p.alt+(1:p.null))]=rho.null
  diag(Sigma)=1
  Xc=mvrnorm(n,mu=rep(0,p),Sigma=Sigma)
  Xt=mvrnorm(n,mu=c(rep(delta1,p.alt1),rep(delta2,p.alt2),rep(0,p.null)),Sigma=Sigma)
  x=t(rbind(Xc,Xt))
  colnames(x)=c(paste("cntrl",1:n,sep=""),paste("trt",1:n,sep=""))
  rownames(x)=c(paste("GeneGrpA",1:(p.alt1)),paste("GeneGrpB",1:(p.alt2)),paste("GeneGrpC",1:p.null,sep=""))
  
  return(x)
}

set.seed(123456)
X=SimDatMVN2(p.alt1=10,p.alt2=15,p.null=25,n=10,rho.alt1=0.5,rho.alt2=0.5,rho.null=0.25,delta1=1,delta2=0)

# visualize the data matrix
require(gplots)
heatmap.2(X,Rowv=NA,Colv=NA,trace="none")


mightyMyT=function(i) t.test(X[i,1:10],X[i,11:20],var.equal=T)$p.value
Pvals=mapply(mightyMyT,1:50)

# let's look at p-values for the truly DE genes:
Pvals[1:10]

# note: with Bonferroni correction
50*Pvals[1:10]

# let's look at p-values..
plot(1:50,Pvals,col=GeneGroup)

# Let's consider a cut-off for significance of 0.05
abline(h=0.05,col="blue")

# factor object for significant genes
SigFactor=factor(Pvals<=0.05, levels=c(T,F))
levels(SigFactor)=c("Significant", "Non-Sig")

# creating a table
table(SigFactor,GeneGroup)

# test association of group with significance using Fisher's Exact Test
fisher.test(SigFactor,GeneGroup)



