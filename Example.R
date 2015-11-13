rm(list=ls());

##load related functions
source("CorSigFunc.R")


###load data
cData=read.table("cData_5363.txt",header=TRUE,sep="\t")
DATA=cData[,3:(which(colnames(cData)=="Systematic.Name")-1)];

##TargetGenes
ind=1;##
corr=cor(x=t(DATA)[,ind], y =t(DATA), use = "all.obs", method = "pearson",
    quick = 0, nThreads = 0, verbose = 0, indent = 0);
#####
tau=0.6;
n=ncol(DATA)
pvalue=PCSignif_func(abs(corr),tau=tau,n=n)
##P values of correlations of a target gene with others
pvalue
