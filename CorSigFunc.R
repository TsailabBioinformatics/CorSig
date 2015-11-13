#################################
####CorSig v1 BY H-Q Wang and CJ TSAI, July 30, 2013
###################################
#############################
####1,The code works in R, and 
####2,the support data file "SD_N_RHO.txt" MUST 
####be copied to your work directory.
########################

###########define related functions
if(!require("WGCNA")) {
install.packages("WGCNA");
library("WGCNA");##load library
}

##z-transformation
fisherszfun<-function(x){0.5*(log(1+x)-log(1-x))}; #Fisher Z transformation;

###
######inputs: 
##r--an observed correlation of interest
##tau--a correlation threshold with default setting of 0.6
##n--Sample size
##sigmaz--standard deviation
######output:
##alpha--p-value

###calculating p-values for a observed r
PCSignif_func<-function(r,tau=0.6,n,sigmaz=NULL)
{
if(is.null(sigmaz))
{
if(n>=10000)
{sigmaz=(n-3)^-0.5;
}else
{
knownDat=read.table("SD_N_RHO.txt",sep="\t");
N=as.numeric(substring(colnames(knownDat),3));N
RHO=as.numeric(substring(rownames(knownDat),5));RHO
x=RHO;
y=as.numeric(knownDat[,which.min(abs(N-n))])##choose that closest to n

##########################################333
if(n%in%N&tau%in%RHO) {
##using loess unappliable to some cases
mod=loess(y~x);sigmaz=predict(mod,tau);sigmaz
} else {
##using linear model appliable to any case
mod=lm(y~x);
}

sigmaz=as.numeric(predict(mod,newdata=data.frame(x=c(tau))));
}
}
###improve func1
Z_1ma=(fisherszfun(r)-fisherszfun(tau))/sigmaz;##Z_ima>...
alpha=1-pnorm(Z_1ma)#
return(alpha);
}

