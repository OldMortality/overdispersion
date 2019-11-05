# this data gives you VERY wide CI without warning
library(VGAM)
x <- seq(0,1,length.out=30)
#y = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,2,2,0,0,0,1,0,0,2)
y = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,2,2,0,0,0,1,0,0,2)
a3<-vglm(y~x,family=negbinomial(parallel=T,zero=NULL))
phihat3<-Confint.nb1(a3)$phi0
# CI for phi
phihat3ci<-Confint.nb1(a3,level=0.90)$CI.phi0
phihat3ci
