library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"

n<-30
phi<-2
beta<-c(1,1)
x<-seq(0,1,length.out=n)

p<-length(beta)
eta<-beta[1]+beta[2]*x
mu<-exp(eta)

nu<-phi-1
mn<-mu/nu
w<-(mn+1)/mn
mx<-log(mu)-0.5*log(w)
sx<-sqrt(log(w))

if (phi==1) { y<-rpois(n,mu) } else
{ y<-rnegbin(n,mu,mn) }                   # Negative binomial
{ y<-rpois(n,rpois(n,mn)*nu) }            # Neyman Type A
{ y<-rpois(n,rlnorm(n,mx,sx)) }           # Poisson lognormal
{ y<-rpois(n,rinvgauss(n,1,1/mn)*mu) }    # Poisson inverse Gaussian

muhat<-fitted(glm(y~x,family="poisson"))
P<-sum((y-muhat)^2/muhat)
sbar<-mean((y-muhat)/muhat)

phihat1<-P/(n-p)
phihat2<-phihat1/(1+sbar)
