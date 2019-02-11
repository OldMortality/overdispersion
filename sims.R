library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"

# calculate CI's for phi, based on chisquare distribution
#   with phihat1 and phihat2, and on Gamma distribution
#   (phihat2 only), and fitting only Negative binomial 
#   for now.
#   Plots first 100 of the simulations, with error rates.
#

for (phi in c(2)) {


{
n<-30

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

N <- 2000
# bounds of CI for chisq and phi1

lower.1 <- vector()
upper.1 <- vector()
# bounds of CI for chisq and phi2
lower.2 <- vector()
upper.2 <- vector()
# bounds of CI for gamma and phi2
gamma.lower <- vector()
gamma.upper <- vector()

# estimates of phi
phihats1 <- vector()
phihats2 <- vector()

# does CI contain true phi?
err.1 <- vector()
err.2 <- vector()
gamma.err <- vector()
ks <- vector()
thetas <- vector()
}

for (sim in 1:N) { 
  if (phi==1) { 
    y<-rpois(n,mu) } 
  else { 
    y<-rnegbin(n,mu,mn) # negative binomial
    
    #y<-rpois(n,rpois(n,mn)*nu)             # Neyman Type A
    #y<-rpois(n,rlnorm(n,mx,sx))            # Poisson lognormal
    #y<-rpois(n,rinvgauss(n,1,1/mn)*mu)     # Poisson inverse Gaussian
  
    }  
  muhat<-fitted(glm(y~x,family="poisson"))
  P<-sum((y-muhat)^2/muhat)
  sbar<-mean((y-muhat)/muhat) 
  phihat1<-P/(n-p)
  phihat2<-phihat1/(1+sbar) 
  phihats1[sim] <- phihat1
  phihats2[sim] <- phihat2
  
  # chisquare CI
  df <- n-p
  lower.1[sim] <- df * phihat1  / qchisq(0.975,df=df)
  upper.1[sim] <- df * phihat1  / qchisq(0.025,df=df)
  lower.2[sim] <- df * phihat2  / qchisq(0.975,df=df)
  upper.2[sim] <- df * phihat2  / qchisq(0.025,df=df)
  
  err.1[sim] <- 0
  if (phi < lower.1[sim] | phi > upper.1[sim]) {
    err.1[sim] <-  1
  }
  err.2[sim] <- 0
  if (phi < lower.2[sim] | phi > upper.2[sim]) {
    err.2[sim] <-  1
  }
  ## Gamma CI. only for phihat2
  e <- (y-muhat)/sqrt(phihat2 * muhat)
  rho3 <- e^3
  rho4 <- e^4 - 3
  alpha.bar <- mean(rho4 + phihat2/muhat - 2 * sqrt(phihat2/muhat) * rho3)
  theta <- 2 + alpha.bar
  k <- df/(2+alpha.bar)
  thetas[sim] <- theta
  ks[sim] <- k
  #gamma.lower[sim] <- df * phihat2  / qgamma(0.975,
  #                                         shape=k,
  #                                         scale=theta)
  ##gamma.upper[sim] <- df * phihat2  / qgamma(0.025,
  #                                         shape=k,
  #                                         scale=theta)
  #gamma.err[sim] <- (phi < gamma.lower[sim]) |
  #             (phi > gamma.upper[sim])
  
}

# error rates
e1 <- sum(err.1)/N
e2 <- sum(err.2)/N
e3 <- sum(gamma.err,na.rm=T)/length(which(!is.na(gamma.err)))
paste(e1,e2,e3,sep='   ')

# plot CIs, show first plotN simulations
plotCIs <- function(low,upp,err,phihat,truePhi,plotN,main,
                    errRate) {
  plot('',xlim=c(0,plotN),ylim=c(0,4*truePhi),
       main=paste(main,'true phi=',truePhi,sep=' '),
       ylab="95% CI")

  abline(h=truePhi,col='red')
  for (i in 1:plotN) {
    col <- 'black'
    if (err[i] == 1) {
      col <- 'red'
    }
    segments(i,low[i],i,upp[i],col=col)
    points(i,phihat[i],pch='.',cex=3)
    text(60,3.5 * truePhi,paste("error rate: ",errRate))
  }  
}  
  
#par(mfrow=c(1,3))
plotCIs(lower.1,upper.1,err.1,phihats1,phi,plotN=100,
        main="Chisq, phihat1",e1)
#plotCIs(lower.2,upper.2,err.2,phihats2,phi,plotN=100,
#        main="Chisq, phihat2",e2)
#plotCIs(gamma.lower,gamma.upper,gamma.err,phihats2,phi,plotN=100,
#        main="Gamma, phihat2",e3)

#par(mfrow=c(1,3))

median(upper.1-lower.1)
median(upper.2-lower.2)
median(gamma.upper-gamma.lower,na.rm=T)

# plot estimates of phihat's


#par(mfrow=c(2,2))
{
  print(paste(phi,
              round(mean((n-p)*phihats1/phi),2),
              round(var((n-p)*phihats1/phi),2),
              round(mean((n-p)*phihats2/phi),2),
              round(var((n-p)*phihats2/phi),2),
              sep='   '))
  
par(mfrow=c(1,1))
main=paste('phi=',phi)
hist((n-p)*phihats2/phi,probability = T,
     main=main,ylim=c(0,0.4),
     xlim=c(0,2 * n),100)
c<- seq(0,qchisq(0.9999,df),0.1)
y1<-dchisq(c,df=n-p)
lines(c,y1,col='red') 
#y2<- dgamma(c,shape=mean(ks),scale=mean(thetas))
#lines(c,y2,col='blue') 
#legend(50,0.08, 
#       c("gamma","chisq"),
#       lty=c(1,1),
#       lwd=c(2.5,2.5),
#       col=c("blue","red"))
}

}

#mean(phihats1)
mean(phihats2)
hist(phihats2)
abline(v=phi,col='red')
abline(v=mean(phihats2,col='blue'))
#par(mfrow=c(2,2))

##(n / df) * phi
#var(phihats1)

# these should be roughly equal:
#(1/(df)) * (sum(y^2/muhat)-sum(muhat))
#phihat2
#


