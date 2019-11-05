library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"

# calculate CI's for phi, based on chisquare distribution
#   with phihat1 and phihat2, and on Gamma distribution
#   (phihat2 only), and fitting only Negative binomial 
#   for now.
#   Plots first 100 of the simulations, with error rates.
#

calcCoverage <- function(population,phi,n,b,N) {

  parms <- paste(population,phi,n,b[1],b[2],N,sep=',')
  print(parms)
  set.seed(10)
  beta<-c(b[1],b[2])
  x<-seq(0,1,length.out=n)
  p<-length(beta)
  eta<-beta[1]+beta[2]*x
  mu<-exp(eta)
  print(phi)
  nu <- phi-1
  
  mn<-mu/nu
  w<-(mn+1)/mn
  mx<-log(mu)-0.5*log(w)
  sx<-sqrt(log(w))
  # bounds of CI for chisq and phi2
  lower.chisq <- vector()
  upper.chisq <- vector()
  # bounds of CI for gamma and phi2
  gamma.lower <- vector()
  gamma.upper <- vector()
    
  # estimates of phi
  phihats2 <- vector()
  # does CI contain true phi?
  err.2 <- vector()
  gamma.err <- vector()
  ks <- vector()
  
  for (sim in 1:N) { 
    # set true population
    if (phi==1) { 
      y<-rpois(n,mu) } 
    else { 
      if (population == "Negbin") {
        y<-rnegbin(n,mu,mn) 
      } 
      if (population == 'Neyman') {
        y<-rpois(n,rpois(n,mn)*nu)            
      }
      if (population == 'Poisson lognormal') {
        y<-rpois(n,rpois(n,mn)*nu)        
      }
      if (population == 'Poisson inverse Gaussian') {
        y<-rpois(n,rinvgauss(n,1,1/mn)*mu)     
      }
    }
      
    muhat<-fitted(glm(y~x,family="poisson"))
    P<-sum((y-muhat)^2/muhat)
    sbar<-mean((y-muhat)/muhat) 
    phihat1<-P/(n-p)
    phihat2<-phihat1/(1+sbar) 
    phihats2[sim] <- phihat2
    
    # chisquare CI
    df <- n-p
    lower.chisq[sim] <- df * phihat2  / qchisq(0.95,df=df)
    upper.chisq[sim] <- df * phihat2  / qchisq(0.05,df=df)
    
    err.2[sim] <- 0
    if (phi < lower.chisq[sim] | phi > upper.chisq[sim]) {
      err.2[sim] <-  1
    }
    
    #
    e <- y-muhat
    alpha3.hat <- (1/df) * sum( (e^3/muhat))
    alpha4.hat <- (1/df) * sum(e^4/muhat) - 3 * muhat* phihat2^2
    tau.i <- (alpha4.hat/phi^2 - 2 * alpha3.hat/phi + phi )
    tau.i <- (1/muhat) * tau.i
    tau = mean(tau.i)
    
    # work out S
    eta <- log(muhat)
    W <- diag(muhat)
    X1 <- rep(1,n)
    X2 <- x
    X <- cbind(X1,X2)
    Q <- X %*% solve(t(X) %*% W %*% X) %*% t(X)
    S <- sum(1/muhat) + n * sum(diag(Q))-sum(Q)
    
    # bias for R = df * phihat/phi
    bias.phi <-  -(alpha3.hat-phihat2^2)*S/n/df
    
    ktheta <- df + df * bias.phi / phihat2
    ktheta2 <- df * (2 + tau) 
    theta <- ktheta2 / ktheta
    k <- ktheta / theta
    
    gamma.lower[sim] <- df * phihat2  / qgamma(0.95,
                                             shape=k,
                                             scale=theta)
    
    gamma.upper[sim] <- df * phihat2  / qgamma(0.05,
                                               shape=k,
                                               scale=theta)
    gamma.err[sim] <- (phi < gamma.lower[sim]) |
                 (phi > gamma.upper[sim])
    
  }
  
  # error rates
  errorChisq <- sum(err.2)/N
  errorGamma <- sum(gamma.err,na.rm=T)/length(which(!is.na(gamma.err)))
  result <- c(1-errorChisq,1-errorGamma)
  return(result)
} # end of calcCoverage

 
pops <- c("Negbin","Neyman","Poisson lognormal","Poisson inverse Gaussian")
phis <- c(1,2,3,5)

betas <- rbind(c(-3,3),c(0.1,2.2),c(2.3,0.7))
ns <- c(30,100,1000)
N <- 10000
results <- vector()
counter <- 0
for (pop in pops) {
  for (phi in phis) {
    for (n in ns) {
      for (b in seq(1,3)) {
        counter <- counter + 1
        print(counter)
        coverage <- calcCoverage(population=pop,
                                 phi=phi,
                                 n=n,
                                 b=betas[b,],
                                 N=1)
        parms <- paste(pop,phi,n,betas[b,1],betas[b,2],N,coverage[1],coverage[2],sep=',')
        results[counter] <- parms
        
      }
    }
  }
}

write.csv(results,'~/Documents/overdispersion/results3.csv',
          row.names = F,quote=F)
f <- read.csv('~/Documents/overdispersion/results3.csv',skip=1)
dim(f)
head(f)



pops <- c("Negbin","Neyman","Poisson lognormal","Poisson inverse Gaussian")
coverage <- calcCoverage(population=pops[1],
                         phi=2,
                         n=0,
                         b=c(-3,3),
                         N=10000)

#tt <- read.csv('~/Documents/overdispersion/results3.csv',skip=1)
#head(tt)
#tt2 <- tt[,c(2:9)]
#tt.colnames=c("Population",'phi','n','b1','b2',"chisq.coverage","gamma.coverage")
#head(tt)
#options(digits=4)
#beta <- c(-3,3)
#beta <- c(0.1,2.2)
#beta <- c(2.3,0.7)

#eta<-beta[1]+beta[2]*x
#mu<-exp(eta)
#range(mu)
#hist(mu)

par(mfrow=c(3,1))
plotCIs(lower.chisq,upper.chisq,err.2,phihats2,phi,plotN=300,main='Chisq',1-result[1])
plotCIs(gamma.lower,gamma.upper,gamma.err,phihats2,phi,plotN=300,main='Gamma',1-result[2])
plotCIs(boots2.low,boots2.upp,boots2.err,phihat.boots2,phi,plotN=300,main='Bootstrap',1-result[3])
