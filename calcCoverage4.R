library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"
library(VGAM)

#library(profvis)
#library(compiler)

getRoots <- function(xx,s,q95,n) {
  A <- sum(s)^2 + n^2 + 2 * n * sum(s) - q95*sum((s-mean(s))^2)
  B <- -2 * sum(xx)*sum(s)-2*n*sum(xx)+
    2*q95*sum((xx-mean(xx))*(s-mean(s)))
  C <- (sum(xx)^2)-q95*sum((xx-mean(xx))^2)
  root1 <- (-B - sqrt(B^2-4*A*C))/(2*A)
  root2 <- (-B + sqrt(B^2-4*A*C))/(2*A)
  result <- c(root1,root2)
  return(result)
}




##
## This is the one with the estimating equations bootstrap
##
##

# for testing:
population = "Negbin"
phi = 1
n=30
b = c(2.3,0.7)
N = 100
#res <- calcCoverage(population,phi,n,b,N)

  
##
## returns the coverage rate.
##

seq.30 <- seq(1,30)
seq.100 <- seq(1,100)
seq.1000 <- seq(1,1000)


calcCoverage <- function(population,phi,n,b,N) {

  
  if (phi==1) {
    if (population != 'Negbin') {
      return(NA)
    }
  }
  
  seq.n <- NA
  if (n==30) {
    seq.n <- seq.30
  } else {
    if (n==100) {
      seq.n <- seq.100
    } else {
      if (n==1000) {
        seq.n <- seq.1000
      }
    }
  }
  
  
  
  
  chis <- vector()
  sumz <- vector()
  
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
  
  
  boots.low <- vector()
  boots.upp <- vector()
  boots.err <- rep(NA,N)

  # estimates of phi
  phihats2 <- vector()
  
  print(date())
  for (sim in 1:N) {
    if (sim %% 100 == 0) { print(paste(sim,date(),sep=' ')) }
    # set true population
    {
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
    } 
    muhat<-fitted(glm(y~x,family="poisson"))
    P<-sum((y-muhat)^2/muhat)
    sbar<-mean((y-muhat)/muhat) 
    phihat1<-P/(n-p)
    phihat2<-phihat1/(1+sbar) 
    phihats2[sim] <- phihat2
    
    # 
    df <- n-p
    
    xx <- (y-muhat)^2/muhat
    s <- (y-muhat)/muhat
    ## phi or phihat?
    z <- xx - phihat2*(1+s)
    z.bar <- mean(z)
    num <- (sum(z))^2
    denom <- sum((z-z.bar)^2)
    sumz[sim] <- sum(z)^2
    
    chis[sim] <- num/denom

    N.boots <- 10000
    q.stars <- vector()
    for (b in 1:N.boots) {
      boots.index <- sampleComp(seq(1,n),replace=T)
      z.star <- z[boots.index]
      z.star.bar <- mean(z.star)
      
      q.stars[b] <- (sum(z.star)^2)/sum((z.star-z.star.bar)^2)
    }
    q95 <- quantile(q.stars,0.95,na.rm=T)
    
    
    # solve equation to get lower and upper limits of CI
    
    roots <- getRoots(xx,s,q95,n)

    boots.low[sim] <- roots[1]
    boots.upp[sim] <- roots[2]
    boots.err[sim] <- ( roots[1] > phi ) |
                      ( roots[2] < phi )
    
    
  }  
  print(date())
  result <- 1-sum(boots.err,na.rm=T)/length(which(!is.na(boots.err)))
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
        print(paste('counter: ',counter))
        coverage <- calcCoverage(population=pop,
                                 phi=phi,
                                 n=n,
                                 b=betas[b,],
                                 N=N)
        parms <- paste(pop,phi,n,betas[b,1],betas[b,2],N,coverage,sep=',')
        results[counter] <- parms
        print(parms)
        
      }
    }
  }
}
write.csv(results,'results-ee.csv',row.names=F,quote=F)






