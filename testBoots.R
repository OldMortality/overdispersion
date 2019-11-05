library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"
library(VGAM)
#
library(proxfvis)   
library(foreach)
library(doParallel)

#
# Calculate CI's for phi, based on 6 methods:
#   Chisquare distribution, based on phihat2
#   Gamma distribution, also based on phihat2
#   Bootstrap, based on phihat1.
#   Bootstrap, based on phihat2.
#   VGAM restricted regression
#   Estimating Equations.
#   
#   
#
#
# Change log
#   <      : parallel version, using 4 CPU cores
#          : sample.index out of the inner loop
#   21 July: confidence level now 90% (was 95%)
#          : Added interval for bootstrap of phi1hat (the one not divided by 1+s.bar)
#          : Vectorised bootstraps (also at EE)
#   08 Aug : Add code to calculate median CI widths, and proportion
#          :    of the CI > 1
#   11 Aug : Calculate and return the NA-rates, i.e. the proportion of 
#          :    times the method does not appear to work.
#
#

# for testing only:##########
population = "Negbin"
phi = 1
n=30
b = c(-3,3)
## number of simulations
N = 1000
#############################


# global constants
CONF_LEVEL <- 0.9
ALPHA_LOWER <- 0.05
ALPHA_UPPER <- 0.95
# number of bootstrap samples
N.BOOTS = 10000




SEQ.30 <- seq(1,30)
SEQ.100 <- seq(1,100)
SEQ.1000 <- seq(1,1000)

# we use 4 CPU cores, one for each value of phi
no_cores= 4
registerDoParallel(makeCluster(no_cores))
##
## 
##
# save some time by setting these up only once
X.30 <- seq(0,1,length.out=30)
X.100 <- seq(0,1,length.out=100)
X.1000 <- seq(0,1,length.out=1000)


create.samples <- function(n) {
  if (n==30) {
    the.sequence <- SEQ.30
  }
  if (n==100) {
    the.sequence <- SEQ.100
  }
  if (n==1000) {
      the.sequence <- SEQ.1000
  }
  return(sample(the.sequence,N.BOOTS*n,replace = T))  
}

sample.index.30 <-   matrix(create.samples(30),  nrow=1,ncol=30*N.BOOTS)
sample.index.100 <-  matrix(create.samples(100), nrow=1,ncol=100*N.BOOTS)
sample.index.1000 <- matrix(create.samples(1000),nrow=1,ncol=1000*N.BOOTS)

# profvis is for profiling (timing)
profvis({
  t <- doCalculations(population,phi,n,b,N)
})

###################################################
#
# doCalculations
#
# Computes everything for 1 combination of population, phi,
#   sample size (n), beta, number of simulations (N)
#
# Returns: a vector of 24. There are 6 methods for working out
#   a confidence interval for phi, so we have:
#   6 coverage rates. 
#   6 proportions with phi.lower > 1
#   6 median widths of the CI for phi
#   6 rates of method not working.
#
doCalculations <- function(population,phi,n,b,N) {
  
  sample.index <- NA
  if (n==30) { 
    #sequ=SEQ.30 
    x <- X.30
    sample.index <- sample.index.30
  }
  if (n==100) { 
    #sequ=SEQ.100
    x <- X.100
    sample.index <- sample.index.100
  }
  if (n==1000) { 
    #sequ=SEQ.1000
    x <- X.1000
    sample.index <- sample.index.1000
  }
   
  
  # skip duplicates: For phi=1, the true population
  #   is always Poisson.
  if (phi==1) {
    if (population != 'Negbin') {
      return(NA)
    }
  }

  parms <- paste(population,phi,n,b[1],b[2],N,sep=',')
  print(parms)
  set.seed(10)
  beta<-c(b[1],b[2])
  p<-length(beta)
  eta<-beta[1]+beta[2]*x
  mu<-exp(eta)
  nu <- phi-1
  mn<-mu/nu
  w<-(mn+1)/mn
  mx<-log(mu)-0.5*log(w)
  sx<-sqrt(log(w))

  # error rates of CIs
  chisq.err <- vector()
  gamma.err <- vector()
  boots1.err <- vector() # based on phi1.hat
  boots2.err <- vector() # based on phi2.hat (phi1.hat / (1+s.bar))
  vgam.err <- vector()	    
  ee.err <- vector() 

  # estimates of phi
  
  phi.hats.boots <- vector()
  # used by EE
  chis <- vector()
  sumz <- vector() 

  
  ## vectors to hold CIs
  {
  chisq.low <- vector()
  chisq.upp <- vector()
  gamma.low <- vector()
  gamma.upp <- vector()
  boots1.low <- vector()
  boots1.upp <- vector()
  boots2.low <- vector()
  boots2.upp <- vector()
  vgam.low <- vector()
  vgam.upp <- vector()
  ee.low <- vector()
  ee.upp <- vector()
  }
  
  #  
  # Start simulations
  #
  for (sim in 1:N) {
    
    # show heartbeat
    if (sim %% 500 == 0) { print(sim) }
    
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
         
      }
    }  
    
    
    # fit Poisson model
    muhat<-fitted(glm(y~x,family="poisson"))
    P<-sum((y-muhat)^2/muhat)
    sbar<-mean((y-muhat)/muhat) 
    phihat1<-P/(n-p)
    phihat2<-phihat1/(1+sbar) 
    #phihats2[sim] <- phihat2
    
    
    ##
    ## Bootstrap
    ##
      
      phihat1.boots <- vector()
      phihat2.boots <- vector()
      e. <- (y-muhat)/sqrt(muhat)
      
                #for (i in 1:N.BOOTS) {
                #  subsample.index <- sample.index[i,] 
                #  e <- e.[subsample.index]
                #  P.b<-sum(e^2)
                #  phi.1.b <- P.b/df
                #  phi.2.b <- phi.1.b / (1+sbar)
                #  phihat1.boots[i] <- phi.1.b
                #  phihat2.boots[i] <- phi.2.b
                #}
      
      e.squared <- e.^2
      # split this into a matrix with 1 row for each bootstrap
      ss <- matrix(
        data=e.squared[sample.index],byrow=T,nrow=N.BOOTS,ncol=n)
      
      # 1 row for each bootstrap 
      phihat1.boots <- rowSums(ss)/df
      phihat2.boots <- phihat1.boots/(1+sbar)
      
  
      # phihat1
      boots1.int <- (df /phihat1) *quantile(phihat1.boots,c(ALPHA_LOWER,ALPHA_UPPER)) 
      boots1.low[sim] <- (df * phihat1 ) / boots1.int[2]
      boots1.upp[sim] <- (df * phihat1 ) / boots1.int[1]
      boots1.err[sim] <- (phi < boots1.low[sim]) |
        (phi > boots1.upp[sim])
      
      # phihat2
      boots2.int <- (df /phihat2) *quantile(phihat2.boots,c(ALPHA_LOWER,ALPHA_UPPER)) 
      boots2.low[sim] <- (df * phihat2 ) / boots2.int[2]
      boots2.upp[sim] <- (df * phihat2 ) / boots2.int[1]
      boots2.err[sim] <- (phi < boots2.low[sim]) |
        (phi > boots2.upp[sim])
      if (is.na(boots2.err[sim])) {
        print(boots2.err[sim])
        print(boots2.int[1])
        print(boots2.int[2])
        quit(save="yes")
      }
    }	  
    
  }  
  which(is.na(boots2.err))
  
  ##
  ## end of simulations
  ##
  
  