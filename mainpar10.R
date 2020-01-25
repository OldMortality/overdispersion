library(MASS)     # for "rnegbin"
#library(rmutil)   # for "rinvgauss"
#library(VGAM)
  
library(foreach)
library(doParallel)
library(purrr)

# Program to simulate CI's for the beta's/
# This one refits the model, for a single population/n,etc
# 
#
# Change log
#   22/11/2019: Created this program.
#   05/12/2019: Corrected bug where did not take sqrt of phi
#   08/12/2019: Keeping CIs for phi.
#   22/01/2020: renamed mainpar10. This one does all combinations in parallel.
#             : bootstrap is through resampling the residuals, not refitting model.
#

rm(list=ls(all=TRUE))
setwd('~/sims')
# global constants
CONF_LEVEL <- 0.9
ALPHA_LOWER <- 0.05
ALPHA_UPPER <- 0.95
# number of bootstrap samples
N.BOOTS = 10000
#
# Returns the proportion of elements in CI,
#   which do not include value. Elements containing
#   NA are disregarded. Returns NaN if all elements are NA,
#   but we hold that for impossible.
calcExcludesValue <- function(value,low,upp) {
  return(sum(low > value | upp < value,na.rm = T)/length(which(!is.na(low))))
}

# returns median width of CIs
calcMedianWidth <- function(lower,upper) {
  return(median(upper-lower,na.rm = T))
}



SEQ.30 <- seq(1,30)
SEQ.100 <- seq(1,100)
SEQ.1000 <- seq(1,1000)

# we use 54 CPU cores.
no_cores= 54
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

#X <- matrix(rep(x,N.BOOTS),byrow=T,nrow=N.BOOTS)

# profvis is for profiling (timing)
#profvis({
#t <- doCalculations(population,phi,n,b,N)

## Get true population.
##    
getY <- function(phi,population,n,mu,mn,nu) {
  
  y <- rep(0,n)
  while (sum(y) == 0) {
    # get a sample, until we get one, which is not all zeros.
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
  return(y)
}


# function to fit Poisson model, used in bootstraps
#   y.x is a vector of both the observations y and 
#   the predictors x. First come the y's, then the x's.
#   
#
fitit <- function(y.x,df,n) {
  model <- glm(y.x[1:n]~y.x[(n+1):(2*n)],family="poisson")
  phihat <- calcPhihat(fitted(model),y.x[1:n],df)
  return(data.frame(phihat))
}


calcPhihat <- function(muhat,y,df) {
  P <-sum((y-muhat)^2/muhat)
  sbar<-mean((y-muhat)/muhat) 
  phihat<-(P/df)/(1+sbar) 
}

###################################################
#
# doCalculations
#
# Computes everything for 1 combination of population, phi,
#   sample size (n), beta, number of simulations (N)
#
# Returns: 


doCalculations <- function(population,phi,n,b1,b2,N) {
  
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

  parms <- paste(population,phi,n,b1,b2,N,sep=',')
  print(parms)
  
  beta<-c(b1,b2)
  p<-length(beta)
  eta<-beta[1]+beta[2]*x
  mu<-exp(eta)
  nu <- phi-1
  mn<-mu/nu
  w<-(mn+1)/mn
  mx<-log(mu)-0.5*log(w)
  sx<-sqrt(log(w))
  
  
    
  ## vectors to hold CIs
  {
    boots.phi.int <- vector()
    
    phi.low <- vector()
    phi.upp <- vector()
    phi.err <- vector()
    
    beta.low <- vector()
    beta.upp <- vector()
    beta.err <- vector() 
  }
  
  #isInInterval(2,c(1,3))
  isInInterval <- function(value,interval) {
    return(value > interval[1] & value < interval[2])
  }
  
  #  
  # Start simulations
  #
  for (sim in 1:N) {
    
    # show heartbeat
    if (sim %% 500 == 0) { print(sim) }
    print(sim)    
    # set true population
    y <- getY(phi,population,n,mu,mn,nu)
    
    # fit Poisson model
    model <- glm(y~x,family="poisson")
    muhat<-fitted(model)
    P<-sum((y-muhat)^2/muhat)
    sbar<-mean((y-muhat)/muhat) 
    phihat<-(P/(n-p))/(1+sbar) 
    df <- n-p
    coef.x <- summary(model)$coefficients[2,1]
    coef.se <- summary(model)$coefficients[2,2]
    
    ##
    ## Bootstrap, resample residuals.
    ##
    do.bootstrap <- T
    
    if (do.bootstrap == T) {
      
      e.squared <- (y-muhat)^2/muhat
      # split this into a matrix with 1 row for each bootstrap
      ss <- matrix(
        data=e.squared[sample.index],byrow=T,nrow=N.BOOTS,ncol=n)
      
      
      
      # 1 row for each bootstrap 
      phihat.boots <- rowSums(ss)/df/(1+sbar)
      
      # 05/12/19 keeping CI for phi bootstrap, as opposed to CI for beta
      boots.phi.int <- phihat^2  / quantile(phihat.boots,c(ALPHA_UPPER,ALPHA_LOWER)) 

      ts <- rnorm(N.BOOTS) * sqrt(phihat/phihat.boots)
      # get t-values 
      t.low <- quantile(ts,ALPHA_LOWER)
      t.upp <- quantile(ts,ALPHA_UPPER)
      
      # ci for beta
      beta.low[sim] <- coef.x - t.upp * sqrt(phihat) * coef.se
      beta.upp[sim] <- coef.x - t.low * sqrt(phihat) * coef.se
      beta.err[sim] <- !isInInterval(beta[2],c(beta.low[sim],beta.upp[sim]))
      # ci for phi
      phi.low[sim] <- boots.phi.int[1]
      phi.upp[sim] <- boots.phi.int[2]
      phi.err[sim] <- !isInInterval(phi,c(phi.low[sim],phi.upp[sim]))
      
    }	  
   
  }
  
  ##
  ## end of simulations
  ##
  
  
  # combine coverage rates
  #
  error.phi <- sum(phi.err)/N
  error.beta <- sum(beta.err,na.rm=T)/length(which(!is.na(beta.err)))
  # convert error rates to coverage
  result.coverage <- c(1-error.phi,1-error.beta)
  
  # combine power calculations
  #
  phi.power <- calcExcludesValue(1,phi.low,phi.upp)
  beta.power <- calcExcludesValue(0,beta.low,beta.upp)
  result.power <- c(phi.power,beta.power)
  
  # combine CI median widths
  #
  phi.width <- calcMedianWidth(phi.low,phi.upp)
  beta.width <- calcMedianWidth(beta.low,beta.upp)
  result.width <- c(phi.width,beta.width)
  
  # gather the rates of NAs, ie. the rate of the method not working
  beta.na <- length(which(is.na(beta.err)))/N
  result.na <- c(beta.na)
  
  
  # combine all results. This one will go to the csv file.
  result <- c(parms, result.coverage,result.power,result.width,result.na)
  return(result)
}

###################################################
##
## end of doCalculations
##
###################################################
# profis })


print('started')
print(Sys.time())
set.seed(10)
N <- 10000
phis <- c(1,3,5)
pops <- c("Negbin" ,"Neyman","Poisson lognormal")
beta1 <- c(-3,0.1,2.3)
ns <-  c(30,100) 
combs <- expand.grid(phis,pops,beta1,ns)
colnames(combs) <- c('phi','pop','beta1','n')
combs$beta2 <- NA
combs[which(combs$beta1==-3),'beta2'] <- 3
combs[which(combs$beta1==0.1),'beta2'] <- 2.2
combs[which(combs$beta1==2.3),'beta2'] <-0.7
combs

# All calculations, for writing to csv file. This will be
#   vector of strings: Each line is the comma separated results for
#   one set of parameters (population,phi,n,beta)

calcs <- foreach(i=1:dim(combs)[1], .combine = rbind, .packages=c('MASS','purrr'))  %dopar% {
  doCalculations(population = combs$pop[i],
                 phi= combs$phi[i],
                 n = combs$n[i],
                 b1 = combs$beta1[i],
                 b2 = combs$beta2[i],
                 N = N)
}  

write.csv(calcs,'results-mainpar10.csv',row.names=F,quote=F)
print('finished')
print(Sys.time())