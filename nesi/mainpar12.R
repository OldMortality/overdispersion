rm(list=ls())

library(MASS)     # for "rnegbin"
library(VGAM)
#
#library(profvis)   
library(foreach)
library(doParallel)

# Program to simulate CI's for beta2 and phi's using 
#   traditional Chisq method and vgam, 
#   fastboot method
#   hybrid method: 
#     "If the lower limit of the chi-squared interval is less than or equal to 1, 
#      use the chi-squared interval; otherwise use the bootstrap interval."
#
# Change log
#   28/07/2020: Created this program, by modifying mainpar10.R and mainpar11.R
#

# for testing only:##########
population = "Negbin"
phi = 3
n=30
b = c(-3,3)
## number of simulations
N = 1000
#############################

# global constants
CONF_LEVEL <- 0.95
ALPHA_LOWER <- 0.025
ALPHA_UPPER <- 0.975
# number of bootstrap samples
N.BOOTS = 10000


# 29/07/20 Replaced calcExcludesZero and calcExcludesOne
#
# Returns the proportion of elements in CI,
#   which do not include the value. Elements containing
#   NA are disregarded. Returns NaN if all elements are NA,
#   but we hold that for impossible.
calcExcludesValue <- function(val,low,upp) {
  return(sum(low > val | upp < val,na.rm = T)/length(which(!is.na(low))))
}

# returns median width of CIs
calcMedianWidth <- function(lower,upper) {
  return(median(upper-lower,na.rm = T))
}

isInInterval <- function(value,low,upp) {
  return(value > low & value < upp)
}

# if lower bound of chisquare interval < 1, use chisquare, else use bootstrap
calcHybridInterval <- function(chi.low,chi.upp,b.low,b.upp) {
  # use bootstrap by default
  h.low <- b.low
  h.upp <- b.upp
  # which rows have chisquare lower bound < 1
  ix <- which(chi.low < 1)
  # for those, use chisquare 
  h.low[ix] <- chi.low[ix]
  h.upp[ix] <- chi.upp[ix]
  return(list(h.low=h.low,h.upp=h.upp))
}


SEQ.30 <- seq(1,30)
SEQ.100 <- seq(1,100)
SEQ.1000 <- seq(1,1000)

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
#profvis({
#  t3 <- doCalculations(population,phi,n,b,N)

## 
getY <- function(phi,population,n,mu,mn,nu)
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
  return(y)
}
##

###################################################
#
# doCalculations
#
# Computes everything for 1 combination of population, phi,
#   sample size (n), beta, number of simulations (N)
#
# Returns: a vector of 12. There are 3 methods for working out
#   a confidence interval for beta, so we have:
#   3 coverage rates. 
#   3 proportions with CI for beta excluding zero
#   3 median widths of the CI for beta
#   3 rates of method not working.
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
  
  beta<-c(b[1],b[2])
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
  chisq.low <- vector()
  chisq.upp <- vector()
  chisq.phi.low <- vector()
  chisq.phi.upp <- vector()
  
  boot.low <- vector()
  boot.upp <- vector()
  boot.phi.low <- vector()
  boot.phi.upp <- vector()
  
  
  hybrid.beta.low <- vector()
  hybrid.beta.upp <- vector()
  
  vgam.low <- vector()
  vgam.upp <- vector()
  vgam.phi.low <- vector()
  vgam.phi.upp <- vector()
  
  # error rates of CIs
  chisq.err <- vector()
  chisq.phi.err <- vector()
  vgam.err <- vector()	    
  vgam.phi.err <- vector()
  boot.err <- vector()
  boot.phi.err <- vector()
  
  }
  
  #  
  # Start simulations
  #
  for (sim in 1:N) {
    
    # show heartbeat
    if (sim %% 500 == 0) { print(sim) }
    
    # set true population, but skip those simulations
    #    in which all observations are zero.
    y <- rep(0,n)
    while (sum(y) == 0) {
      y <- getY(phi,population,n,mu,mn,nu)
    }
    
    # fit Poisson model
    model <- glm(y~x,family="poisson")
    muhat<-fitted(model)
    P<-sum((y-muhat)^2/muhat)
    sbar<-mean((y-muhat)/muhat) 
    phihat<-(P/(n-p))/(1+sbar) 
    
    ##
    ## chisquare 
    ##
    df <- n-p
    summ <- summary(model)
    coef.x <- summ$coefficients[2,1]
    coef.se <- summ$coefficients[2,2]
    chisq.low[sim] <- coef.x - qt(ALPHA_UPPER,df=df) * phihat * coef.se
    chisq.upp[sim] <- coef.x - qt(ALPHA_LOWER,df=df) * phihat * coef.se
    
    # 29/07/20 using isInInterval()
    chisq.err[sim] <- !isInInterval(value=beta[2],low=chisq.low[sim],upp=chisq.upp[sim])
    
    phi.int <- df * phihat  / qchisq(df=df,c(ALPHA_UPPER,ALPHA_LOWER))
    chisq.phi.low[sim] <- phi.int[1]
    chisq.phi.upp[sim] <- phi.int[2]
    
    # 29/07/20 using isInInterval()
    chisq.phi.err[sim] <-  !isInInterval(value=phi,low=chisq.phi.low[sim],upp=chisq.phi.upp[sim])
    
    
    
    ## Fastboot
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
      boot.low[sim] <- coef.x - t.upp * sqrt(phihat) * coef.se
      boot.upp[sim] <- coef.x - t.low * sqrt(phihat) * coef.se
      boot.err[sim] <- !isInInterval(beta[2],boot.low[sim],boot.upp[sim])
      
      # ci for phi
      boot.phi.low[sim] <- boots.phi.int[1]
      boot.phi.upp[sim] <- boots.phi.int[2]
      boot.phi.err[sim] <- !isInInterval(phi,boot.phi.low[sim],boot.phi.upp[sim])
      
    }	  
    
    
    ##
    ## VGAM
    ##
    
    do.RR <- T

    if (do.RR) {
      phihat3ci <- c(NA,NA)
      ci <- c(NA,NA)
      tryCatch(
        {
          a3 <- vglm(y~x,family=negbinomial(parallel=T,zero=NULL))
          # method = "wald
          # method = "profile"
          ci <- confint(a3,level = CONF_LEVEL)[3,]
          
          # estimate for phi
          #phihat3<-Confint.nb1(a3)$phi0
          # CI for phi
          phihat3ci<-Confint.nb1(a3,level=CONF_LEVEL)$CI.phi0
          
          
        },
        error=function(cond) {
          #message(paste('rr',y,x,cond,sep=' '))
          return(NA)
        } ,
        warning=function(cond) {
          #ci <- confint(a3)[3,]
          return(NA)
        }
      
      )
      # beta
      vgam.low[sim] <- ci[1]
      vgam.upp[sim] <- ci[2]
      
      # phi
      vgam.phi.low[sim] <- phihat3ci[1]
      vgam.phi.upp[sim] <- phihat3ci[2]
      
      
      # this will be NA if vgam.upp[sim]==NA
      # 29/07/20 replaced by function isInInterval
      vgam.err[sim] <- !isInInterval(value=beta[2],low=ci[1],upp=ci[2])
      vgam.phi.err[sim] <- !isInInterval(value=phi,low=phihat3ci[1],upp=phihat3ci[2])
      
    }
   
  }
  
  ##
  ## end of simulations
  ##
  
  # calculate hybrid intervals
  h <- calcHybridInterval(chisq.low,chisq.upp,boot.low,boot.upp) 
  h.low <- h$h.low
  h.upp <- h$h.upp
  h.phi <- calcHybridInterval(chisq.phi.low,chisq.phi.upp,boot.phi.low,boot.phi.upp)
  h.phi.low <- h.phi$h.low
  h.phi.upp <- h.phi$h.upp
  
  # combine power calculations
  
  # beta
  chi.p1 <- calcExcludesValue(val=0,low=chisq.low,upp=chisq.upp)
  vga.p1 <- calcExcludesValue(val=0,low=vgam.low, upp=vgam.upp)
  boot.p1 <- calcExcludesValue(val=0,low=boot.low, upp=boot.upp)
  h.p1 <- calcExcludesValue(val=0,low=h.low, upp=h.upp)
  
  
  # phi
  chi.phi.p1 <- calcExcludesValue(val=1,low=chisq.phi.low,upp=chisq.phi.upp)
  vga.phi.p1 <- calcExcludesValue(val=1,low=vgam.phi.low,upp=vgam.phi.upp)
  boot.phi.p1 <- calcExcludesValue(val=1,low=boot.phi.low,upp=boot.phi.upp)
  h.phi.p1 <- calcExcludesValue(val=1,low=h.phi.low,upp=h.phi.upp)
  
  
  result.prop1 <- c(chi.p1,
                    vga.p1,
                    boot.p1,
                    h.p1,
                    chi.phi.p1,
                    vga.phi.p1,
                    boot.phi.p1,
                    h.p1)
  
  # combine coverage rates
  #
  errorChisq <- sum(chisq.err)/N
  errorChisq.phi <- sum(chisq.err)/N
  
  errorVGAM  <- sum(vgam.err,na.rm=T)/length(which(!is.na(vgam.err)))
  errorVGAM.phi  <- sum(vgam.phi.err,na.rm=T)/length(which(!is.na(vgam.phi.err)))
  
  errorBoot <- sum(boot.err)/N
  errorBoot.phi <- sum(boot.phi.err/N)
  
  error.h <- calcExcludesValue(b[2],h.low,h.upp)
  error.h.phi <- calcExcludesValue(phi,h.phi.low,h.phi.upp)
  
  # convert error rates to coverage
  # 29/07/20 added 1- for the final two parameters
  result.coverage <- c(1-errorChisq,
                       1-errorVGAM,
                       1-errorChisq.phi,
                       1-errorVGAM.phi,
                       1-errorBoot,
                       1-errorBoot.phi,
                       1-error.h,
                       1-error.h.phi)
  
  # combine CI median widths
  #
  chi.width <- calcMedianWidth(chisq.low,chisq.upp)
  chi.width.phi <- calcMedianWidth(chisq.phi.low,chisq.phi.upp)
  
  vga.width <- calcMedianWidth(vgam.low,vgam.upp)
  vga.width.phi <- calcMedianWidth(vgam.phi.low,vgam.phi.upp)
  
  boot.width <- calcMedianWidth(boot.low,boot.upp)
  boot.width.phi <- calcMedianWidth(boot.phi.low,boot.phi.upp)
  
  h.width <- calcMedianWidth(h.low,h.upp)
  h.width.phi <- calcMedianWidth(h.phi.low,h.phi.upp)
  
  result.width <- c(chi.width,
                    chi.width.phi,
                    vga.width,
                    vga.width.phi,
                    boot.width,
                    boot.width.phi,
                    h.width,
                    h.width.phi)
  
  # gather the rates of NAs, ie. the rate of the method not working
  chi.na <- 0 # always works
  vga.na <- length(which(is.na(vgam.err)))/N
  result.na <- c(chi.na, vga.na)
  
  
  # combine all results. This one will go to the csv file.
  result <- c(result.coverage,result.prop1,result.width,result.na)
  return(result)
}

###################################################
##
## end of doCalculations
##
###################################################
# profis })


print('started mainpar11.R')
print(Sys.time())
set.seed(10)
N <- 10000
phis <- c(1,3,5)
pops <- c("Negbin" ,"Neyman","Poisson lognormal")
beta1 <- c(-3,0.1,2.3)
ns <-  c(30,100,1000)
combs <- expand.grid(phis,pops,beta1,ns)
colnames(combs) <- c('phi','pop','beta1','n')
combs$beta2 <- NA
combs[which(combs$beta1==-3),'beta2'] <- 3
combs[which(combs$beta1==0.1),'beta2'] <- 2.2
combs[which(combs$beta1==2.3),'beta2'] <-0.7
combs

# we use a CPU cores for every combination
no_cores= dim(combs)[1]
registerDoParallel(makeCluster(no_cores))


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


write.csv(results,'~/results-mainpar12.csv',row.names=F,quote=F)
print('finished')
print(Sys.time())

