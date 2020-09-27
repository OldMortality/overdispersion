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
#   29/07/2020: Complete rewrite

# for testing only:##########
population = "Negbin"
phi = 3
n=100
b = c(-3,3)
## number of simulations
N = 1000
#############################

# global constants
CONF_LEVEL <- 0.95
ALPHA_LOWER <- 0.025
ALPHA_UPPER <- 0.975
# number of bootstrap samples
N.BOOTS = 10





SEQ.30 <- seq(1,30)
 

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
  
  return(sample(the.sequence,N.BOOTS*n,replace = T))  
}

sample.index.30 <-   matrix(create.samples(30),  nrow=1,ncol=30*N.BOOTS)

# profvis is for profiling (timing)
#profvis({

intervals.vgam <- function(x,y) {
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
  return(c(ci[1],ci[2]))
} 


intervals.vgam2 <- function(x,y) {
      a3 <- vglm(y~x,family=negbinomial(parallel=T,zero=NULL))
      # method = "wald
      # method = "profile"
      ci <- confint(a3,level = CONF_LEVEL)[3,]
      
      # estimate for phi
      #phihat3<-Confint.nb1(a3)$phi0
      # CI for phi
      phihat3ci<-Confint.nb1(a3,level=CONF_LEVEL)$CI.phi0
      
      
  
  return(c(ci[1],ci[2]))
} 


###################################################
#
# doCalculations
#
# Computes everything for 1 combination of population, phi,
#   sample size (n), beta, number of simulations (N)
#
#
doCalculations <- function(population,phi,n,b1,b2,N) {
  
  x <- seq(0,1,length.out=30)
  
  p <- 2
  eta <- b1 + b2*x
  mu <- exp(eta)
  nu <- phi-1
  mn <- mu/nu
  w <-  (mn+1)/mn
  mx <- log(mu)-0.5*log(w)
  sx <- sqrt(log(w))

  y <- rnegbin(n,mu,mn) 
  intervals <- intervals.vgam2(x=x,y=y)
  #intervalsT <- intervals.T(sim)
    
  # this does not help 
  vg.low <- intervals[1]
  vg.low <- unname(intervals[1])
  result <- c(vg.low)
  result.works <- -0.6908811
  return(c(result, result.works))
}




#set.seed(10)
N <- 10
phis <- c(3,5) # 5
pops <- c("Negbin" )
beta1 <- c(-3)
ns <-  c(30)
combs <- expand.grid(phis,pops,beta1,ns)
colnames(combs) <- c('phi','pop','beta1','n')
combs$beta2 <- NA
combs[which(combs$beta1==-3),'beta2'] <- 3
combs

# we use a CPU cores for every combination
no_cores= dim(combs)[1]
registerDoParallel(makeCluster(no_cores))


# All calculations, for writing to csv file. This will be
#   vector of strings: Each line is the comma separated results for
#   one set of parameters (population,phi,n,beta)
set.seed(10)

calcs <- foreach(i=1, .combine = rbind, .verbose=T ,
                 .packages=c('MASS','purrr','VGAM') )  %dopar% {
    doCalculations(population = combs$pop[i],
                 phi= combs$phi[i],
                 n = combs$n[i],
                 b1 = combs$beta1[i],
                 b2 = combs$beta2[i],
                 N = 1)
}
calcs



