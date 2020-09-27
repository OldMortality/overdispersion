rm(list=ls())

library(MASS)     # for "rnegbin"
library(VGAM)
#
#library(profvis)   

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
#   12/08/2020: fixed bugs with beta parameters

# for testing only:##########
population = "Negbin"
phi = 3
n=100
b = c(-3,3)
## number of simulations
N = 5000
#############################

# global constants
CONF_LEVEL <- 0.95
ALPHA_LOWER <- 0.025
ALPHA_UPPER <- 0.975
# number of bootstrap samples
N.BOOTS = 10


# 29/07/20 Replaced calcExcludesZero and calcExcludesOne
#
# Returns the proportion of elements in CI,
#   which do not include the value. Elements containing
#   NA are disregarded. Returns NaN if all elements are NA,
#   but we hold that for impossible.
calcExcludesValue <- function(val,low,upp) {
  return(sum(low > val | upp < val,na.rm = T)/length(which(!is.na(low))))
}

calcIncludesValue <- function(val,low,upp) {
  return(1 - calcExcludesValue(val,low,upp))
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
##
## Set true population, but skip those simulations
##    in which all observations are zero.
##
getY <- function(phi,population,n,mu,mn,nu) {
  y <- rep(0,n)
  while (sum(y) == 0) {
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
##

 

intervals.chisq <- function(coef.x,coef.se,df,phihat) {
  beta.low <- coef.x - qt(ALPHA_UPPER,df=df) * phihat * coef.se
  beta.upp <- coef.x - qt(ALPHA_LOWER,df=df) * phihat * coef.se
  
  phi.int <- df * phihat  / qchisq(df=df,c(ALPHA_UPPER,ALPHA_LOWER))
  phi.low <- phi.int[1]
  phi.upp <- phi.int[2]
  return(c(beta.low,beta.upp,phi.low,phi.upp))
  
}

intervals.fastboot <- function(y,muhat,sample.index,n,df,sbar,phihat,coef.x,coef.se) {
  e.squared <- (y-muhat)^2/muhat
  # split this into a matrix with 1 row for each bootstrap
  ss <- matrix(
    data=e.squared[sample.index],byrow=T,nrow=N.BOOTS,ncol=n)
  # 1 row for each bootstrap 
  phihat.boots <- rowSums(ss)/df/(1+sbar)
  
  boots.phi.int <- phihat^2  / quantile(phihat.boots,c(ALPHA_UPPER,ALPHA_LOWER)) 
  ts <- rnorm(N.BOOTS) * sqrt(phihat/phihat.boots)
  # get t-values 
  t.low <- quantile(ts,ALPHA_LOWER)
  t.upp <- quantile(ts,ALPHA_UPPER)
  
  # ci for beta
  beta.low <- coef.x - t.upp * sqrt(phihat) * coef.se
  beta.upp <- coef.x - t.low * sqrt(phihat) * coef.se
  # ci for phi
  phi.low <- boots.phi.int[1]
  phi.upp <- boots.phi.int[2]
  return(c(beta.low,beta.upp,phi.low,phi.upp))
}

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
  return(c(ci[1],ci[2],phihat3ci[1],phihat3ci[2]))
} 
  

getSampleIndex <- function(n) {
  x = NA
  sample.index = NA
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
  return(list(x=x,sample.index=sample.index))
}

###########################doCal########################
#
# doCalculations
#
# Computes everything for 1 combination of population, phi,
#   sample size (n), beta, number of simulations (N)
#
#
doCalculations <- function(population,phi,n,b1,b2,N) {
  
  sample.index.x <- getSampleIndex(n)
  sample.index = sample.index.x$sample.index
  x = sample.index.x$x
  
  # skip duplicates: For phi=1, the true population
  #   is always Poisson.
  if (phi==1) {
    if (population != 'Negbin') {
      return(NA)
    }
  }

  parms <- paste(population,phi,n,b1,b2,N,sep=',')
  print(parms)
  
  p <- 2
  eta <- b1 + b2*x
  mu <- exp(eta)
  nu <- phi-1
  mn <- mu/nu
  w <-  (mn+1)/mn
  mx <- log(mu)-0.5*log(w)
  sx <- sqrt(log(w))

  
  
  ## vectors to hold CIs
  {
  ch.low <- vector()
  ch.upp <- vector()
  ch.phi.low <- vector()
  ch.phi.upp <- vector()
  
  fb.low <- vector()
  fb.upp <- vector()
  fb.phi.low <- vector()
  fb.phi.upp <- vector()
  
  hy.low <- vector()
  hy.upp <- vector()
  hy.phi.low <- vector()
  hy.phi.upp <- vector()
  
  vg.low <- vector()
  vg.upp <- vector()
  vg.phi.low <- vector()
  vg.phi.upp <- vector()
  }
  
  #  
  # Start simulations
  #
  for (sim in 1:N) {
    
    # show heartbeat
    if (sim %% 500 == 0) { print(sim) }
    
    y <- getY(phi,population,n,mu,mn,nu)
    
    # fit Poisson model
    model <- glm(y~x,family="poisson")
    muhat<-fitted(model)
    P<-sum((y-muhat)^2/muhat)
    sbar<-mean((y-muhat)/muhat) 
    phihat<-(P/(n-p))/(1+sbar) 
    
    ##
    df <- n-p
    summ <- summary(model)
    coef.x <- summ$coefficients[2,1]
    coef.se <- summ$coefficients[2,2]

    ## chisquare 
    ##
    intervals <- intervals.chisq(coef.x,coef.se,df,phihat) 
    ch.low[sim] <- intervals[1]
    ch.upp[sim] <- intervals[2]
    ch.phi.low[sim] <- intervals[3]
    ch.phi.upp[sim] <- intervals[4]
    
    ## Fastboot
    ##
    intervals <- intervals.fastboot(y=y,
                                    muhat=muhat,
                                    sample.index=sample.index,
                                    n=n,
                                    df=df,
                                    sbar=sbar,
                                    phihat=phihat,
                                    coef.x=coef.x,
                                    coef.se=coef.se) 
    fb.low[sim] <- intervals[1]
    fb.upp[sim] <- intervals[2]
    fb.phi.low[sim] <- intervals[3]
    fb.phi.upp[sim] <- intervals[4]
    
    
    intervals <- intervals.vgam(x=x,y=y) 
    vg.low[sim] <- intervals[1]
    vg.upp[sim] <- intervals[2]
    vg.phi.low[sim] <- intervals[3]
    vg.phi.upp[sim] <- intervals[4]
    
  }
  
  #print("---vg.low---")
  #print(vg.low)
  #print("------------------")
  ##
  ## end of simulations
  ##
  
  # calculate hybrid intervals
  h <- calcHybridInterval(ch.low,ch.upp,fb.low,fb.upp) 
  hy.low <- h$h.low
  hy.upp <- h$h.upp
  h.phi <- calcHybridInterval(ch.phi.low,ch.phi.upp,fb.phi.low,fb.phi.upp)
  hy.phi.low <- h.phi$h.low
  hy.phi.upp <- h.phi$h.upp
  
  # combine power calculations
  
  # beta
  ch.pow <- calcExcludesValue(val=0,low=ch.low,upp=ch.upp)
  fb.pow <- calcExcludesValue(val=0,low=fb.low, upp=fb.upp)
  vg.pow <- calcExcludesValue(val=0,low=vg.low, upp=vg.upp)
  hy.pow <- calcExcludesValue(val=0,low=hy.low, upp=hy.upp)
  
  # phi
  ch.phi.pow <- calcExcludesValue(val=1,low=ch.phi.low,upp=ch.phi.upp)
  fb.phi.pow <- calcExcludesValue(val=1,low=fb.phi.low,upp=fb.phi.upp)
  vg.phi.pow <- calcExcludesValue(val=1,low=vg.phi.low,upp=vg.phi.upp)
  hy.phi.pow <- calcExcludesValue(val=1,low=hy.phi.low,upp=hy.phi.upp)
  
  result.pow   <- c(ch.pow,
                    fb.pow,
                    vg.pow,
                    hy.pow,
                    ch.phi.pow,
                    fb.phi.pow,
                    vg.phi.pow,
                    hy.phi.pow)
  
  # combine coverage rates
  #
  ch.cov  <- calcIncludesValue(b2,ch.low,ch.upp)
  ch.cov.phi <- calcIncludesValue(phi,ch.phi.low,ch.phi.upp)
  
  fb.cov <- calcIncludesValue(b2,fb.low,fb.upp)
  fb.cov.phi <- calcIncludesValue(phi,fb.phi.low,fb.phi.upp)
  
  vg.cov  <- calcIncludesValue(b2,vg.low,vg.upp)
  vg.cov.phi  <- calcIncludesValue(phi,vg.phi.low,vg.phi.upp)
  
  hy.cov <- calcIncludesValue(b2,hy.low,hy.upp)
  hy.cov.phi <- calcIncludesValue(phi,hy.phi.low,hy.phi.upp)
  
  result.cov  <- c(ch.cov,
                   ch.cov.phi,
                   fb.cov,
                   fb.cov.phi,
                   vg.cov,
                   vg.cov.phi,
                   hy.cov,
                   hy.cov.phi)
  
  # combine CI median widths
  #
  ch.width <- calcMedianWidth(ch.low,ch.upp)
  ch.width.phi <- calcMedianWidth(ch.phi.low,ch.phi.upp)
  
  fb.width <- calcMedianWidth(fb.low,fb.upp)
  fb.width.phi <- calcMedianWidth(fb.phi.low,fb.phi.upp)
  
  result.wid <- c(ch.width,
                  ch.width.phi,
                  fb.width,
                  fb.width.phi,
                  vg.width,
                  vg.width.phi,
                  hy.width,
                  hy.width.phi)
  
  # gather the rates of NAs, ie. the rate of the method not working
  ch.na <- length(which(is.na(ch.low)))/N
  fb.na <- length(which(is.na(fb.low)))/N
  vg.na <- length(which(is.na(vg.low)))/N
  hy.na <- length(which(is.na(hy.low)))/N
  result.na <- c(ch.na,fb.na, vg.na,hy.na)

  # combine all results. This one will go to the csv file.
  result <- c(as.character(population),phi,n,b1,b2,N,result.cov,result.pow,result.wid,result.na)
  return(result)
}

###################################################
##
## end of doCalculations
##
###################################################

getColnames <- function() {
  return(c('population',
           'phi',
           'n',
           'b1',
           'b2',
           'N',
           'ch cov b','ch cov phi','fb cov b','fb cov phi', 'vg cov b','vg cov phi', 'hy cov b', 'hy cov phi',
           'ch pow b','ch pow phi','fb pow b','fb pow phi', 'vg pow b','vg pow phi', 'hy pow b', 'hy pow phi',
           'ch wid b','ch wid phi','fb wid b','fb wid phi', 'vg wid b','vg wid phi', 'hy wid b', 'hy wid phi',
           'ch err', 'fb.err','vg err','hy err'))
}


print('started mainpar14.R')
print(Sys.time())
set.seed(10)
N <- 2500
phis <- c(1,3,5) # 5
pops <- c("Negbin" ,"Neyman","Poisson lognormal")
beta1 <- c(-3,0.1,2.3)
ns <-  c(30,100,1000)
combs <- expand.grid(phis,pops,beta1,ns)
colnames(combs) <- c('phi','pop','beta1','n')
combs$beta2 <- NA
combs[which(combs$beta1==-3),'beta2'] <- 3
combs[which(combs$beta1==0.1),'beta2'] <- 2.2
combs[which(combs$beta1==2.3),'beta2'] <-0.7

# we use a CPU cores for every combination


# All calculations, for writing to csv file. This will be
#   vector of strings: Each line is the comma separated results for
#   one set of parameters (population,phi,n,beta)


i = 3   
calcs <- doCalculations(population = combs$pop[i],
                  phi= combs$phi[i],
                  n = combs$n[i],
                  b1 = combs$beta1[i],
                  b2 = combs$beta2[i],
                  N = 500 )
dim(calcs)
#dropm <- which(is.na(calcs[,1]))
#calcs <- calcs[-dropm,]
#dim(calcs)
#write.csv(calcs,'results-mainpar14.csv',row.names=F,quote=F)
#print('finished')
#print(Sys.time())




