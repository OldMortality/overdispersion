library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"
library(VGAM)
#
library(profvis)   
library(foreach)
library(doParallel)

# Program to simulate CI's for the beta's/
#
# 
# Calculate CI's for phi, based on 6 methods:
#   Chisquare distribution, based on phihat
#   Bootstrap, based on phihat.
#   VGAM restricted regression
#   
#   
#   
#
#
# Change log
#   04/11/2019: Created this program.
#

# for testing only:##########
population = "Negbin"
phi = 1
n=100
b = c(-3,3)
## number of simulations
N = 100
#############################


# global constants
CONF_LEVEL <- 0.9
ALPHA_LOWER <- 0.05
ALPHA_UPPER <- 0.95
# number of bootstrap samples
N.BOOTS = 10000


#
# Returns the proportion of elements in CI,
#   which do not include zero. Elements containing
#   NA are disregarded. Returns NaN if all elements are NA,
#   but we hold that for impossible.
calcExcludesZero <- function(low,upp) {
  return(sum(low > 0 | upp < 0,na.rm = T)/length(which(!is.na(low))))
}

# returns median width of CIs
calcMedianWidth <- function(lower,upper) {
  return(median(upper-lower,na.rm = T))
}



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
#profvis({
#  t <- doCalculations(population,phi,n,b,N)

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
#   a confidence interval for phi, so we have:
#   3 coverage rates. 
#   3 proportions with phi.lower > 1
#   3 median widths of the CI for phi
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
  boots2.err <- vector() 
  vgam.err <- vector()	    
  
  # estimates of phi
  
  ## vectors to hold CIs
  {
  chisq.low <- vector()
  chisq.upp <- vector()
  boots2.low <- vector()
  boots2.upp <- vector()
  vgam.low <- vector()
  vgam.upp <- vector()
  
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
    
    chisq.err[sim] <- 0
    if (beta[2] < chisq.low[sim] | beta[2] > chisq.upp[sim]) {
      chisq.err[sim] <-  1
    }
    
    
    ##
    ## Bootstrap
    ##
    do.bootstrap <- T
    
    if (do.bootstrap == T) {
      
      
      #phihat.boots <- vector()
      e. <- (y-muhat)/sqrt(muhat)
      e.squared <- e.^2
      # split this into a matrix with 1 row for each bootstrap
      ss <- matrix(
        data=e.squared[sample.index],byrow=T,nrow=N.BOOTS,ncol=n)
      
      # 1 row for each bootstrap 
      
      phihat.boots <- (rowSums(ss)/df)/(1+sbar)
      ts <- rnorm(N.BOOTS) * sqrt(phihat/phihat.boots)
      #hist(ts) 
      t.low <- quantile(ts,ALPHA_LOWER)
      t.upp <- quantile(ts,ALPHA_UPPER)
      
      # phihat
      
      boots2.low[sim] <- coef.x - t.upp * phihat * coef.se
      boots2.upp[sim] <- coef.x - t.low * phihat * coef.se
      boots2.err[sim] <- (beta[2] < boots2.low[sim]) |
        (beta[2] > boots2.upp[sim])
    }	  
    
    ##
    ## VGAM
    ##
    
    do.RR <- T

    if (do.RR) {
    
      ci <- tryCatch(
        {
          a3 <- vglm(y~x,family=negbinomial(parallel=T,zero=NULL))
          # method = "wald
          # method = "profile"
          ci <- confint(a3,level = CONF_LEVEL)[3,]
          # estimate for phi
          
        },
        error=function(cond) {
          message(paste('rr',y,x,cond,sep=' '))
          return(NA)
        } ,
        warning=function(cond) {
          #ci <- confint(a3)[3,]
          return(NA)
        }
      
      )
      vgam.low[sim] <- ci[1]
      vgam.upp[sim] <- ci[2]
      
      # this will be NA if vgam.upp[sim]==NA
      vgam.err[sim] <- (beta[2] < ci[1]) | (beta[2] > ci[2])

    }
    
   
  }
  
  ##
  ## end of simulations
  ##
  
  # combine power calculations
  #
  chi.p1 <- calcExcludesZero(chisq.low,chisq.upp)
  bo2.p1 <- calcExcludesZero(boots2.low,boots2.upp)
  vga.p1 <- calcExcludesZero(vgam.low,vgam.upp)
  result.prop1 <- c(chi.p1,bo2.p1,vga.p1)
  
  # combine coverage rates
  #
  errorChisq <- sum(chisq.err)/N
  errorboots2 <- sum(boots2.err,na.rm=T)/length(which(!is.na(boots2.err)))
  errorVGAM  <- sum(vgam.err,na.rm=T)/length(which(!is.na(vgam.err)))
  # convert error rates to coverage
  result.coverage <- c(1-errorChisq,1-errorboots2,1-errorVGAM)
  
  # combine CI median widths
  #
  chi.width <- calcMedianWidth(chisq.low,chisq.upp)
  bo2.width <- calcMedianWidth(boots2.low,boots2.upp)
  vga.width <- calcMedianWidth(vgam.low,vgam.upp)
  result.width <- c(chi.width,bo2.width,vga.width)
  
  # gather the rates of NAs, ie. the rate of the method not working
  chi.na <- 0 # always works
  bo2.na <- length(which(is.na(boots2.err)))/N
  vga.na <- length(which(is.na(vgam.err)))/N
  result.na <- c(chi.na, bo2.na,vga.na)
  
  
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


print('started')
print(Sys.time())
pops <- c("Negbin","Neyman","Poisson lognormal")
betas <- rbind(c(-3,3),c(0.1,2.2),c(2.3,0.7))
ns <-  c(30,100,1000)
N <- 100 #00

# All calculations, for writing to csv file. This will be
#   vector of strings: Each line is the comma separated results for
#   one set of parameters (population,phi,n,beta)
results <- vector()

counter <- 0
phis <- c(1,3,5)
for (pop in pops) {
    for (n in ns) {
      for (coef in seq(1,3)) {
        
        print(paste(Sys.time(),'counter',counter,sep=': '))
         calcs <- foreach(ph = phis, .verbose=TRUE, .combine = rbind, .packages=c('VGAM','MASS'))  %dopar% {
                   doCalculations(population = pop,
                                 phi= ph,
                                 n = n,
                                 b = betas[coef,],
                                 N = N) 
        }
       
        # 
        for (j in 1:length(phis)) {
          counter <- counter + 1
          calcs.phi <- paste(pop,phis[j],n,betas[b,1],betas[b,2],N,sep=',')
          for (k in 1:12) {
            calcs.phi <- paste(calcs.phi,calcs[j,k],sep=',')
          }
          
          #print(calcs.phi)
          # store in results, which becomes a csv file
          results[counter] <- calcs.phi
        }
      }
    }
  
}
write.csv(results,'~/Documents/overdispersion/results-mainpar7.csv',row.names=F,quote=F)
print('finished')
print(Sys.time())