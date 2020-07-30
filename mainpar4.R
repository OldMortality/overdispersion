library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"
library(VGAM)
#library(profvis)   
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
N = 10000
#############################


# global constants
CONF_LEVEL <- 0.9
ALPHA_LOWER <- 0.05
ALPHA_UPPER <- 0.95
# number of bootstrap samples
N.BOOTS = 10000


# get the roots of the quadratic, for Estimating Equations
#   roots are returned in ascending order
getRoots <- function(xx,s,q95,n) {
  # the square term
  A <- sum(s)^2 + n^2 + 2 * n * sum(s) - q95*sum((s-mean(s))^2)
  # the linear term
  B <- -2 * sum(xx)*sum(s)-2*n*sum(xx)+
    2*q95*sum((xx-mean(xx))*(s-mean(s)))
  # the constant
  C <- (sum(xx)^2)-q95*sum((xx-mean(xx))^2)
  # the usual formula 
  discriminant <- sqrt(B^2-4*A*C)
  two.A <- 2 * A
  root1 <- (-B - discriminant)/two.A
  root2 <- (-B + discriminant)/two.A
  result <- c(root1,root2)
  # return in ascending order
  return(sort(result))
}

#
# Returns the proportion of elements in vector,
#   which are greater than 1. We will use this with
#   v as the lower limit of the CIs
calcPropGreaterOne <- function(v) {
  return(sum(v > 1,na.rm=T)/length(which(!is.na(v))))
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
#  t <- doCalculations(population,phi,n,b,500)
#})

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
    ## chisquare 
    ##
    df <- n-p
    chisq.low[sim] <- df * phihat2  / qchisq(ALPHA_UPPER,df=df)
    chisq.upp[sim] <- df * phihat2  / qchisq(ALPHA_LOWER,df=df)
    
    chisq.err[sim] <- 0
    if (phi < chisq.low[sim] | phi > chisq.upp[sim]) {
      chisq.err[sim] <-  1
    }
    
    ###
    ### Gamma
    ###
    
    e <- y-muhat
    alpha3.hat <- (1/df) * sum( (e^3/muhat))
    alpha4.hat <- (1/df) * sum(e^4/muhat) - 3 * muhat* phihat2^2
    tau.i <- (alpha4.hat/phihat2^2 - 2 * alpha3.hat/phihat2 + phihat2 )
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
    
    gamma.low[sim] <- df * phihat2  / qgamma(ALPHA_UPPER,
                                             shape=k,
                                             scale=theta)
    
    gamma.upp[sim] <- df * phihat2  / qgamma(ALPHA_LOWER,
                                               shape=k,
                                               scale=theta)
    gamma.err[sim] <- (phi < gamma.low[sim]) |
                 (phi > gamma.upp[sim])


    ##
    ## Bootstrap
    ##
    do.bootstrap <- T
    
    if (do.bootstrap == T) {
      
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
      ss <- matrix(data=e.squared[sample.index],byrow=T,nrow=N.BOOTS,ncol=n)
      
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
    }	  
    
    ##
    ## VGAM
    ##
    
    do.RR <- T

    if (do.RR) {
    
      ## RR regression
      ci <- tryCatch(
        {
          a3<-vglm(y~x,family=negbinomial(parallel=T,zero=NULL))
          # estimate for phi
          phihat3<-Confint.nb1(a3)$phi0
          # CI for phi
          phihat3ci<-Confint.nb1(a3,level=CONF_LEVEL)$CI.phi0
          c(phihat3ci[1],phihat3ci[2])
        },
        error=function(cond) {
          message(paste('rr',y,x,cond,sep=' '))
          return(NA)
        } ,
        warning=function(cond) {
          return(NA)
        }
      
      )
      vgam.low[sim] <- ci[1]
      vgam.upp[sim] <- ci[2]
      
      # this will be NA if vgam.upp[sim]==NA
      vgam.err[sim] <- (phi < ci[1]) | (phi > ci[2])

    }
    
    ##
    ## Estimating equations
    ## 
    do.EE <- T
    if (do.EE) {
    
      xx <- (y-muhat)^2/muhat
      s <- (y-muhat)/muhat
      z <- xx - phihat2*(1+s)
      
              # code before vectorisation
              # z.bar <- mean(z)
              # num <- (sum(z))^2
              # denom <- sum((z-z.bar)^2)
              # sumz[sim] <- sum(z)^2
              # chis[sim] <- num/denom
              # q.stars.old <- vector()
              # for (bt in 1:N.BOOTS) {
              #   boots.index <- sample.index[seq((bt-1)*n+1,bt*n)]
              #  z.star <- z[boots.index]
              #  z.star.bar <- mean(z.star)
              #  q.stars.old[bt] <- (sum(z.star)^2)/sum((z.star-z.star.bar)^2)
              # }
              # qConf <- quantile(q.stars.old,CONF_LEVEL,na.rm=T)
              # 
      ss <- matrix(data=z[sample.index],byrow=T,nrow=N.BOOTS,ncol=n)
      ss.mean.by.row <- apply(ss,1,mean)
      z.minus.zbar <- sweep(ss,1,ss.mean.by.row)
      # square(sum)/sum(square(z-z.bar))
      q.stars <- (rowSums(ss)^2)/ rowSums(z.minus.zbar^2)
      qConf <- quantile(q.stars,CONF_LEVEL,na.rm=T)
    
    
      # solve equation to get lower and upper limits of CI
    
      roots <- getRoots(xx,s,qConf,n)
    
      ee.low[sim] <- roots[1]
      ee.upp[sim] <- roots[2]
      ee.err[sim] <- ( ee.low[sim] > phi ) |
        ( ee.upp[sim] < phi )
    } # if do.EE
  }
  
  ##
  ## end of simulations
  ##
  
  # combine power calculations
  #
  chi.p1 <- calcPropGreaterOne(chisq.low)
  gam.p1 <- calcPropGreaterOne(gamma.low)
  bo1.p1 <- calcPropGreaterOne(boots1.low)
  bo2.p1 <- calcPropGreaterOne(boots2.low)
  vga.p1 <- calcPropGreaterOne(vgam.low)
  ee.p1  <- calcPropGreaterOne(ee.low)
  result.prop1 <- c(chi.p1,gam.p1,bo1.p1,bo2.p1,vga.p1,ee.p1)
  
  # combine coverage rates
  #
  errorChisq <- sum(chisq.err)/N
  errorGamma <- sum(gamma.err,na.rm=T)/length(which(!is.na(gamma.err)))
  errorboots1 <- sum(boots1.err,na.rm=T)/length(which(!is.na(boots1.err)))
  errorboots2 <- sum(boots2.err,na.rm=T)/length(which(!is.na(boots2.err)))
  errorVGAM  <- sum(vgam.err,na.rm=T)/length(which(!is.na(vgam.err)))
  errorEE    <- sum(ee.err,na.rm=T)/length(which(!is.na(ee.err)))
  # convert error rates to coverage
  result.coverage <- c(1-errorChisq,1-errorGamma,1-errorboots1,1-errorboots2,1-errorVGAM,1-errorEE)
  
  # combine CI median widths
  #
  chi.width <- calcMedianWidth(chisq.low,chisq.upp)
  gam.width <- calcMedianWidth(gamma.low,gamma.upp)
  bo1.width <- calcMedianWidth(boots1.low,boots1.upp)
  bo2.width <- calcMedianWidth(boots2.low,boots2.upp)
  vga.width <- calcMedianWidth(vgam.low,vgam.upp)
  ee.width  <- calcMedianWidth(ee.low,ee.upp)
  result.width <- c(chi.width,gam.width,bo1.width,bo2.width,vga.width,ee.width)
  
  # gather the rates of NAs, ie. the rate of the method not working
  chi.na <- 0 # always works
  gam.na <- length(which(is.na(gamma.err)))/N
  bo1.na <- length(which(is.na(boots1.err)))/N
  bo2.na <- length(which(is.na(boots2.err)))/N
  vga.na <- length(which(is.na(vgam.err)))/N
  ee.na  <- length(which(is.na(ee.err)))/N
  result.na <- c(chi.na,gam.na,bo1.na,bo2.na,vga.na,ee.na)
  
  # combine all results. This one will go to the csv file.
  result <- c(result.coverage,result.prop1,result.width,result.na)
  return(result)
}

###################################################
##
## end of doCalculations
##
###################################################



print('started')
print(Sys.time())
pops <- c("Negbin","Neyman","Poisson lognormal")
betas <- rbind(c(-3,3),c(0.1,2.2),c(2.3,0.7))
ns <-  c(30,100,1000)
N <- 10000

# All calculations, for writing to csv file. This will be
#   vector of strings: Each line is the comma separated results for
#   one set of parameters (population,phi,n,beta)
results<- vector()

counter <- 0
phis <- c(1,2,3,5)
for (pop in pops) {
    for (n in ns) {
      for (b in seq(1,3)) {
        
        print(paste(Sys.time(),'counter',counter,sep=': '))
         calcs <- foreach(ph = phis, .combine = rbind, .packages=c('VGAM','MASS'))  %dopar% {
                   doCalculations(population=pop,
                                 phi=ph,
                                 n=n,
                                 b=betas[b,],
                                 N=N) 
        }
       
        # 
        for (j in 1:length(phis)) {
          counter <- counter + 1
          calcs.phi <- paste(pop,phis[j],n,betas[b,1],betas[b,2],N,sep=',')
          for (k in 1:24) {
            calcs.phi <- paste(calcs.phi,calcs[j,k],sep=',')
          }
          
          #print(calcs.phi)
          # store in results, which becomes a csv file
          results[counter] <- calcs.phi
        }
      }
    }
  
}
write.csv(results,'~/Documents/overdispersion/results-mainpar5.csv',row.names=F,quote=F)
print('finished')
print(Sys.time())