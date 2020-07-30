library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"
library(VGAM)
#library(profvis)   
library(foreach)
library(doParallel)

# Hacked version of mainpar3.R, used to show CIs.
#
# calculate CI's for phi, based on chisquare distribution
#   with phihat1 and phihat2, and on Gamma distribution
#   (phihat2 only), and fitting only Negative binomial 
#   for now.
#
# parallel version
#   sample.index out of the inner loop
# change log
#   21 July: confidence level now 90% (was 95%)
#          : Added interval for bootstrap of phi1hat (the one not divided by 1+s.bar)
#          : Vectorised bootstrap (also at EE)
#





CONF_LEVEL <- 0.9
ALPHA_LOWER <- 0.05
ALPHA_UPPER <- 0.95


# get the roots of the quadratic, for Estimating Equations
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

SEQ.30 <- seq(1,30)
SEQ.100 <- seq(1,100)
SEQ.1000 <- seq(1,1000)
no_cores= 4
registerDoParallel(makeCluster(no_cores))
##
## 
##
X.30 <- seq(0,1,length.out=30)
X.100 <- seq(0,1,length.out=100)
X.1000 <- seq(0,1,length.out=1000)


# number of bootstrap samples
N.BOOTS = 100


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

sample.index.30 <- matrix(create.samples(30),nrow=1,ncol=30*N.BOOTS)
sample.index.100 <- matrix(create.samples(100),nrow=1,ncol=100*N.BOOTS)
sample.index.1000 <- matrix(create.samples(1000),nrow=1,ncol=1000*N.BOOTS)

chi.low <- vector()
chi.upp <- vector()
chi.hat <- vector()
gam.low <- vector()
gam.upp <- vector()
gam.hat <- vector()
bo2.low <- vector()
bo2.upp <- vector()
bo2.hat <- vector()


# for testing only:
population = "Negbin"
phi = 5
n=100
b = c(-3,3)
## number of simulations
N = 100
calcCoverage(population,phi,n,b,500)
#profvis({
#  calcCoverage(population,phi,n,b,500)
#})


calcCoverage <- function(population,phi,n,b,N) {
  
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
  #x<-seq(0,1,length.out=n)
  p<-length(beta)
  eta<-beta[1]+beta[2]*x
  mu<-exp(eta)
  nu <- phi-1
  mn<-mu/nu
  w<-(mn+1)/mn
  mx<-log(mu)-0.5*log(w)
  sx<-sqrt(log(w))

  # errors of CIs
  chisq.err <- vector()
  gamma.err <- vector()
  boots1.err <- vector() # based on phi1.hat
  boots2.err <- vector() # based on phi2.hat (phi1.hat / (1+s.bar))
  vgam.err <- vector()	    
  ee.err <- vector()
  

  # estimates of phi
  #phihats2 <- vector()
  phi.hats.boots <- vector()
  # used by EE
  chis <- vector()
  sumz <- vector() 

  #  
  # start simulations
  #
  for (sim in 1:N) {
    
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
    chisq.low <- df * phihat2  / qchisq(ALPHA_UPPER,df=df)
    chisq.upp <- df * phihat2  / qchisq(ALPHA_LOWER,df=df)
    
    chi.low[sim] <<- chisq.low
    chi.upp[sim] <<- chisq.upp
    chi.hat[sim] <<- phihat2
    
    chisq.err[sim] <- 0
    if (phi < chisq.low | phi > chisq.upp) {
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
    
    gamma.low <- df * phihat2  / qgamma(ALPHA_UPPER,
                                             shape=k,
                                             scale=theta)
    
    gamma.upp <- df * phihat2  / qgamma(ALPHA_LOWER,
                                               shape=k,
                                               scale=theta)
    gamma.err[sim] <- (phi < gamma.low) |
                 (phi > gamma.upp)

    gam.low[sim] <<- gamma.low
    gam.upp[sim] <<- gamma.upp
    gam.hat[sim] <<- phihat2

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
      boots1.low <- (df * phihat1 ) / boots1.int[2]
      boots1.upp <- (df * phihat1 ) / boots1.int[1]
      boots1.err[sim] <- (phi < boots1.low) |
        (phi > boots1.upp)
      # phihat2
      boots2.int <- (df /phihat2) *quantile(phihat2.boots,c(ALPHA_LOWER,ALPHA_UPPER)) 
      boots2.low <- (df * phihat2 ) / boots2.int[2]
      boots2.upp <- (df * phihat2 ) / boots2.int[1]
      boots2.err[sim] <- (phi < boots2.low) |
        (phi > boots2.upp)
      
      bo2.low[sim] <<- boots2.low
      bo2.upp[sim] <<- boots2.upp
      bo2.hat[sim] <<- phihat2.boots
      print(bo2.low[sim])
    }	  
    
    ##
    ## VGAM
    ##
    
    do.RR <- F

    if (do.RR) {
    
      ## RR regression
      tryCatch(
        {
          a3<-vglm(y~x,family=negbinomial(parallel=T,zero=NULL))
          # estimate for phi
          phihat3<-Confint.nb1(a3)$phi0
          # CI for phi
          phihat3ci<-Confint.nb1(a3,level=CONF_LEVEL)$CI.phi0
          vgam.err[sim] <- (phi < phihat3ci[1]) |
            (phi > phihat3ci[2])
        
        },
        error=function(cond) {
          message(paste('rr',y,x,cond,sep=' '))
          return(NA)
        } 
      
      )

    }
    
    ##
    ## Estimating equations
    ## 
    do.EE <- F
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
    
      ee.err[sim] <- ( roots[1] > phi ) |
        ( roots[2] < phi )
    } # if do.EE
  }
  
  ##
  ## end of simulations
  ##
  
  # combine error rates
  
  errorChisq <- sum(chisq.err)/N
  errorGamma <- sum(gamma.err,na.rm=T)/length(which(!is.na(gamma.err)))
  errorboots1 <- sum(boots1.err,na.rm=T)/length(which(!is.na(boots1.err)))
  errorboots2 <- sum(boots2.err,na.rm=T)/length(which(!is.na(boots2.err)))
  errorVGAM  <- sum(vgam.err,na.rm=T)/length(which(!is.na(vgam.err)))
  errorEE    <- sum(ee.err,na.rm=T)/length(which(!is.na(ee.err)))
  # convert error rates to coverage
  result <- c(1-errorChisq,1-errorGamma,1-errorboots1,1-errorboots2,1-errorVGAM,1-errorEE)
  return(result)
}


plotCI <-function(l,u,h,phi,plotN=100,main) {
  sims <- length(l)
  plot('',xlim=c(1,plotN),ylim = c(0,20),main=main)
  
  for (i in 1:plotN) {
    if (is.na(l[i]) | is.na(u[i] | is.na(h[i]))) {
      # do nothing
      next
    } 
    c = 'black'
    if (phi<l[i] | phi>u[i]) {
      #if (i %% 2 == 0) {
      c = 'red'
    } else {
      c = 'black'
    }
    segments(i,l[i],i,u[i],col=c)
    #points(i,h[i],col=c)
  }
  abline(h=phi,col='blue')
}

par(mfrow=c(3,1))
plotCI(chi.low,chi.upp,chi.hat,phi,plotN=100,main='Chisq')
plotCI(gam.low,gam.upp,gam.hat,phi,plotN=100,main='Gamma')
plotCI(bo2.low,bo2.upp,bo2.hat,phi,plotN=100,main='B2')


betas <- rbind(c(-3,3),c(0.1,2.2),c(2.3,0.7))
x <- c(0,1)
betas[1,] %x% x
