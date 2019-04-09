library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"
library(VGAM)
   

# calculate CI's for phi, based on chisquare distribution
#   with phihat1 and phihat2, and on Gamma distribution
#   (phihat2 only), and fitting only Negative binomial 
#   for now.
#   Plots first 100 of the simulations, with error rates.
#



# for testing:
population = "Negbin"
phi = 2
n=30
b = c(-3,3)
N = 10000
#calcCoverage(population,phi,n,b,N)

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
  # bounds of CI for chisq
  lower.chisq <- vector()
  upper.chisq <- vector()
  err.chisq <- vector()
    
  # bounds of CI for gamma 
  gamma.lower <- vector()
  gamma.upper <- vector()
  gamma.err <- vector()
    
  # bounds of CI for bootstrap
  boots2.low <- vector()
  boots2.upp <- vector()  
  boots2.err <- vector()
    
  # bounds of CI for vgam
  vgam.low <- vector()
  vgam.upp <- vector()
  vgam.err <- vector()	    

  b2.hat.low.chisq <-  vector()
  b2.hat.upp.chisq <-  vector()
  
  b2.hat.low.gamma <-  vector()
  b2.hat.upp.gamma <-  vector()

  err.b2.chisq <- vector()
  err.b2.gamma <- vector()
  
  # estimates of phi
  phihats2 <- vector()
  phi.hats.boots2 <- vector()


  # does CI contain true phi?
  ks <- vector()
  
  for (sim in 1:N) {
    if ( sim %% 100 == 0) { print(sim) }
    # set true population
    {
    if (phi==1) { 
      y <- rpois(n,mu) 
    } else 
      { 
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
    m <- glm(y~x,family="poisson")
    muhat<-fitted(m)
    P<-sum((y-muhat)^2/muhat)
    sbar<-mean((y-muhat)/muhat) 
    phihat1<-P/(n-p)
    phihat2<-phihat1/(1+sbar) 
    phihats2[sim] <- phihat2
    
    # chisquare CI
    df <- n-p
    lower.chisq[sim] <- df * phihat2  / qchisq(0.95,df=df)
    upper.chisq[sim] <- df * phihat2  / qchisq(0.05,df=df)
    
    err.chisq[sim] <- 0
    if (phi < lower.chisq[sim] | phi > upper.chisq[sim]) {
      err.chisq[sim] <-  1
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

    N2 <- 10000
    num <- rnorm(N)
    denom1 <- sqrt(rchisq(N2,df=df)/df)
    denom2 <- sqrt(rgamma(N2,shape=k,scale=theta)/df)
    t1 <- num/denom1
    t2 <- num/denom2
    #hist(t,30)
    qt1 <- quantile(t1,0.975)
    
    qt2 <- NA
    if (k>0 && theta>0) {
      qt2 <- quantile(t2,0.975)
    }
    
    b2.hat.low.chisq[sim] <- summary(m)$coefficients[2,1] -
              qt1 * summary(m)$coefficients[2,2]                           
    b2.hat.upp.chisq[sim] <-  summary(m)$coefficients[2,1] +
              qt1 * summary(m)$coefficients[2,2]    
    err.b2.chisq[sim] <- (beta[2] < b2.hat.low.chisq[sim]) |
      (beta[2] > b2.hat.upp.chisq[sim])
    
    
    b2.hat.low.gamma[sim] <- summary(m)$coefficients[2,1] -
      qt2 * summary(m)$coefficients[2,2]                           
    b2.hat.upp.gamma[sim] <-  summary(m)$coefficients[2,1] +
      qt2 * summary(m)$coefficients[2,2]    
    err.b2.gamma[sim] <- (beta[2] < b2.hat.low.gamma[sim]) |
      (beta[2] > b2.hat.upp.gamma[sim])
    
                               
              
              
                                         
                                         
    
    

    # bootstrap

    do.bootstrap <- F
    
    if (do.bootstrap == T) {

      
      B = 10000
      phihat.boots2 <- vector()
    
      e. <- (y-muhat)/sqrt(muhat)
      for (i in 1:B) {
        #print(i)
        subsample.index <- sample(seq(1,n),replace = T)
        e <- e.[subsample.index]
    
        P.b<-sum(e^2)
        phi.1.b <- P.b/df
        phi.2.b <- phi.1.b / (1+sbar)
        phihat.boots2[i] <- phi.2.b
      }
  
      boots2.int <- (df /phihat2) *quantile(phihat.boots2,c(0.025,0.975)) 
      #phi.hats.boots2[sim] <- mean(phihat.boots2)
      boots2.low[sim] <- (df * phihat2 ) / boots2.int[2]
      boots2.upp[sim] <- (df * phihat2 ) / boots2.int[1]
  
      boots2.err[sim] <- (phi < boots2.low[sim]) |
        (phi > boots2.upp[sim])
    }	  
    
    ## RR regression
    tryCatch(
      {
        a3<-vglm(y~x,family=negbinomial(parallel=T,zero=NULL))
        phihat3<-Confint.nb1(a3)$phi0
        phihat3ci<-Confint.nb1(a3)$CI.phi0
        vgam.low[sim] <- phihat3ci[1] 
        vgam.upp[sim] <- phihat3ci[2] 
        vgam.err[sim] <- (phi < vgam.low[sim]) |
          (phi > vgam.upp[sim])
        
      },
      error=function(cond) {
        message(cond)
        return(NA)
      }
      
    )
    


  }
  # error rates
  errorChisq <- sum(err.chisq)/N
  errorGamma <- sum(gamma.err,na.rm=T)/length(which(!is.na(gamma.err)))
  errorBoots2 <- sum(boots2.err,na.rm=T)/length(which(!is.na(boots2.err)))
  errorVGAM <- sum(vgam.err,na.rm=T)/length(which(!is.na(vgam.err)))
  error.b2.chisq <- sum(err.b2.chisq)/length(err.b2.chisq)
  error.b2.gamma <- sum(err.b2.gamma,na.rm=T)/length(which(!is.na(err.b2.gamma)))
  
  result <- c(1-errorChisq,1-errorGamma,1-errorBoots2,1-errorVGAM,error.b2.chisq,
              error.b2.gamma)
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
        parms <- paste(pop,phi,n,betas[b,1],betas[b,2],N,coverage[1],coverage[2],coverage[3],coverage[4],sep=',')
        results[counter] <- parms
        print(parms)
        
      }
    }
  }
}
write.csv(results,'results.csv',row.names=F,quote=F)




