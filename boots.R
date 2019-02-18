library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"

# calculate CI's for phi, based on chisquare distribution
#   with phihat1 and phihat2, and on Gamma distribution
#   (phihat2 only), and fitting only Negative binomial 
#   for now.
#   Plots first 100 of the simulations, with error rates.
#

phi <- 5
  
  
  {
    n<-30
    
    beta<-c(1,1)
    x<-seq(0,1,length.out=n)
    p<-length(beta)
    eta<-beta[1]+beta[2]*x
    mu<-exp(eta)
    nu<-phi-1
    mn<-mu/nu
    w<-(mn+1)/mn
    mx<-log(mu)-0.5*log(w)
    sx<-sqrt(log(w))
    
    N <- 500
    # bounds of CI for chisq and phi1
    
    lower.1 <- vector()
    upper.1 <- vector()
    # bounds of CI for chisq and phi2
    lower.2 <- vector()
    upper.2 <- vector()
    # bounds of CI for gamma and phi2
    gamma.lower <- vector()
    gamma.upper <- vector()
    
    boots2.low <- vector()
    boots2.upp <- vector()
    
    # estimates of phi
    phihats1 <- vector()
    phihats2 <- vector()
    phi.hats.boots2 <- vector()     
    # does CI contain true phi?
    err.1 <- vector()
    err.2 <- vector()
    gamma.err <- vector()
    boots.err <- vector()
    boots2.err <- vector()
    ks <- vector()
    thetas <- vector()
  }
  
  for (sim in 1:N) { 
    print(sim)
    if (phi==1) { 
      y<-rpois(n,mu) } 
    else { 
      y<-rnegbin(n,mu,mn) # negative binomial
      
      #y<-rpois(n,rpois(n,mn)*nu)             # Neyman Type A
      #y<-rpois(n,rlnorm(n,mx,sx))            # Poisson lognormal
      #y<-rpois(n,rinvgauss(n,1,1/mn)*mu)     # Poisson inverse Gaussian
      
    }  
    muhat<-fitted(glm(y~x,family="poisson"))
    P<-sum((y-muhat)^2/muhat)
    sbar<-mean((y-muhat)/muhat) 
    phihat1<-P/(n-p)
    phihat2<-phihat1/(1+sbar) 
    phihats1[sim] <- phihat1
    phihats2[sim] <- phihat2
    
    # chisquare CI
    df <- n-p
    lower.1[sim] <- df * phihat1  / qchisq(0.975,df=df)
    upper.1[sim] <- df * phihat1  / qchisq(0.025,df=df)
    lower.2[sim] <- df * phihat2  / qchisq(0.975,df=df)
    upper.2[sim] <- df * phihat2  / qchisq(0.025,df=df)
    
    err.1[sim] <- 0
    if (phi < lower.1[sim] | phi > upper.1[sim]) {
      err.1[sim] <-  1
    }
    err.2[sim] <- 0
    if (phi < lower.2[sim] | phi > upper.2[sim]) {
      err.2[sim] <-  1
    }
    #
    e <- y-muhat
    alpha3.hat <- (1/df) * sum( (e^3/muhat))
    alpha4.hat <- (1/df) * sum((e^4 - 3 * e^2)/muhat)
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
    
    ktheta <- df - (alpha3.hat/phihat2-phihat2)*S/n
    ktheta2 <- 2 * ktheta #df * (2 + tau) 
    
    theta <- ktheta2 / ktheta
    k <- ktheta / theta
    
    
    ## Gamma CI. only for phihat2
    #e <- (y-muhat)/sqrt(phihat2 * muhat)
    #rho3 <- e^3
    #rho4 <- e^4 - 3
    #alpha.bar <- mean(rho4 + phihat2/muhat - 2 * sqrt(phihat2/muhat) * rho3)
    #theta <- 2 + alpha.bar
    #k <- df/(2+alpha.bar)
    thetas[sim] <- theta
    ks[sim] <- k
    gamma.lower[sim] <- df * phihat2  / qgamma(0.975,
                                             shape=k,
                                             scale=theta)
    gamma.upper[sim] <- df * phihat2  / qgamma(0.025,
                                             shape=k,
                                             scale=theta)
    gamma.err[sim] <- (phi < gamma.lower[sim]) |
                 (phi > gamma.upper[sim])
    
    ## bootstrap
    phihat.boots <- vector()
    for (i in 1:0) {
      subsample.index <- sample(seq(1,n),replace = T)
      y.s <- y[subsample.index]
      x.s <- x[subsample.index]
      m <- glm(y.s~x.s,family="poisson")
      summary(m)
      muhat<-fitted(m)
      P<-sum((y-muhat)^2/muhat)
      sbar<-mean((y-muhat)/muhat) 
      phihat1.boots<-P/(n-p)
      phihat2.boots<-phihat1.boots/(1+sbar) 
      phihat.boots[i] <- phihat2.boots
    }
    boots.int <- quantile(phihat.boots,c(0.025,0.975))
    boots.err[sim] <- (phi < boots.int[1]) |
      (phi > boots.int[2])
  

  # boots2
  B = 1
  phihat.boots2 <- vector()
  for (i in 1:B) {
    subsample.index <- sample(seq(1,n),replace = T)
    e. <- (y-muhat)^2
    e <- e.[subsample.index]/muhat[subsample.index]
    
    P<-sum(e)
    phi.1.b <- P/df
    phi.2.b <- phi.1.b / (1+sbar)
    
    
    sbar<-mean((y-muhat)/muhat) 
    phihat1<-P/(n-p)
    phihat2<-phihat1/(1+sbar)
    phihat.boots2[i] <- phihat2
  }
  
  boots2.int <- quantile(phihat.boots2,c(0.025,0.975))
  phi.hats.boots2[sim] <- mean(phihat.boots2)
  boots2.low[sim] <- boots2.int[1]
  boots2.upp[sim] <- boots2.int[2]
  
  boots2.err[sim] <- (phi < boots2.int[1]) |
    (phi > boots2.int[2])
  
  
  }
  
  # error rates
  e1 <- sum(err.1)/N
  e2 <- sum(err.2)/N
  e3 <- sum(gamma.err,na.rm=T)/length(which(!is.na(gamma.err)))
  paste(e1,e2,e3,sep='   ')
  sum(boots.err)/N
  e5 <- sum(boots2.err)/N
  
  
  
  # plot CIs, show first plotN simulations
  plotCIs <- function(low,upp,err,phihat,truePhi,plotN,main,
                      errRate) {
    plot('',xlim=c(0,plotN),ylim=c(0,4*truePhi),
         main=paste(main,'true phi=',truePhi,sep=' '),
         ylab="95% CI")
    
    abline(h=truePhi,col='red')
    for (i in 1:plotN) {
      col <- 'black'
      if (err[i] == 1) {
        col <- 'red'
      }
      segments(i,low[i],i,upp[i],col=col)
      points(i,phihat[i],pch='.',cex=3)
      text(60,3.5 * truePhi,paste("error rate: ",errRate))
    }  
  }  
  
  par(mfrow=c(1,3))
  #plotCIs(lower.1,upper.1,err.1,phihats1,phi,plotN=100,
  #        main="Chisq, phihat1",e1)
  plotCIs(lower.2,upper.2,err.2,phihats2,phi,plotN=100,
          main="Chisq, phihat2",e2)
  plotCIs(gamma.lower,gamma.upper,gamma.err,phihats2,phi,plotN=100,
          main="Gamma, phihat2",e3)
  plotCIs(boots2.low,boots2.upp,boots2.err,phi.hats.boots2,phi,plotN=100,
          main="bootstrap2",e5)
  
  
  par(mfrow=c(3,1))
  
  median(upper.1-lower.1)
  median(upper.2-lower.2)
  median(gamma.upper-gamma.lower,na.rm=T)
  
  # plot estimates of phihat's
  
  
  #par(mfrow=c(2,2))
  {
    print(paste(phi,
                round(mean((n-p)*phihats1/phi),2),
                round(var((n-p)*phihats1/phi),2),
                round(mean((n-p)*phihats2/phi),2),
                round(var((n-p)*phihats2/phi),2),
                sep='   '))
    
    par(mfrow=c(1,1))
    main=paste('phi=',phi)
    hist((n-p)*phihats2/phi,probability = T,
         main=main,ylim=c(0,0.4),
         xlim=c(0,2 * n),100)
    c<- seq(0,qchisq(0.9999,df),0.1)
    y1<-dchisq(c,df=n-p)
    lines(c,y1,col='red') 
    
    y2<- dgamma(c,shape=mean(ks),scale=mean(thetas))
    lines(c,y2,col='blue') 
    
  
  
}

#mean(phihats1)
df * mean(phihats2)/ phi
2 * mean(ks)
df

hist(phihats2)
abline(v=phi,col='red')
abline(v=mean(phihats2,col='blue'))
#par(mfrow=c(2,2))

##(n / df) * phi
#var(phihats1)

# these should be roughly equal:
#(1/(df)) * (sum(y^2/muhat)-sum(muhat))
#phihat2
#


