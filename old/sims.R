library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"

# calculate CI's for phi, based on chisquare distribution
#   with phihat1 and phihat2, and on Gamma distribution
#   (phihat2 only), and fitting only Negative binomial 
#   for now.
#   Plots first 100 of the simulations, with error rates.
#

phi <- 10
  
  set.seed(10)
  
  {
    n <- 100
    
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
    
    N <- 10000
    # bounds of CI for chisq and phi1
    
    # bounds of CI for chisq and phi2
    lower.2 <- vector()
    upper.2 <- vector()
    # bounds of CI for gamma and phi2
    gamma.lower <- vector()
    gamma.upper <- vector()
    
    # estimates of phi
    phihats2 <- vector()
    
    # does CI contain true phi?
    err.1 <- vector()
    err.2 <- vector()
    gamma.err <- vector()
    ks <- vector()
    thetas <- vector()
    kthetas <- vector()
    biases <- vector()
  }
  
  for (sim in 1:N) { 
    print(sim)
    if (phi==1) { 
      y<-rpois(n,mu) } 
    else { 
      #y<-rnegbin(n,mu,mn) # negative binomial
      
      #y<-rpois(n,rpois(n,mn)*nu)             # Neyman Type A
      #y<-rpois(n,rlnorm(n,mx,sx))            # Poisson lognormal
      y<-rpois(n,rinvgauss(n,1,1/mn)*mu)     # Poisson inverse Gaussian
      
    }  
    muhat<-fitted(glm(y~x,family="poisson"))
    P<-sum((y-muhat)^2/muhat)
    sbar<-mean((y-muhat)/muhat) 
    phihat1<-P/(n-p)
    phihat2<-phihat1/(1+sbar) 
    phihats2[sim] <- phihat2
    
    # chisquare CI
    df <- n-p
    lower.2[sim] <- df * phihat2  / qchisq(0.95,df=df)
    upper.2[sim] <- df * phihat2  / qchisq(0.05,df=df)
    
    err.2[sim] <- 0
    if (phi < lower.2[sim] | phi > upper.2[sim]) {
      err.2[sim] <-  1
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
    
    biases[sim] <- bias.phi
    
    # phihat. <- phihat2 - bias.phi
    # bias.phi <- -(alpha3.hat-phihat.^2)*S/n/df
    # phihat. <- phihat2 - bias.phi
    # bias.phi <- -(alpha3.hat-phihat.^2)*S/n/df
    # phihat. <- phihat2 - bias.phi
    # bias.phi <- -(alpha3.hat-phihat.^2)*S/n/df 
    
    
    # plus or minus?
    ktheta <- df + df * bias.phi / phihat2
    ktheta2 <- df * (2 + tau) 
    
    theta <- ktheta2 / ktheta
    k <- ktheta / theta
    
    kthetas[sim] <- ktheta
    
    
    thetas[sim] <- theta
    ks[sim] <- k 
    gamma.lower[sim] <- df * phihat2  / qgamma(0.95,
                                             shape=k,
                                             scale=theta)
    
    gamma.upper[sim] <- df * phihat2  / qgamma(0.05,
                                               shape=k,
                                               scale=theta)
                                               
    gamma.err[sim] <- (phi < gamma.lower[sim]) |
                 (phi > gamma.upper[sim])
    
  }
  
  # error rates
  #e1 <- sum(err.1)/N
  e2 <- sum(err.2)/N
  e3 <- sum(gamma.err,na.rm=T)/length(which(!is.na(gamma.err)))
  paste(e2,e3,sep='   ')
  
  # plot CIs, show first plotN simulations
  plotCIs <- function(low,upp,err,phihat,truePhi,plotN,main,
                      errRate) {
    plot('',xlim=c(0,plotN),ylim=c(0,4*truePhi),
         main=paste(main,'true phi=',truePhi,sep=' '),
         ylab="95% CI")
    
    abline(h=truePhi,col='red')
    for (i in 1:plotN) {
      col <- 'black'
      if (!is.na(err[i]) && err[i] == 1) {
        col <- 'red'
      }
      segments(i,low[i],i,upp[i],col=col)
      points(i,phihat[i],pch='.',cex=3)
      text(60,3.5 * truePhi,paste("error rate: ",errRate))
    }  
  }  
  
  par(mfrow=c(1,2))
  plotCIs(lower.2,upper.2,err.2,phihats2,phi,plotN=100,
          main="Chisq, phihat2",e2)
  plotCIs(gamma.lower,gamma.upper,gamma.err,phihats2,phi,plotN=100,
          main="Gamma, phihat2",e3)
  
  #par(mfrow=c(1,3))
  
  median(upper.2-lower.2)
  median(gamma.upper-gamma.lower,na.rm=T)
  
  # plot estimates of phihat's
  
  
  par(mfrow=c(2,1))
  
    print(paste(phi,
                round(mean((n-p)*phihats1/phi),2),
                round(var((n-p)*phihats1/phi),2),
                round(mean((n-p)*phihats2/phi),2),
                round(var((n-p)*phihats2/phi),2),
                sep='   '))
    
    
    main=paste('phi=',phi)
    hist((n-p)*phihats2/phi,probability = T,
         main=main,ylim=c(0,0.10),
         xlim=c(0,2 * n),100)
    c<- seq(0,qchisq(0.9999,df),0.1)
    y1<-dchisq(c,df=n-p)
    lines(c,y1,col='black',lwd=3)
    
    #par(mfrow=c(1,1))
    main=paste('phi=',phi)
    hist((n-p)*phihats2/phi,probability = T,
         main=main,ylim=c(0,0.10),
         xlim=c(0,2 * n),100)
    c<- seq(0,qchisq(0.9999,df),0.1)
    y1<-dchisq(c,df=n-p)
    lines(c,y1,col='black',lwd=3) 
    
    gs <- 10
    for (i in 1:gs) {    
      y2<- dgamma(c,shape=ks[i],scale=2)
      lines(c,y2,col=rainbow(gs)[i]) 
    }
    y1<-dchisq(c,df=n-p)
    lines(c,y1,col='black',lwd=3) 
    


#mean(phihats1)
df * mean(phihats2)/ phi
2 * mean(ks)
df

hist(phihats2)
abline(v=phi,col='red')
abline(v=mean(phihats2,col='blue'))



paste(e2,e3,sep='   ')
length(which(!is.na(gamma.err)))

zz <- phi * rchisq(1000,df=df)/df
plot(density(phihats2))
lines(density(zz),col='red')
