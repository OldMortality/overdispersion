library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"

# calculate CI's for phi, based on chisquare distribution
#   with phihat1 and phihat2, and on Gamma distribution
#   (phihat2 only), and fitting only Negative binomial 
#   for now.
#   Plots first 100 of the simulations, with error rates.
#

par(mfrow=c(1,1))

phi <- 2
  
  
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
    
    N <- 2000
    
    # estimates of phi
    phihats2 <- vector()
    ks <- vector()
    thetas <- vector()
  }
  
  alpha3.hats <- vector()
  bias.true <- vector()
  bias.hat <- vector()
  thetas <- vector()
  ks <-  vector()
  alpha3.true <- phi*(2*phi-1)
  for (sim in 1:N) { 
    print(sim)
    y<-rnegbin(n,mu,mn) # negative binomial
      
    muhat<-fitted(glm(y~x,family="poisson"))
    P<-sum((y-muhat)^2/muhat)
    sbar<-mean((y-muhat)/muhat) 
    phihat1<-P/(n-p)
    phihat2<-phihat1/(1+sbar) 
    phihats2[sim] <- phihat2
    
    # chisquare CI
    df <- n-p
    #
    e <- y-muhat
    alpha3.hat <- (1/df) * sum( (e^3/muhat))
    alpha3.hats[sim] <- alpha3.hat
    
    # work out S
    
    W <- diag(muhat)
    W.true <- diag(mu)
    X1 <- rep(1,n)
    X2 <- x
    X <- cbind(X1,X2)
    Q <- X %*% solve(t(X) %*% W %*% X) %*% t(X)
    Q.true <- X %*% solve(t(X) %*% W.true %*% X) %*% t(X)
    S <- sum(1/muhat) + n * sum(diag(Q))-sum(Q)
    
    S.true <- sum(1/mu) + n * sum(diag(Q.true))-sum(Q.true)
    bias.true[sim] <- (alpha3.true/phi-phi)*S.true/n
    
    bias.hat[sim] <- (alpha3.hat/phihat2-phihat2)*S/n
    
    ktheta <- df - bias.hat[sim]
    ktheta2 <- 2 * ktheta 
    theta <- ktheta2 / ktheta
    k <- ktheta / theta
    thetas[sim] <- theta
    ks[sim] <- k
        
  }

  
  mean(ktheta)
  mean(bias.hat)
  hist(bias.hat,60,
       
       main=paste('phi=',phi,'n=',n,sep=' '))
  abline(v=bias.true[1],col='red')
  

par(mfrow=c(1,1))
for (i in 1:9) { 
hist((n-p)*phihats2/phi,probability = T,
     ylim=c(0,0.4),
     xlim=c(0,2 * n),100)
c<- seq(0,qchisq(0.9999,df),0.1)
y1<-dchisq(c,df=n-p)
lines(c,y1,col='red') 
y2 <- dgamma(c,shape=ks[112],scale=thetas[112])
lines(c,y2,col='blue') 
}
par(mfrow=c(1,1))
hist(alpha3.hats,60)
abline(v=alpha3.true,col='red')
mean(alpha3.hats)
alpha3.true

hist(bias.true,60)  
  
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


