library(MASS)




### Negative binomial

# offsets
log.c <- log(seq(0.1,10,0.1))

true.mu= exp(log.c)
n <- length(true.mu)
n

# size = 1/dispersion 
# lam = (mu^2/(var-lam))
# so mu=10, var=20, phi=2 => lam=10
true.var <- 2 * true.mu
lam = (true.mu^2/(true.var-true.mu))

#true.var <- true.mu * (1+true.mu/lam)
true.phi <- true.var / true.mu

phi1 <- numeric()
phi2 <- numeric()
low.wald1 <- numeric()
upp.wald1 <- numeric()
low.wald2 <- numeric()
upp.wald2 <- numeric()
low.sc1 <- numeric()
upp.sc1 <- numeric()
low.sc2 <- numeric()
upp.sc2 <- numeric()
low.boot <- numeric()
upp.boot <- numeric()


par(mfrow=c(2,1))
n.sims = 100
# number of parameters
p <- 1

for (sim in 1:n.sims) {
  print(sim)
  y <- rnbinom(n,mu=true.mu,size=lam)
  
  m <- glm(y~log.c,family="quasipoisson")
  summary(m)
  
  mu <- predict(m,type='response')
  
  # variance equals the mean in the Poisson glm
  v <- mu
  #phi1 <- (sum(residuals(m,type='response')^2)/v)/(n-1)
  #phi1 <- sum(residuals(m,type="pearson")^2)/(n-1)
  # or equivalently
  phi1[sim] <- summary(m)$dispersion
  
  e <- residuals(m,type='response')/sqrt(phi1[sim]*v)
  rho3 <- e^3
  rho4 <- e^4 - 3
  
  
  X <- matrix(log.c,ncol=1)
  W <- matrix(0,nrow=n,ncol=n)
  for (index in 1:n){
    W[index,index] <- mu[index]
  }
  Q <- X %*% ( solve( t(X) %*% (W %*% X) ) %*% t(X) )
  q <- mu * rowSums(Q)
  
  
  #
  r <- sqrt(phi1[sim]/v)
  a <- rho4 + q * (phi1[sim] /mu - 2 * r * rho3)
  a.bar <- mean(a)
  
  
  nu <- ( 2*phi1[sim]^2/(n-p) ) * (1+0.5*a.bar)
  
  # wald interval
  low.wald1[sim] <- phi1[sim]-1.96*sqrt(nu)
  upp.wald1[sim] <- phi1[sim]+1.96*sqrt(nu)
  
  # skew corrected interval
  b <- r * q
  tau.bar <- mean(e^6 - 3 * b * e^5 - 15*e^4 - (b^3 - 3 * b^2 - 6*b - 12) * e^3 - 3 * b^2 - 6)
  k <- ( 8 * phi1[sim]^3 / (n-1)^2 ) * (1 + tau.bar/8)
  
  tl <- (3/k)*((1+k*(1.96-k/6))^(1/3)-1)
  tu <- (3/k)*((1+k*(-1.96-k/6))^(1/3)-1)
  
  low.sc1[sim] <- phi1[sim]-tl*sqrt(nu)
  upp.sc1[sim] <- phi1[sim]-tu*sqrt(nu)
  
  ###
  ### phi2
  ###
  
  #   q is different, and we are now using phi2
  s <- mean(residuals(m,type='response')/v)
  phi2[sim] <- phi1[sim]/(1+s)
  q <- rep(1,n)
  r <- sqrt(phi2[sim]/v)
  a <- rho4 + q * (phi2[sim] /mu - 2 * r*rho3)
  a.bar <- mean(a)
  
  nu <- 2*phi2[sim]^2/(n-p) * (1+0.5*a.bar)
  
  # wald interval
  low.wald2[sim] <- phi2[sim]-1.96*sqrt(nu)
  upp.wald2[sim] <- phi2[sim]+1.96*sqrt(nu)
  
  # skewness corrected interval.
  
  b <- r * q
  tau.bar <- mean(e^6 - 3 * b * e^5 - 15*e^4 - (b^3 - 3 * b^2 - 6*b - 12) * e^3 - 3 * b^2 - 6)
  k <- ( 8 * phi2[sim]^3 / (n-1)^2 ) * (1 + tau.bar/8)
  
  tl <- (3/k)*((1+k*(1.96-k/6))^(1/3)-1)
  tu <- (3/k)*((1+k*(-1.96-k/6))^(1/3)-1)
  
  low.sc2[sim] <- phi2[sim]-tl*sqrt(nu)
  upp.sc2[sim] <- phi2[sim]-tu*sqrt(nu)
  
  # bootstrap interval
  n.boot = 10000
  phi.b <- numeric()
  for (b in 1:n.boot) {
    s <- sample.int(n,replace = T)
    y.s <- y[s]
    x.s <- log.c[s]
    m <- glm(y.s~x.s,family="quasipoisson")
    phi.b[b] <- summary(m)$dispersion
  }
  low.boot[sim] <- quantile(phi.b,0.025)
  upp.boot[sim] <- quantile(phi.b,0.975)
}


###
### Plot the intervals
###

par(mfrow=c(3,2))
##
## phi1 Walt interval
##
plot(0,pch='',xlim=c(0,n.sims),ylim=c(min(low.wald1),max(upp.wald1)),main='phi1, Wald')
errs <- 0
for (i in 1:n.sims) {
  col = 'black'
  # true.phi is the same for all data, so we just take the first
  if (true.phi[1] < low.wald1[i] | true.phi[1] > upp.wald1[i]) {
    errs <- errs+1
    col = 'red'
  }
  segments(i,low.wald1[i],i,upp.wald1[i],col=col)
}
abline(h=true.phi,col='blue')
errs/n.sims

##
## phi1 Skewness corrected interval
##
plot(0,pch='',xlim=c(0,n.sims),ylim=c(min(low.sc1),max(upp.sc1)),main='phi1, Skew')
errs <- 0
for (i in 1:n.sims) {
  col = 'black'
  # true.phi is the same for all data
  if (true.phi[1] < low.sc1[i] | true.phi[1] > upp.sc1[i]) {
    errs <- errs+1
    col = 'red'
  }
  segments(i,low.sc1[i],i,upp.sc1[i],col=col)
}
abline(h=true.phi,col='blue')
errs/n.sims

####
#### phi2, Wald interval
####
plot(0,pch='',xlim=c(0,n.sims),ylim=c(min(low.wald2),max(upp.wald2)),main='phi2, Wald')
errs <- 0
for (i in 1:n.sims) {
  col = 'black'
  # true.phi is the same for all data
  if (true.phi[1] < low.wald2[i] | true.phi[1] > upp.wald2[i]) {
    errs <- errs+1
    col = 'red'
  }
  segments(i,low.wald2[i],i,upp.wald2[i],col=col)
}
abline(h=true.phi,col='blue')
errs/n.sims

####
#### phi2, Skewness corrected interval
####
plot(0,pch='',xlim=c(0,n.sims),ylim=c(min(low.sc2),max(upp.sc2)),main='phi2, Skew')
errs <- 0
for (i in 1:n.sims) {
  col = 'black'
  # true.phi is the same for all data
  if (true.phi[1] < low.sc2[i] | true.phi[1] > upp.sc2[i]) {
    errs <- errs+1
    col = 'red'
  }
  segments(i,low.sc2[i],i,upp.sc2[i],col=col)
}
abline(h=true.phi,col='blue')
errs/n.sims

###
### bootstrap interval
###
plot(0,pch='',xlim=c(0,n.sims),ylim=c(min(low.sc2),max(upp.sc2)),main='phi2, Skew')
errs <- 0
for (i in 1:n.sims) {
  col = 'black'
  # true.phi is the same for all data
  if (true.phi[1] < low.boot[i] | true.phi[1] > upp.boot[i]) {
    errs <- errs+1
    col = 'red'
  }
  segments(i,low.boot[i],i,upp.boot[i],col=col)
}
abline(h=true.phi,col='blue')
errs/n.sims


hist(phi1,30)
hist(phi2,30)


mean(phi1)
var(phi1)
mean(phi2)
var(phi2)
summary(m)





