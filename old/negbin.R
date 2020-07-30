library(MASS)

# prob of failure
p <- 1/6
r <- 1

x <- rnbinom(250,size=r,p=p)
mean(x)
r*(1-p)/p
var(x)
r * (1-p)/p^2
var(x)

# 1/p
var(x)/mean(x)

hist(x)
m <- glm(x~1,family="quasipoisson")
summary(m)


mu=10
# size = 1/dispersion 
lam = 10
x <- rnbinom(100,mu=mu,size=lam)
mean(x)
# lam = (mu^2/(var-lam))
# so mu=10, var=20, phi=2 => lam=10
var(x)






par(mfrow=c(2,1))
# small size, large variance
plot(rnbinom(100,mu=mu,size=2),ylim=c(0,15),main="size 2")
# large size, small variance
plot(rnbinom(100,mu=mu,size=10),ylim=c(0,15),main="size 10")




mean(x)
var(x)
# should be var(x)
mu*(1+mu/lam)

x <- rnbinom(100000,mu=mu,size=2)
mean(x)
var(x)
library(MASS)
m2 <- glm.nb(x~1)
summary(m2)
m2$theta
exp(1.1)
# the size, not the dispersion!
summary(m2)


x <- rnbinom(100000,p=2/5,size=2)
mean(x)
var(x)
library(MASS)
m2 <- glm.nb(x~1,link="log")
summary(m2)



N <- 100
phi1 <- numeric()
phi2 <- numeric()
n.sims <- 1
0
for (i in 1:n.sims) {
  for (phi in seq(1,1,0.1)) {
    x <- rnbinom(N,mu=1,size=1/phi)
    m <- glm(x~1,family="quasipoisson")
    s <- sum(residuals(m, type = "pearson")^2)
    phi1[sims] <- s /(N-1)
    #
    # is V'== 1 ?
    mu.hat <- exp(coefficients(m))[[1]] 
    phi2[sim] <- phi1[sim]/(1+mean((x-mu.hat)/mu.hat))
  }
}
par(mfrow=c(2,1))
hist(phi1,30)
hist(phi2,30)



m10 <- glm.nb(rnbinom(1000,mu=mu,size=10)~1)

m10$theta

summary(m10)



true.mu=10
# size = 1/dispersion 
# lam = (mu^2/(var-lam))
# so mu=10, var=20, phi=2 => lam=10
lam = 10

true.var <- true.mu * (1+true.mu/lam)
true.phi <- true.var / mu

phi1 <- numeric()
phi2 <- numeric()
low <- numeric()
upp <- numeric()
nu <- numeric()
par(mfrow=c(2,1))
n.sims = 1000
n <- 50
for (sim in 1:n.sims) {
  x <- rnbinom(n,mu=true.mu,size=lam)
  mean(x)
  var(x)
  m <- glm(x~1,family="quasipoisson")
  
  mu <- exp(coefficients(m)[[1]])
  v <- mu
  #phi1 <- (sum(residuals(m,type='response')^2)/v)/(n-1)
  #phi1 <- sum(residuals(m,type="pearson")^2)/(n-1)
  # or equivalently
  phi1[sim] <- summary(m)$dispersion
  
  e <- residuals(m,type='response')/sqrt(phi1[sim]*v)
  rho3 <- e^3
  rho4 <- e^4
  
  
  
  
  #W <- matrix(0,nrow=n,ncol=n)
  #for (m in 1:n) {
  #  W[m,m] <- mu
  #}
  #Q <- solve(W) 
  
  q <- rep(1,n)
  r <- sqrt(phi1[sim]/v)
  a <- rho4 + q * (phi1[sim] * 1 - 2 * r*rho3)
  a.bar <- mean(a)
  nu[sim] <-2*phi1[sim]^2*(1+a.bar/2)/(n-1)
  
  # wald interval
  low[sim] <- phi1[sim]-1.96*sqrt(nu[sim])
  upp[sim] <- phi1[sim]+1.96*sqrt(nu[sim])
  
  ###
  ### phi2
  ###
  # variance is the mean is exp(linear predictor)
  v <- exp(coefficients(m)[[1]])
  s <- mean(residuals(m,type='response')/v)
  phi2[sim] <- phi1[sim]/(1+s)
}


plot(0,pch='',xlim=c(0,n.sims),ylim=c(min(low),max(upp)))
errs <- 0
for (i in 1:n.sims) {
  col = 'black'
  if (true.phi < low[i] | true.phi > upp[i]) {
    errs <- errs+1
    col = 'red'
  }
  segments(i,low[i],i,upp[i],col=col)
}
abline(h=true.phi,col='blue')
errs/n.sims



hist(phi1,30)
hist(phi2,30)

summary(m)


### not all mu-i the same


log.c <- log(seq(0.1,10,0.1))

true.mu= exp(log.c)
n <- length(true.mu)

# size = 1/dispersion 
# lam = (mu^2/(var-lam))
# so mu=10, var=20, phi=2 => lam=10
true.var <- 2 * true.mu
lam = (true.mu^2/(true.var-true.mu))

#true.var <- true.mu * (1+true.mu/lam)
true.phi <- true.var / true.mu

phi1 <- numeric()
phi2 <- numeric()
low <- numeric()
upp <- numeric()
low2 <- numeric()
upp2 <- numeric()

nu <- numeric()
par(mfrow=c(2,1))
n.sims = 1000

for (sim in 1:n.sims) {
  y <- rnbinom(n,mu=true.mu,size=lam)
  
  m <- glm(y~log.c,family="quasipoisson")
  summary(m)
  
  mu <- predict(m,type='response')
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
  q <- mu * colSums(Q)
  
  
  #
  r <- sqrt(phi1[sim]/v)
  a <- rho4 + q * (phi1[sim] * (1/mu) - 2 * r * rho3)
  a.bar <- mean(a)
  
  nu[sim] <- 2*phi1[sim]^2/(n-p)* (1+0.5*a.bar)
  
  # wald interval
  low[sim] <- phi1[sim]-1.96*sqrt(nu[sim])
  upp[sim] <- phi1[sim]+1.96*sqrt(nu[sim])
  
  
  ###
  ### phi2
  ###
  # variance is the mean is exp(linear predictor)
  
  s <- mean(residuals(m,type='response')/v)
  phi2[sim] <- phi1[sim]/(1+s)
  q <- rep(1,n)
  r <- sqrt(phi2[sim]/v)
  a <- rho4 + q * (phi2[sim] * (1/mu) - 2 * r*rho3)
  a.bar <- mean(a)
  nu[sim] <- 2*phi2[sim]^2/(n-p) * (1+0.5*a.bar)
  
  # wald interval
  low2[sim] <- phi2[sim]-1.96*sqrt(nu[sim])
  upp2[sim] <- phi2[sim]+1.96*sqrt(nu[sim])
  
}


plot(0,pch='',xlim=c(0,n.sims),ylim=c(min(low),max(upp)))
errs <- 0
for (i in 1:n.sims) {
  col = 'black'
  # true.phi is the same for all data
  if (true.phi[1] < low[i] | true.phi[1] > upp[i]) {
    errs <- errs+1
    col = 'red'
  }
  segments(i,low[i],i,upp[i],col=col)
}
abline(h=true.phi,col='blue')
errs/n.sims
# oops, min(low)<0
min(low)
max(low)
min(upp)
max(upp)


plot(0,pch='',xlim=c(0,n.sims),ylim=c(min(low2),max(upp2)))
errs <- 0
for (i in 1:n.sims) {
  col = 'black'
  if (true.phi[1] < low2[i] | true.phi[1] > upp2[i]) {
    errs <- errs+1
    col = 'red'
  }
  segments(i,low2[i],i,upp2[i],col=col)
}
abline(h=true.phi,col='blue')
errs/n.sims



hist(phi1,30)
hist(phi2,30)

summary(m)





