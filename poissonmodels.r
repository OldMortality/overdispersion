library(MASS)

mu = 10
N = 5000
n.sim <- 1000


qs <- numeric()
disp <- numeric()
for (i in 1:n.sim) {
  d <- rpois(N,mu)
  qs[i] <- sum((d-mu)^2/mu)
  disp[i] <- qs[i]/(N-1)
}

plot(density(qs))
curve(dchisq(x,df=N-1),from=18000,to=22000,col='red',add=T)

plot(density(disp))

mean(qs)
N-1
var(qs)
2 * (N-1)

# so we see that the sum of the squared Pearson residuals ~ chisq(df=N-1)
# qs ~ chisq(df=N-1), so mean(qs) = N-1, and the var(qs)=2*(N-1)
#

# estimate for dispersion parameter.
mean(qs)/(N-1)

## same thing, with overdispersed data
qs <- numeric()
disp <- numeric()
vars <- numeric()
means <- numeric()
mu <- 5
N <- 25
s.bar <- numeric()
disp2 <- numeric()

for (i in 1:n.sim) {
  d <- 2 * rpois(N,mu)
  vars[i] <- var(d)
  means[i] <- mean(d)
  
  m <- glm(d~1,family='quasipoisson')
  m.hat <- exp(coefficients(m)[1])
  qs[i] <- sum((d-m.hat)^2/m.hat)
  disp[i] <- qs[i]/(N-1)
  s.bar[i] <- mean((d-m.hat)/m.hat)
  disp2[i] <- disp[i]/(1+s.bar[i])
}
plot(vars~means)
plot(density(vars/means))

#mean(qs)
#mean(qs)/(N-1)
par(mfrow=c(2,1))
plot(density(vars/means))

plot(density(disp))
plot(density(disp2))
#plot(density(disp),add=T,col='red')
mean(disp)

par(mfrow=c(1,1))
# these are the same
plot(vars/means~disp)
head(vars/means)
head(disp)


N <- 5
n.sim <- 1000
vars <- numeric()
for (i in 1:n.sim) {
  x <- rnorm(N,0,1)
  vars[i] <- var(x)
}
plot(density(vars))
mean(vars)

N <- 100
x <- runif(N,10,20)
mu <- 1 + 3*x
boots <- numeric()
# simulated phi hat's
phis <- numeric()
# true phi
true.phi <- 5

n.sim <- 5000
for (i in 1:n.sim) {
  y <- true.phi*rpois(N,mu)

  m <- glm(y~x,family='quasipoisson')
  phis[i] <- summary(m)$dispersion
}
plot(phis)
par(mfrow=c(2,2))
# asymptotically chi-sq?
chi <- (N-2)*phis/true.phi
plot(density(chi),main='simulated chi=(N-2)*phi.hat/phi')
curve(dchisq(x,df=N-2),add=T,col='red')
mean(phis)

plot(density(phis))
mean(phis)
sqrt(var(phis))

# bootstrap interval for final simulation

# subsampling
for (i in 1:1000) {
  d <- data.frame(x,y)
  s <- d[sample(1:nrow(d),N,replace=T),]
  
  m <- glm(y~x,family='quasipoisson',data=s)
  boots[i] <- summary(m)$dispersion
  
}
chi <- (N-2)*boots/true.phi
plot(chi)
plot(density(chi),main='bootstrap')
curve(dchisq(x,df=N-2),add=T,col='red')


par(mfrow=c(2,2))
plot(density(chi))
curve(dchisq(x,df=N-2),add=T,col='red')
plot(density(boots))
plot(density(phis))

mean(boots)
var(boots)
var(phis)



z <- 2 * rpois(1000,5)
mean(z)
var(z)
m <- glm(z~1, family='quasipoisson')
summary(m)$dispersion
