N <- 10
var.true <- 10
N.sims <- 10000
s2 <- vector()
bias <- vector()
for (sim in 1:N.sims) {
  x <- rnorm(N,0,sqrt(var.true))
  x.bar <- mean(x)
  s2[sim] <- (1/N)*sum((x-x.bar)^2)
  bias[sim] <- (s2[sim]+(s2[sim]+s2[sim]/N)/N)/N
}
hist(s2-var.true)
abline(v=var.true/N,col='red')

par(mfrow=c(2,1))
hist(bias,probability = T)
abline(v=5/N,col='red')
mean(bias)

hist(s2 + bias,probability = T)
mean(s2)
mean(s2 + bias)

var(bias)



c<- seq(0.01,0.5,0.01)
y2 <- dgamma(c,shape=0.05,scale=2)
y3 <- rgamma(10000,shape=0.05,scale=2)
mean(y3)
var(y3)
hist(y3,probability = T,50)
lines(c,y2,col='red') 


mean(s2-1)

