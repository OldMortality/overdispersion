# compares CI's for phihat2 for the famous ecklonia dataset
# compares t-value and generalised t-value (Gamma instead of chisq)

f <- read.csv('~/Documents/overdispersion/ecklonia.txt',sep='\t',
              header=T)
head(f)
m <- glm(eck~kina,family="poisson",data=f)
summary(m)
 

dim(f)

y <- f$eck
x <- f$kina
df <- m$df.residual

muhat<-fitted(glm(y~x,family="poisson"))
P<-sum((y-muhat)^2/muhat)
sbar<-mean((y-muhat)/muhat) 
phihat1<-P/df
phihat2<-phihat1/(1+sbar) 

# chisquare CI

lower.chisq <- df * phihat2  / qchisq(0.95,df=df)
upper.chisq <- df * phihat2  / qchisq(0.05,df=df)
round(c(lower.chisq,upper.chisq),1)

#
e <- y-muhat
alpha3.hat <- (1/df) * sum( (e^3/muhat))
alpha4.hat <- (1/df) * sum(e^4/muhat) - 3 * muhat* phihat2^2

#tau.i <- (alpha4.hat/phi^2 - 2 * alpha3.hat/phi + phi )
tau.i <- (alpha4.hat/phihat2^2 - 2 * alpha3.hat/phihat2 + phihat2 )

tau.i <- (1/muhat) * tau.i
tau = mean(tau.i)

# work out S
eta <- log(muhat)
W <- diag(muhat)
n <- dim(f)[1]
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

gamma.lower <- df * phihat2  / qgamma(0.95,
                                           shape=k,
                                           scale=theta)

gamma.upper <- df * phihat2  / qgamma(0.05,
                                           shape=k,
                                           scale=theta)

round(c(gamma.lower,gamma.upper),1)
round(c(lower.chisq,upper.chisq),1)

N <- 10000
num <- rnorm(N)
denom1 <- sqrt(rchisq(N,df=df)/df)
denom2 <- sqrt(rgamma(N,shape=k,scale=theta)/df)
t1 <- num/denom1
t2 <- num/denom2

c(quantile(t1,0.025),quantile(t1,0.975))
c(quantile(t2,0.025),quantile(t2,0.975))

par(mfrow=c(2,1))
hist(denom1)
hist(denom2)
