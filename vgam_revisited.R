library(MASS)     # for "rnegbin"
library(rmutil)   # for "rinvgauss"
library(VGAM)

# sample size
n <- 30
#beta <- c(-3,3) # sparse
beta <- c(2.3,0.7)
x <- seq(0,1,length.out=30)
eta<-beta[1]+beta[2]*x
mu<-exp(eta)

vgam.upp <- vector()
vgam.low <- vector()

phi <- 1

N = 100 # saves time, you find pretty much the same thing as
         #   with N=1e4

count.warnings <- 0
for (i in 1:N) {
  theCI <- tryCatch(
    {
      print('=======')
      print(i)
      print(y)
      print('========')
      #y<-rpois(n,mu) 
      a3<-vglm(y~x,family=negbinomial(parallel=T,zero=NULL))
      
      #plot(y~x)
      phihat3<-Confint.nb1(a3)$phi0
      # CI for phi
      phihat3ci<-Confint.nb1(a3,level=0.90)$CI.phi0
      vgam.low[i] <- phihat3ci[1]
      vgam.upp[i] <- phihat3ci[2]
      c(phihat3ci[1],phihat3ci[2])
    }
    ,
    error=function(cond) {
      
      message(paste('rr',y,x,cond,sep=' '))
      return(NA)
    },
    warning=function(cond) {
      print(i)
      print('warning.............')
      count.warnings <<- count.warnings + 1
      return(NA)
    }
    
  ) # tryCatch
  print(theCI)
  if (!is.na(theCI[1])) {
    vgam.low[i] <- theCI[1]
    vgam.upp[i] <- theCI[2]
  } else {
    vgam.low[i] <- NA
    vgam.upp[i] <- NA
  }
} # for i=

# example with mad result without warning
y = 
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,2,2,0,0,0,1,0,0,2)



print(count.warnings)
sum(is.na(vgam.upp))
#sum(vgam.low == 1,na.rm=T)
sum(vgam.upp > 1000,na.rm=T)
dropm <- which(vgam.upp > 1000)
length(dropm)
vgam.low2 <- vgam.low[-dropm]
vgam.upp2 <- vgam.upp[-dropm]
coverage <- sum(vgam.low2 <= phi & vgam.upp2 >= phi)
min(vgam.low2)
coverage




#
# phi != 3
#
phi = 5
nu <- phi-1
mn<-mu/nu

for (i in 1:N) {
  
  tryCatch(
    {
    # exactly the same as before, except this line
    y<-rnegbin(n,mu,mn) 
    #
    a3<-vglm(y~x,family=negbinomial(parallel=T,zero=NULL))
    length(warnings)
    # estimate for phi
    #plot(y~x)
    phihat3<-Confint.nb1(a3)$phi0
    # CI for phi
    phihat3ci<-Confint.nb1(a3,level=0.90)$CI.phi0
    vgam.low[i] <- phihat3ci[1]
    vgam.upp[i] <- phihat3ci[2]
    }
    ,
    error=function(cond) {
      message(paste('rr',y,x,cond,sep=' '))
      return(NA)
    } 
    
  )
}

sum(is.na(vgam.upp))
sum(vgam.low == 1,na.rm=T)
sum(vgam.upp > 1e5,na.rm=T)
dropm <- which(vgam.upp > 1000)
length(dropm)
if (length(dropm)>1) {
  vgam.low2 <- vgam.low[-dropm]
  vgam.upp2 <- vgam.upp[-dropm]
} else {
  vgam.low2 <- vgam.low
  vgam.upp2 <- vgam.upp
}
min(vgam.low2)
coverage <- sum(vgam.low2 <= phi & vgam.upp2 >= phi)
coverage/(N-length(dropm))
phi

# pretty mad
table(vgam.low2)
table(vgam.upp2)



