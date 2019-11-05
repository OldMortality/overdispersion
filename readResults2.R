##
## Changes:
##  - added boots2
##  - reading 3 files (for n=30, n=100, n=1000)
##

results30 <- read.csv('~/Documents/overdispersion/results30.csv',
               stringsAsFactors = F,header=F)
results100 <- read.csv('~/Documents/overdispersion/results100.csv',
                      stringsAsFactors = F,header=F)
results1000 <- read.csv('~/Documents/overdispersion/results1000.csv',
                      stringsAsFactors = F,header=F)
results <- rbind(results30,results100,results1000)
dim(results)

head(results)
dim(results)
results[1,]

colnames(results)=c("Population",'phi','n','b1','b2',"N",
                    "cov.chisq","cov.gamma","cov.boots1","cov.boots2","cov.vgam","cov.ee",
                    "pow.chisq","pow.gamma","pow.boots1","pow.boots2","pow.vgam","pow.ee",
                    "med.chisq","med.gamma","med.boots1","med.boots2","med.vgam","med.ee",
                    "err.chisq","err.gamma","err.boots1","err.boots2","err.vgam","err.ee"
                    )
head(results)
dim(results)



## add mu lower and upper, calculated from the betas
results$mu.lower <- NA
results$mu.upper <- NA

for (i in 1:dim(results)[1]) {
    x<-seq(0,1,length.out=results$n[i])
    results$mu.lower[i]  <- range(exp(results$b1[i] + results$b2[i] * x))[1]
    results$mu.upper[i]  <- range(exp(results$b1[i] + results$b2[i] * x))[2]
}
head(results)



## drop duplicates (phi == 1)
dropm <- which(is.na(results$ee))
if (length(dropm) > 0) {
  results <- results[-dropm,]
}
