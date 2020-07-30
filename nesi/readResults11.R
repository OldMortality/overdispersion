##
## Changes:
##  - added boots2
##  - reading 3 files (for n=30, n=100, n=1000)
##
rm(list=ls())
setwd('~/Documents/overdispersion/nesi')
results <- read.csv('results-mainpar11.csv',
               stringsAsFactors = F,header=F)


dim(results)

head(results)
dim(results)
results[1,]

colnames(results)=c("Population",'phi','n','b1','b2',"N",
                    "cover.beta.chisq","cover.beta.vgam","cover.phi.chisq","cover.phi.vgam",
                    "power.beta.chisq","power.beta.vgam","power.phi.chisq","power.phi.vgam",
                    "width.beta.chisq","width.beta.vgam","width.phi.chisq","width.phi.vgam",
                    "chisq.na","vgam.na")
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
