results <- read.csv('~/Documents/overdispersion/results.csv',
               stringsAsFactors = F,header=F)
head(results)

colnames(results)=c("Population",'phi','n','b1','b2',"N","chisq","gamma",
               "bootstrap","vgam")
head(results)
dim(results)


resultsGamma <- read.csv('~/Documents/overdispersion/resultsgamma.csv',
                    stringsAsFactors = F,header=F)
colnames(resultsGamma)=c("Population",'phi','n','b1','b2',"N","chisq","gamma",
                    "bootstrap","vgam")
dim(resultsGamma)

results$gamma <- resultsGamma$gamma
head(results[,"gamma"])
head(results[,"chisq"])



## add mu lower and upper, calculated from the betas
results$mu.lower <- NA
results$mu.upper <- NA

for (i in 1:dim(results)[1]) {
    x<-seq(0,1,length.out=results$n[i])
    results$mu.lower[i]  <- range(exp(results$b1[i] + results$b2[i] * x))[1]
    results$mu.upper[i]  <- range(exp(results$b1[i] + results$b2[i] * x))[2]
}
head(results)



results.ee <- read.csv('~/Documents/overdispersion/results-ee2.csv',
                         stringsAsFactors = F,header=F)
colnames(results.ee)=c("Population",'phi','n','b1','b2',"N","ee")
dim(results.ee)
head(results.ee)

results$ee <- results.ee$ee
head(results[,"ee"])

## take chisq and gamma from resultsgamma2, because the previous ones
##  used 90% CI, not 95%

resultsGamma2 <- read.csv('~/Documents/overdispersion/resultsgamma2.csv',
                         stringsAsFactors = F,header=F)
colnames(resultsGamma2)=c("Population",'phi','n','b1','b2',"N","chisq","gamma",
                         "bootstrap","vgam")
dim(resultsGamma2)

results$chisq <- resultsGamma2$chisq
results$gamma <- resultsGamma2$gamma
head(results[,"gamma"])
head(results[,"chisq"])


## drop duplicates (phi == 1)
dropm <- which(is.na(results$ee))
results <- results[-dropm,]
