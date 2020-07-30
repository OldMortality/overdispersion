##
## bootstrap done by refitting model
##

setwd('~/overdispersion')
results <- read.csv('results-mainpar10.csv',
               stringsAsFactors = F,header=F)


dim(results)

head(results)
dim(results)
results[1,]

colnames(results)=c("Population",'phi','n','b1','b2',"N",
                    "cover.phi","cover.beta",
                    "pow.phi","pow.beta",
                    "med.phi","med.beta",
                    "err"
                    )
head(results)
dim(results)

# results$cover.chisq <- NA
# results$pow.chisq <- NA
# results$med.chisq <- NA
# results$err.chisq <- NA

## drop duplicates (phi == 1)
dropm <- which(is.na(results$N))
if (length(dropm) > 0) {
  results <- results[-dropm,]
}


## add mu lower and upper, calculated from the betas
results$mu.lower <- NA
results$mu.upper <- NA

for (i in 1:dim(results)[1]) {
    x<-seq(0,1,length.out=results$N[i])
    results$mu.lower[i]  <- range(exp(results$b1[i] + results$b2[i] * x))[1]
    results$mu.upper[i]  <- range(exp(results$b1[i] + results$b2[i] * x))[2]
}
head(results)

#results[which(results$Population==1),"Population"] <- "Negbin"
#results[which(results$Population==2),"Population"] <- "Neyman"
#results[which(results$Population==3),"Population"] <- "Pois lognorm"


