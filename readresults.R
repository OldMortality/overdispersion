tt <- read.csv('~/Documents/overdispersion/results.csv',
               stringsAsFactors = F,header=F)
head(tt)

colnames(tt)=c("Population",'phi','n','b1','b2',"N","chisq","gamma")
head(tt)


tt$mu.lower <- NA
tt.mu.upper <- NA

for (i in 1:dim(tt)[1]) {
    x<-seq(0,1,length.out=tt$n[i])
    tt$mu.lower[i]  <- range(exp(tt$b1[i] + tt$b2[i] * x))[1]
    tt$mu.upper[i]  <- range(exp(tt$b1[i] + tt$b2[i] * x))[2]
}
head(tt)




