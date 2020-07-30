dfs<-c(5,10,20)
tau.bars<-c(0,1,2)

par(mfcol=c(length(dfs),length(tau.bars)))

for (i in 1:length(dfs)){

  df<-dfs[i]
  x<-seq(0,qchisq(0.9999,df),0.01)
  
  for (j in 1:length(tau.bars)){
  
    tau.bar<-tau.bars[j]
    
    theta<-2+tau.bar
    k<-df/theta
    
    y1<-dchisq(x,df)
    y2<-dgamma(x,shape=k,scale=theta)
    y3<-dlnorm(x,log(df),sqrt(theta/df))
    
    plot(x,y1,typ="l",xlab="",ylab="",yaxt="n",main=paste("n-p = ",df,", tau.bar = ",tau.bar))
    lines(x,y2,col="blue")
    lines(x,y3,col="red")  
    
  }
  
  
}

