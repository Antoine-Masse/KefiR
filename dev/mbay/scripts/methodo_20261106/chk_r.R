.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib",.libPaths()))
suppressMessages(library(BayesFactor)); set.seed(20261106)
for(cfg in list(c(10,20,30),c(15,25))){
k<-length(cfg); y<-unlist(lapply(seq_len(k),function(i)rnorm(cfg[i],(i-1)*.4,1))); f<-factor(rep(seq_len(k),cfg))
b<-anovaBF(y~f,data=data.frame(y=y,f=f),rscaleFixed=.5,progress=FALSE,iterations=50000)
cat("R",cfg,"logBF",extractBF(b)$bf|>log(),"\n"); write.csv(data.frame(y=y,f=as.integer(f)),paste0("d",k,".csv"),row.names=FALSE)}
