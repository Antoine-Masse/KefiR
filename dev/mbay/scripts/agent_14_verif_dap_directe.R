.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); set.seed(909)
k<-3;n<-10;r<-.5;P<-eigen(diag(k)-1/k)$vectors[,1:(k-1)]; g<-factor(rep(1:k,each=n)); X<-model.matrix(~g-1)%*%P
sim<-function(){h<-rbinom(1,1,.5);gg<-1/rgamma(1,.5,r^2/2);b<-if(h) rnorm(k-1,0,sqrt(gg)) else rep(0,k-1);y<-X%*%b+rnorm(k*n)
 F<-anova(lm(y~g))[1,4];bf<-exp(oneWayAOV.Fstat(F,n,k,rscale=r)$bf);c(h=h,p=bf/(1+bf))}
R<-t(replicate(2000,sim()));cat("N=2000 moyenne P=",mean(R[,2])," se=",sd(R[,2])/sqrt(2000),"\n")
b<-cut(R[,2],c(0,.2,.4,.6,.8,1),include.lowest=TRUE);print(round(cbind(annonce=tapply(R[,2],b,mean),freq=tapply(R[,1],b,mean),n=table(b)),3))
