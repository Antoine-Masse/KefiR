.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); options(BFprogress=FALSE)
nT<-6;ns<-20
run<-function(Sig,R=100){L<-chol(Sig);z<-t(replicate(R,{Y<-matrix(rnorm(ns*nT),ns)%*%L;dl<-data.frame(y=c(Y),T=factor(rep(1:nT,each=ns)),id=factor(rep(1:ns,nT)))
 b<-extractBF(lmBF(y~T+id,dl,whichRandom="id",iterations=2000)/lmBF(y~id,dl,whichRandom="id",iterations=2000))$bf
 Yc<-Y-rowMeans(Y);Yc<-t(t(Yc)-colMeans(Y));sse<-sum(Yc^2);sst<-ns*sum((colMeans(Y)-mean(Y))^2);F<-(sst/(nT-1))/(sse/((nT-1)*(ns-1)));c(b/(1+b),1-pf(F,nT-1,(nT-1)*(ns-1)))}))
 c(oui=mean(z[,1]>.5),fort=mean(z[,1]>.95),medP=median(z[,1]),p05=mean(z[,2]<.05))}
set.seed(71)
print(round(rbind(CS=run(matrix(.7,nT,nT)+diag(.3,nT)),AR95=run(.95^abs(outer(1:nT,1:nT,"-"))),AR80=run(.8^abs(outer(1:nT,1:nT,"-")))),3))
