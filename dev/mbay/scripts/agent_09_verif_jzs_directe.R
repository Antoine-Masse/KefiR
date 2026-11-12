.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor))
packageVersion("BayesFactor")
set.seed(20261105)
# 1. JZS deux groupes, ma propre formule : int (1+Ng)^-1/2 (1+t2/((1+Ng)v))^-(v+1)/2 / (1+t2/v)^-(v+1)/2 * pi(g), g ~ IG(1/2, r^2/2)
a<-rnorm(20,0,1); b<-rnorm(25,.6,1)
tt<-t.test(a,b,var.equal=TRUE); t2<-unname(tt$statistic)^2; v<-tt$parameter; N<-20*25/45; r<-sqrt(2)/2
f<-function(g) (1+N*g)^-.5*((1+t2/((1+N*g)*v))/(1+t2/v))^(-(v+1)/2)*(r/sqrt(2*pi))*g^-1.5*exp(-r^2/(2*g))
mine<-integrate(f,0,Inf,rel.tol=1e-10)$value
bf<-exp(ttestBF(a,b,rscale=r)@bayesFactor$bf)
cat("JZS moi",mine," ttestBF",bf," ecart rel",abs(mine/bf-1),"\n")
# 2. Inclusion : recalcul a la main sur A*B, tous modeles, a priori uniforme
d<-expand.grid(A=factor(1:3),B=factor(1:2),r=1:10); d$y<-rnorm(60)+as.numeric(d$A)*.5
m<-anovaBF(y~A*B,data=d,whichModels="all",progress=FALSE,iterations=50000)
x<-exp(m@bayesFactor$bf); names(x)<-rownames(m@bayesFactor); x<-c(x,Intercept=1)
print(round(x,4))
inc<-function(t) { has<-sapply(strsplit(sub(" .*","",""),""),function(z)NA); }
terms<-lapply(names(x),function(n) if(n=="Intercept") character() else strsplit(n," \+ ")[[1]])
for(t in c("A","B","A:B")){ h<-sapply(terms,function(z) t %in% z); cat(t,"P(inclus|D)=",sum(x[h])/sum(x)," BFincl=",(sum(x[h])/sum(h))/(sum(x[!h])/sum(!h)),"\n") }
# controle frequentiste a cote
print(anova(lm(y~A*B,d)))
