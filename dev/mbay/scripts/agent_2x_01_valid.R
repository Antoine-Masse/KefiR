source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/agent_2x_commun.R")
set.seed(20261101)
print(oneWayAOV.Fstat); 
for (cfg in list(c(3,10),c(3,30),c(2,10),c(2,30))) { k<-cfg[1]; n<-cfg[2]; f<-gl(k,n); dd<-c(); dr<-c()
 for(i in 1:100){ al<-ralpha(k); y<-rnorm(n*k)+rep(al,each=n); F<-Fstat(y,n,k); s<-lbf(y,f); q<-lbf_F(F,n,k)
   N<-n*k;Q<-k-1; dd<-c(dd,s-q); dr<-c(dr, lbf_repli(Q*F/(Q*F+N-1-Q),N,Q,n)-q)}
 cat("k",k,"n",n,"max|anovaBF/ttestBF - voie rapide| =",max(abs(dd),na.rm=TRUE)," NA:",sum(is.na(dd)),
     " max|repli - rapide| =",max(abs(dr),na.rm=TRUE),"\n")
 cat("  aov p vs pF :", {y<-rnorm(n*k); paov(y,f)-pF(Fstat(y,n,k),n,k)},"\n")
}
f<-gl(3,10); y<-rnorm(30); cat("r=1 :",lbf_anova(y,f,1),lbf_F(Fstat(y,10,3),10,3,1)," ; k2 r=.25 :",{f2<-gl(2,10);y2<-rnorm(20);c(lbf_anova(y2,f2,.25),lbf_F(Fstat(y2,10,2),10,2,.25))},"\n")
print(system.time(for(i in 1:2000){y<-rnorm(30)+rep(ralpha(3),each=10);l<-lbf_F(Fstat(y,10,3),10,3)}))
print(system.time(for(i in 1:2000){y<-rnorm(20)+rep(ralpha(2),each=10);l<-lbf_F(Fstat(y,10,2),10,2)}))
cat("replis :",.rep$n,"\n")
