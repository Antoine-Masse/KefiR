source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
set.seed(1); mk <- function(n,m,s){y<-unlist(Map(function(a,b,c){v<-rnorm(a);b+c*v},n,m,s));list(y=y,g=rep(seq_along(n),n))}
for (cfg in list(list(k=20,sp=3,M=3),list(k=20,sp=1,M=0),list(k=10,sp=3,M=3),list(k=5,sp=3,M=3))) {
 k<-cfg$k; d <- mk(rep(10,k),seq(0,cfg$M,length.out=k),exp(seq(-cfg$sp,cfg$sp,length.out=k))); cat("k=",k," log sd +-",cfg$sp,"\n")
 res <- t(sapply(1:6,function(s){tb<-tab_of(d$y,d$g,draws=4000L,seed=s,laws="normal_het"); c(tb$log_ml,tb$rel_se,tb$ess)}))
 colnames(res)<-c("H0","H1","se0","se1","ess0","ess1"); print(round(res,3)); cat("sd reelle H0",round(sd(res[,1]),3)," H1",round(sd(res[,2]),3),"; logBF",round(res[,2]-res[,1],2),"\n"); flush.console()
}
