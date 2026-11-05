source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/agent_2x_commun.R")
set.seed(20261105)
cat("BayesFactor", as.character(packageVersion("BayesFactor")), "R", R.version.string, "\n")
# controle du codage
for (a in 2:4) { P <- as.matrix(proj(a)); cat("a=",a," ||P'P-I||=",max(abs(crossprod(P)-diag(a-1))),
   " ||PP'-(I-1/a)||=",max(abs(tcrossprod(P)-(diag(a)-1/a))),"\n") }
d <- data.frame(y=rnorm(30),A=gl(3,10)); dt <- c(A="fixed")
X <- as.matrix(BayesFactor:::fullDesignMatrix(y~A,d,dt)); X2 <- as.matrix(model.matrix(~A-1,d)) %*% as.matrix(proj(3))
cat("X paquet = indicatrices x P :", max(abs(X-X2)), "\n")
d2 <- expand.grid(rep=1:2,A=gl(2,1),B=gl(3,1)); dt2 <- c(A="fixed",B="fixed")
X <- as.matrix(BayesFactor:::oneDesignMatrix("A:B",d2,dt2))
XA <- as.matrix(model.matrix(~A-1,d2))%*%as.matrix(proj(2)); XB <- as.matrix(model.matrix(~B-1,d2))%*%as.matrix(proj(3))
Xk <- do.call(cbind,lapply(1:2,function(j) XA[,1]*XB[,j]))
cat("X_AB paquet vs produit ligne a ligne :", max(abs(X-Xk)), "\n")
# loi de |b| quand q=1 : mediane = r (Cauchy d'echelle r)
b <- replicate(1e5, rnorm(1,0,sqrt(rg(.5)))); cat("mediane |b| (r=.5, attendu .5) :", median(abs(b)), "\n")
# equivalences de calcul
d3 <- data.frame(y=rnorm(20)+rep(c(0,1),10),A=gl(2,1,20))
cat("ttestBF vs anovaBF k=2 :", lbf_ttest(d3$y,d3$A), lbf_anova(d3$y,d3$A), "\n")
tm <- function(e, R=100) { t0<-Sys.time(); for(i in 1:R) eval(e, parent.frame()); as.numeric(Sys.time()-t0,units="secs")/R*1000 }
for (n in c(10,30,100)) {
  f3<-gl(3,n); f2<-gl(2,n)
  cat("n=",n," anovaBF k=3 ms:",tm(quote(lbf_anova(rnorm(3*n),f3))),
      " ttestBF k=2 ms:",tm(quote(lbf_ttest(rnorm(2*n),f2))),
      " aov ms:",tm(quote(paov(rnorm(3*n),f3))),"\n")
}
d4 <- expand.grid(i=1:10,A=gl(2,1),B=gl(3,1))
t0<-Sys.time(); for(i in 1:10){d4$y<-rnorm(60); o<-anovaBF(y~A*B,data=d4,whichModels="withmain",progress=FALSE)}
cat("A*B withmain s/jeu :",as.numeric(Sys.time()-t0,units="secs")/10,"\n")
print(o@bayesFactor[,c("bf","error")])
print(rownames(o@bayesFactor))
cat("cauchy extreme : ", {y<-rnorm(30)+rep(c(-3000,0,3000),each=10); lbf_anova(y,gl(3,10))},
    {y<-rnorm(20)+rep(c(-3000,3000),each=10); lbf_ttest(y,gl(2,10))},"\n")
print(system.time({ for(i in 1:200) o2<-tryCatch(lbf_anova(rnorm(30)+rep(c(-50,0,50),each=10),gl(3,10)),error=function(e) NA)}))
