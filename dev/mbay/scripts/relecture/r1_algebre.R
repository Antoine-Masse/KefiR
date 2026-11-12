.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); options(BFprogress=FALSE)
Qs <- function(a) { e <- eigen(diag(a)-1/a); e$vectors[,seq_len(a-1),drop=FALSE] }
Xf <- function(f){f<-factor(f); scale(model.matrix(~0+f)%*%Qs(nlevels(f)),TRUE,FALSE)}
ld <- function(M) as.numeric(determinant(M)$modulus)
logBFg <- function(Y,X,g){ # §6 : log BF(g) a g commun
  Y<-scale(as.matrix(Y),TRUE,FALSE); n<-nrow(Y); d<-ncol(Y); q<-ncol(X)
  A<-diag(1/g,q)+crossprod(X); S<-crossprod(Y)-t(crossprod(X,Y))%*%solve(A,crossprod(X,Y))
  (d/2)*(-ld(A)-q*log(g)) - (n-1)/2*(ld(S)-ld(crossprod(Y)))
}
digam <- function(g,r) (r^2/2)^.5/gamma(.5)*g^(-1.5)*exp(-r^2/(2*g))
BFint <- function(Y,X,r=.5) integrate(Vectorize(function(u){g<-exp(u);exp(logBFg(Y,X,g))*digam(g,r)*g}),-25,25,rel.tol=1e-7,subdivisions=500)$value

cat("=== 1. Derivation independante BF(g) par esperance sous posterieur de M0 (d=2, n=12)\n")
set.seed(1); n<-12; d<-2; f<-gl(3,4); X<-Xf(f); Y<-matrix(rnorm(n*d),n)%*%chol(matrix(c(1,.5,.5,1),2)); Y[,1]<-Y[,1]+.8*(as.integer(f)-2)
Yc<-scale(Y,TRUE,FALSE); g<-.7
# BF(g)=E_{Sigma|Y,M0}[ p(Y|Sigma,g,M1)/p(Y|Sigma,M0) ], Sigma|Y,M0 ~ IW(n-1, Y'Y)
M<-4e5; W<-rWishart(M,n-1,solve(crossprod(Yc)))
Pn<-diag(n)-1/n; Vm<-X%*%t(X)*g   # cov ligne sous M1 dans l'espace contraste : I+gXX'
iV<-solve(diag(n)+Vm); ldV<-ld(diag(n)+Vm)
T1<-t(Yc)%*%iV%*%Yc; T0<-crossprod(Yc)
ll<-sapply(1:M,function(i) -(d/2)*ldV-.5*sum(W[,,i]*(T1-T0)))
est<-log(mean(exp(ll-max(ll))))+max(ll); w<-exp(ll-max(ll))
cat(sprintf(" forme fermee %.5f ; MC independant %.5f (se %.4f)\n",logBFg(Y,X,g),est,sd(w)/sqrt(M)/mean(w)))

cat("\n=== 2. Forme equilibree, valeurs propres\n")
E0<-crossprod(Yc); H<-t(crossprod(X,Yc))%*%solve(crossprod(X),crossprod(X,Yc)); th<-Re(eigen(solve(E0,H))$values); nc<-4; wg<-nc*g/(1+nc*g)
cat(sprintf(" formule %.6f ; forme propre %.6f\n",logBFg(Y,X,g), -2*d/2*log(1+nc*g)-(n-1)/2*sum(log(1-wg*th))))

cat("\n=== 3. Reduction d=1 contre anovaBF\n")
set.seed(2); ff<-factor(rep(1:3,c(8,12,20))); y<-rnorm(40)+.4*(as.integer(ff)-2)
for(r in c(.25,.5,1)) cat(sprintf(" r=%.2f integre=%.5f anovaBF=%.5f\n",r,log(BFint(y,Xf(ff),r)),extractBF(anovaBF(y~ff,data.frame(y,ff),rscaleFixed=r))$bf|>log()))

cat("\n=== 4. Sec.6 pt 4 : 'n-1-Sq<d : l'integrale peut diverger' ?\n")
cat(" Exposant asymptotique g->inf : g^(-3/2) * g^(-(q-k)(d-k)/2), k=q+d-(n-1) valeurs theta=1\n")
for(cfg in list(c(n=6,q=2,d=4),c(n=6,q=2,d=5),c(n=8,q=3,d=5),c(n=9,q=2,d=8),c(n=4,q=3,d=1),c(n=7,q=3,d=6))){
  nn<-cfg["n"];q<-cfg["q"];dd<-cfg["d"];a<-q+1;cn<-nn%/%a; ff<-factor(rep(1:a,length.out=nn))
  set.seed(5);Y0<-matrix(rnorm(nn*dd),nn); X0<-Xf(ff)
  v<-sapply(c(1e2,1e4,1e6,1e8),function(g) logBFg(Y0,X0,g))
  k<-q+dd-(nn-1)
  I<-tryCatch(BFint(Y0,X0,.5),error=function(e)NA)
  cat(sprintf(" n=%d q=%d d=%d k=%d : logBF(g=1e2..1e8)=%s ; int=%.4g\n",nn,q,dd,k,paste(round(v,2),collapse=" "),I))
}

cat("\n=== 5. Masses a priori : f^2 = (r^2 q/a) F(q,1)\n")
cat(sprintf(" %-3s %-8s %-10s %-10s %-10s %-10s\n","a","P(f<.05)","med f","P(f>.4)","P(f>1)","P(ROPE)π=.5"))
for(a in c(2,3,4,5,6,10,20)){q<-a-1;r<-.5; F<-function(x) pf(x^2*a/(r^2*q),q,1)
 medf<-sqrt(r^2*q/a*qf(.5,q,1)); cat(sprintf(" %-3d %-8.4f %-10.3f %-10.3f %-10.3f %-10.3f\n",a,F(.05),medf,1-F(.4),1-F(1),.5+.5*F(.05)))}
cat(" verif MC a=3 : "); set.seed(3); g<-.25/rchisq(1e6,1); b<-matrix(rnorm(2e6),ncol=2)*sqrt(g); fm<-sqrt(rowSums(b^2)/3); cat(mean(fm<.05),median(fm),"\n")
cat(" Paire : var(alpha_i-alpha_j)/g pour a=3,5,10 (doit dependre de a si Sec.5.2 a raison) : ")
for(a in c(3,5,10)){Q<-Qs(a); cat(round(sum((Q[1,]-Q[2,])^2),4)," ")}; cat("\n")

cat("\n=== 6. Partitions et Bell ; saut Option A -> B\n")
B<-c(1,1,2,5,15,52,203,877,4140,21147); bell<-function(k) B[k+1]
for(k in 3:9){ pa<-.5+.5*(bell(k-1)-1)/(bell(k)-1); m<-k*(k-1)/2; pb<-.5^(1/m); pt<-.5^(1/(k-1))
 cat(sprintf(" k=%d P_partition(egales)=%.3f ; pi0^(1/m)=%.3f ; pi0^(1/(k-1))=%.3f ; cotes H1:H0 %.3f / %.3f\n",k,pa,pb,pt,(1-pa)/pa,(1-pb)/pb))}
cat(" Independance de 3 paires : P(1=2,2=3)=pi^2 et P(1=2,2=3,1=3)=pi^3 alors que le 3e evenement est implique : impossible si pi<1\n")
cat(" k=5 -> 6 facteur sur cotes a priori :", ((1-.637)/.637)/((1-.5^(1/15))/.5^(1/15)),"\n")

cat("\n=== 7. Plafond P(H0|D) a F=0 et correspondance p<->P, un facteur 3 groupes r=.5\n")
BFF <- function(Fv,nc,a=3,r=.5){ # balanced
 q<-a-1;n<-a*nc;R2<-Fv*q/(Fv*q+n-a)
 integrate(function(u){g<-exp(u);exp(-q/2*log(1+nc*g)-(n-1)/2*log(1-nc*g/(1+nc*g)*R2))*digam(g,r)*g},-30,30,subdivisions=1000)$value}
for(nc in c(3,5,10,15,30,100,300,1000)){ fa<-qf(.95,2,3*nc-3); bf0<-BFF(0,nc); bf5<-BFF(fa,nc); fb<-uniroot(function(x) log(BFF(x,nc)),c(.01,40))$root
 cat(sprintf(" n/gr=%4d : F=0 => BF=%.3f, P(H0)max=%.3f ; p=.05 => P(H1)=%.3f ; P(H1)=.5 <=> p=%.4f\n",nc,bf0,1/(1+bf0),bf5/(1+bf5),1-pf(fb,2,3*nc-3)))}

cat("\n=== 8. Effet de l'unite si la covariable est 'centree' seulement (Sec.2 litteral)\n")
set.seed(4); n<-40;x<-rnorm(n);y<-.4*x+rnorm(n)
lit<-function(x,y,r=sqrt(2)/4) BFint(y,scale(cbind(x),TRUE,FALSE),r)
cat(sprintf(" Sec.2 litteral, x en unite 1 : BF=%.3f ; x*1000 : BF=%.3f ; x/1000 : BF=%.3f ; lmBF (unite 1 / *1000) : %.3f / %.3f\n",lit(x,y),lit(x*1000,y),lit(x/1000,y),extractBF(lmBF(y~x,data.frame(x,y)))$bf,extractBF(lmBF(y~x,data.frame(x=x*1000,y)))$bf))
cat(" ROPE : facteur 2 niveaux code +/-.5 : sd(x)=.5 ; ROPE covariable |beta*.5/sigma|<.1 <=> |delta|<.2 ; ROPE facteur |delta|<.1\n")

cat("\n=== 9. f^2 Sec.4 vs Sec.6 a d=1, a=2 : delta=1\n")
cat(" Sec.4 f=",.5," ; Sec.6 f^2=tr(S^-1 B'B)/q avec B=beta*=delta/sqrt2 => f=",1/sqrt(2),"\n")
cat("\nFIN\n")
