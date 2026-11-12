.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); options(BFprogress=FALSE)
Qs <- function(a) { e <- eigen(diag(a)-1/a); e$vectors[,seq_len(a-1),drop=FALSE] }
digam <- function(g,r) (r^2/2)^.5/gamma(.5)*g^(-1.5)*exp(-r^2/(2*g))
u<-seq(-14,16,length.out=400); gg<-exp(u); du<-u[2]-u[1]
ld <- function(M) as.numeric(determinant(M)$modulus)
Xf <- function(f){f<-factor(f); scale(model.matrix(~0+f)%*%Qs(nlevels(f)),TRUE,FALSE)}
# BF MANOVA generale (un terme) par valeurs propres de X'X et algebre
lBFg <- function(Yc,X,g){d<-ncol(Yc);n<-nrow(Yc);q<-ncol(X);A<-diag(1/g,q)+crossprod(X);XY<-crossprod(X,Yc);S<-crossprod(Yc)-t(XY)%*%solve(A,XY); (d/2)*(-ld(A)-q*log(g))-(n-1)/2*(ld(S)-ld(crossprod(Yc)))}
BFman <- function(Y,X,r=.5){Yc<-scale(Y,TRUE,FALSE); gv<-gg[seq(1,400,by=4)]; v<-sapply(gv,function(g) lBFg(Yc,X,g)); mx<-max(v+log(digam(gv,r)*gv)); mx+log(sum(exp(v+log(digam(gv,r)*gv)-mx))*du*4)} # log BF
BFbal <- function(Y,f,r=.5){ # equilibre : rapide par valeurs propres
  Yc<-scale(Y,TRUE,FALSE); n<-nrow(Y);d<-ncol(Y);a<-nlevels(f);nc<-n/a;q<-a-1
  mm<-apply(Yc,2,function(z) tapply(z,f,mean)); H<-nc*crossprod(mm); E0<-crossprod(Yc); th<-pmin(Re(eigen(solve(E0,H),only.values=TRUE)$values),1-1e-14); th<-th[th>1e-14]
  w<-nc*gg/(1+nc*gg); lw<- -q*d/2*log(1+nc*gg)-(n-1)/2*rowSums(sapply(th,function(t) log1p(-w*t))); v<-lw+log(digam(gg,r)*gg); mx<-max(v); mx+log(sum(exp(v-mx))*du)}
pil<-function(Y,f) if(ncol(Y)==1) summary(aov(Y[,1]~f))[[1]][1,5] else summary(manova(Y~f),test="Pillai")$stats[1,6]
cat("=== 0. BFbal vs BFman :"); set.seed(1); f<-gl(3,10);Y<-matrix(rnorm(90),30); cat(BFbal(Y,f),BFman(Y,Xf(f)),"\n")

cat("\n=== 1. Calibration sous l'a priori exact (Sigma=I par invariance), a=3 ; pi=1/2\n")
calib<-function(d,nc,R=2000,r=.5){a<-3;q<-2;n<-a*nc;f<-gl(a,nc);Q<-Qs(a);res<-matrix(NA,R,2)
 for(i in 1:R){h<-rbinom(1,1,.5);Y<-matrix(rnorm(n*d),n); if(h){g<-r^2/rchisq(1,1);B<-matrix(rnorm(q*d,0,sqrt(g)),q,d);Y<-Y+(Q%*%B)[as.integer(f),,drop=FALSE]}
  res[i,]<-c(plogis(BFbal(Y,f,r)),h)}
 b<-cut(res[,1],c(0,.1,.3,.5,.7,.9,1),include.lowest=TRUE);cat(sprintf(" d=%d n/gr=%d : moyenne P=%.3f (se %.3f)\n",d,nc,mean(res[,1]),sd(res[,1])/sqrt(R)));print(round(rbind(annonce=tapply(res[,1],b,mean),realise=tapply(res[,2],b,mean)),3))}


cat("\n=== 2. Effet sur une seule DV (moyennes -.5,0,.5 ; f=.408), DV ajoutees sans rapport, 300 jeux\n")
cat(" d  n/gr | H1: med P, P>.5, P>.95, P<.25, Pillai<.05 | H0: med P, P>.5, P<.05, Pillai<.05 | med Mahalanobis f a priori\n")
for(nc in c(15,30)) for(d in c(1,2,3,5,10,20)){ if(3*nc-3<d+3) next; f<-gl(3,nc); n<-3*nc; set.seed(100+d+nc)
 s<-function(eff){t(replicate(300,{Y<-matrix(rnorm(n*d),n);Y[,1]<-Y[,1]+eff*(as.integer(f)-2);c(plogis(BFbal(Y,f)),pil(Y,f))}))}
 z1<-s(.5);z0<-s(0); medM<-median(sqrt(.25/rchisq(1e5,1)*rchisq(1e5,2*d)/3))
 cat(sprintf(" %-2d %-3d | %.3f %.3f %.3f %.3f %.3f | %.3f %.3f %.3f %.3f | %.2f\n",d,nc,median(z1[,1]),mean(z1[,1]>.5),mean(z1[,1]>.95),mean(z1[,1]<.25),mean(z1[,2]<.05),median(z0[,1]),mean(z0[,1]>.5),mean(z0[,1]<.05),mean(z0[,2]<.05),medM))}

cat("\n=== 3. Sigma non commune sous H0, desequilibre (10,10,40) ; Sigma groupe 1,2 = 9*Sigma3 ; d=3 ; 300 jeux\n")
fu<-factor(rep(1:3,c(10,10,40)));Xu<-Xf(fu);set.seed(3)
for(sc in list(c(3,3,1),c(1,1,3),c(1,1,1))){z<-t(replicate(300,{Y<-matrix(rnorm(60*3),60)*sc[fu];c(plogis(BFman(Y,Xu)),pil(Y,fu))}));cat(sprintf(" sd=%s : P>.5 %.3f ; P>.95 %.3f ; Pillai p<.05 %.3f\n",paste(sc,collapse=","),mean(z[,1]>.5),mean(z[,1]>.95),mean(z[,2]<.05)))}

cat("\n=== 4. Sensibilite r a d=3 et d=10, H1 ci-dessus (un jeu moyen) ; n/gr=15\n")
set.seed(4);f<-gl(3,15);for(d in c(3,10)){Y<-matrix(rnorm(45*d),45);Y[,1]<-Y[,1]+.8*(as.integer(f)-2);cat(sprintf(" d=%d p Pillai=%.2g ; P(r=.25,.5,1)=%s\n",d,pil(Y,f),paste(round(plogis(sapply(c(.25,.5,1),function(r) BFbal(Y,f,r))),3),collapse=" ")))}

cat("\n=== 5. Partitions\n")
parts<-function(k){ if(k==1) return(list(1L)); out<-list(); for(p in parts(k-1)){ m<-max(p); for(j in 1:(m+1)) out[[length(out)+1]]<-c(p,j)}; out}
lBF1<-function(y,cl){ if(length(unique(cl))==1) return(0); X<-Xf(cl); yc<-y-mean(y); n<-length(y); XtX<-crossprod(X); Xty<-crossprod(X,yc); q<-ncol(X); T<-sum(yc^2); ev<-eigen(XtX,symmetric=TRUE); bb<-crossprod(ev$vectors,Xty)
  gv<-gg[seq(1,400,by=5)]; lv<-sapply(gv,function(g){ -.5*sum(log1p(g*ev$values))-(n-1)/2*log1p(-sum(bb^2/(ev$values+1/g))/T)}); v<-lv+log(digam(gv,.5)*gv); mx<-max(v); mx+log(sum(exp(v-mx))*du*5)}
analyse<-function(y,f){k<-nlevels(f);P<-parts(k);lb<-sapply(P,function(p) lBF1(y,p[as.integer(f)]));nul<-which(sapply(P,function(p) length(unique(p))==1));pr<-rep(.5/(length(P)-1),length(P));pr[nul]<-.5;w<-pr*exp(lb);w<-w/sum(w)
  full<-which(sapply(P,function(p) length(unique(p))==k)); pair<-matrix(NA,k,k);for(i in 1:k)for(j in 1:k) pair[i,j]<-sum(w[sapply(P,function(p) p[i]==p[j])])
  list(Pglob_part=1-w[nul],Pglob_terme=plogis(lb[full]),pair=pair,w=w,P=P)}
set.seed(5);f<-gl(4,15);y<-rnorm(60)+c(0,0,0,.9)[f];a<-analyse(y,f);cat(sprintf(" k=4, 1 groupe devie .9 : P(terme global, Sec.5.2)=%.3f ; 1-P(toutes egales, Sec.9)=%.3f\n",a$Pglob_terme,a$Pglob_part));print(round(a$pair,3))
cat(" Vs lmBF full :",round(extractBF(anovaBF(y~f,data.frame(y,f)))$bf/(1+extractBF(anovaBF(y~f,data.frame(y,f)))$bf),3),"\n")
glob<-t(replicate(150,{y<-rnorm(60)+c(0,0,0,.6)[f];a<-analyse(y,f);c(a$Pglob_terme,a$Pglob_part)}));cat(sprintf(" Sur 150 jeux (ecart .6): med P terme=%.3f ; med P partition=%.3f ; ecart abs med=%.3f ; ecart max=%.3f ; desaccord autour de .5 : %.3f\n",median(glob[,1]),median(glob[,2]),median(abs(glob[,1]-glob[,2])),max(abs(glob[,1]-glob[,2])),mean((glob[,1]>.5)!=(glob[,2]>.5))))
cat("\n Multiplicite sur paires nulles, k=5 : groupes 1-4 egaux, groupe 5 decale 1.5 ; n=20 ; 150 jeux\n")
f5<-gl(5,20);set.seed(6);z<-t(replicate(150,{y<-rnorm(100)+c(0,0,0,0,1.5)[f5];a<-analyse(y,f5);pn<-1-a$pair[1:4,1:4][upper.tri(diag(4))];c(max(pn),mean(pn))}))
cat(sprintf(" P(diff) a priori conditionnelle ; max sur 6 paires nulles : mediane %.3f ; >.5 : %.3f ; >.75 : %.3f ; moyenne %.3f\n",median(z[,1]),mean(z[,1]>.5),mean(z[,1]>.75),mean(z[,2])))
f5<-gl(5,20);z<-t(replicate(150,{y<-rnorm(100);a<-analyse(y,f5);pn<-1-a$pair[upper.tri(a$pair)];c(max(pn),a$Pglob_part)}))
cat(sprintf(" H0 globale, k=5, n=20 : med max P(diff)=%.3f ; >.5 : %.3f ; P(toutes pas egales)>.5 : %.3f\n",median(z[,1]),mean(z[,1]>.5),mean(z[,2]>.5)))
cat("\n Saut A/B : meme paire (groupes 1 et 2, ecart .8, n=20), autres groupes au centre\n")
set.seed(7); base<-function(k){f<-gl(k,20);y<-rnorm(20*k);y<-y-ave(y,f);y<-y+c(-.4,.4,rep(0,k-2))[f];list(y=y,f=f)}
for(k in c(3,4,5)){s<-base(k);a<-analyse(s$y,s$f);cat(sprintf(" k=%d partition P(mu1!=mu2)=%.3f\n",k,1-a$pair[1,2]))}
s<-base(6);i<-s$f%in%1:2;bfp<-extractBF(anovaBF(y~f,data.frame(y=s$y[i],f=droplevels(s$f[i]))))$bf;for(k in 6:8){m<-k*(k-1)/2;pi<-.5^(1/m);o<-(1-pi)/pi;cat(sprintf(" k=%d option B P(mu1!=mu2)=%.3f (BF paire=%.2f, cotes %.3f)\n",k,o*bfp/(1+o*bfp),bfp,o))}
cat("\nFIN\n")
