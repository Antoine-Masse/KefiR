.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); options(BFprogress=FALSE)
digam <- function(g,r) (r^2/2)^.5/gamma(.5)*g^(-1.5)*exp(-r^2/(2*g))
Qs <- function(a) { e <- eigen(diag(a)-1/a); e$vectors[,seq_len(a-1),drop=FALSE] }
# posterieur d'un facteur equilibre sous H1 ; renvoie BF, P(f<w|H1)
post1 <- function(means, s2w, a, nc, r=.5, w=.05, D=4000){ # means : moyennes de groupe, s2w : SSW
  n<-a*nc; q<-a-1; Q<-Qs(a); mu<-mean(means); H<-nc*sum((means-mu)^2); T<-H+s2w
  bh<-drop(t(Q)%*%(means-mu)) # beta chapeau (X'X=nc I)
  u<-seq(-14,14,length.out=300); g<-exp(u); wg<-nc*g/(1+nc*g)
  lw<- -q/2*log(1+nc*g)-(n-1)/2*log(1-wg*H/T)
  dens<-exp(lw)*digam(g,r)*g; bf<-sum(dens)*(u[2]-u[1])
  idx<-sample(length(g),D,TRUE,prob=dens); gd<-g[idx]; wd<-wg[idx]
  S<-T-wd*H; s2<-S/rchisq(D,n-1)
  bd<-matrix(rnorm(D*q),D,q)*sqrt(s2/(nc+1/gd))+outer(wd,bh)
  f<-sqrt(rowSums(bd^2)/a/s2)
  c(bf=bf,inR=mean(f<w),medf=median(f))
}
cat("=== 0. Validation contre posterior() de BayesFactor, a=3, nc=30\n")
set.seed(21); a<-3;nc<-30;ff<-gl(a,nc);y<-rnorm(a*nc)+.2*(as.integer(ff)-2); dd<-data.frame(y,ff)
m<-tapply(y,ff,mean); ssw<-sum((y-ave(y,ff))^2); pr<-post1(m,ssw,a,nc)
b<-anovaBF(y~ff,dd); ch<-posterior(b,iterations=20000); eff<-as.matrix(ch)[,2:4]; fch<-sqrt(rowMeans(eff^2)/as.matrix(ch)[,"sig2"])
cat(sprintf(" BF %.4f / %.4f ; P(f<.05|H1) %.4f / %.4f ; med f %.4f / %.4f\n",pr["bf"],extractBF(b)$bf,pr["inR"],mean(fch<.05),pr["medf"],median(fch)))

cat("\n=== 1. ROPE : verite exactement nulle (f=0), puis f=.03 (negligeable), puis f=.2 (hors ROPE). 40 jeux/cellule, medianes\n")
cat(" a   n/gr  f_vrai  P(H1|D)  P(f<.05|D,H1)  P(ROPE|D)=P0+P1*.   med f post | H1\n")
for(fv in c(0,.03,.2)) for(a in c(2,3,5,10)) for(nc in c(30,100,500,2000,10000)){
  set.seed(a*nc+round(fv*100))
  z<-t(replicate(40,{ al<-rep(c(-1,1),length.out=a); if(a%%2==1) al[a]<-0; al<-al-mean(al); al<-al*fv/sqrt(mean(al^2)+(fv==0));
     mm<-al+rnorm(a,0,1/sqrt(nc)); ssw<-rchisq(1,a*(nc-1)); p<-post1(mm,ssw,a,nc,D=1500); P1<-p["bf"]/(1+p["bf"]); c(P1,p["inR"],(1-P1)+P1*p["inR"],p["medf"])}))
  cat(sprintf(" %-3d %-5d %-6.2f  %-7.3f  %-13.4f  %-9.3f  %-7.3f\n",a,nc,fv,median(z[,1]),median(z[,2]),median(z[,3]),median(z[,4])))
}
cat("\n=== 2. Biais du posterieur de f : moyenne attendue sous H1 pour verite nulle, a=5\n")
for(nc in c(10,30,100,500)) {set.seed(nc); z<-replicate(40,{mm<-rnorm(5,0,1/sqrt(nc));p<-post1(mm,rchisq(1,5*(nc-1)),5,nc,D=1500);p["medf"]}); cat(sprintf(" n/gr=%d med f|H1=%.3f ; sqrt(q/(a*nc))=%.3f\n",nc,median(z),sqrt(4/5/nc)))}

cat("\n=== 3. Un lot deviant parmi a lots (r=.5), n/lot=10 ; 200 jeux, medianes\n")
BFdev<-function(a,nc,D,r=.5,R=200){out<-t(replicate(R,{mm<-c(D,rep(0,a-1))+rnorm(a,0,1/sqrt(nc));ssw<-rchisq(1,a*(nc-1));H<-nc*sum((mm-mean(mm))^2);n<-a*nc;q<-a-1;
   u<-seq(-14,14,length.out=300);g<-exp(u);lw<- -q/2*log(1+nc*g)-(n-1)/2*log(1-nc*g/(1+nc*g)*H/(H+ssw));bf<-sum(exp(lw)*digam(g,r)*g)*(u[2]-u[1]);c(P=bf/(1+bf),p=1-pf((H/q)/(ssw/(n-a)),q,n-a))}));c(median(out[,1]),median(out[,2]),mean(out[,1]<.5))}
set.seed(31); for(a in c(3,5,10,20,40)) for(D in c(1.5,2,3)){v<-BFdev(a,10,D); cat(sprintf(" a=%-2d Delta=%.1f sigma : P(H1) med=%.3f ; p aov med=%.2g ; frac 'plutot non' (P<.5)=%.2f\n",a,D,v[1],v[2],v[3]))}
cat(" r=.25 / r=1 a a=20, Delta=2 :",round(BFdev(20,10,2,r=.25)[1],3),round(BFdev(20,10,2,r=1)[1],3),"\n")

cat("\n=== 4. Non-normalite ; une branche unique gaussienne (niveau 1)\n")
Pa<-function(y,f){b<-extractBF(anovaBF(y~f,data.frame(y,f),progress=FALSE))$bf;c(P=b/(1+b),p=summary(aov(y~f))[[1]][1,5],kw=kruskal.test(y,f)$p.value)}
set.seed(41); f<-gl(3,15)
y<-rnorm(45)+1*(as.integer(f)-2); cat(" Sans aberrant, moyennes -1,0,1 :",round(Pa(y,f),4),"\n"); y2<-y; y2[1]<-30; cat(" Une valeur aberrante +30 (groupe 1) :",round(Pa(y2,f),4),"\n")
rt<-function(g,R=200) {z<-t(replicate(R,g())); c(oui=mean(z[,1]>.5),fort=mean(z[,1]>.95),non=mean(z[,1]<.25),medP=median(z[,1]),p05=mean(z[,2]<.05),kw05=mean(z[,3]<.05))}
set.seed(42)
print(round(rbind(
 lognorm_effet_raw=rt(function(){y<-exp(rnorm(45,0,1.5)+.6*(as.integer(f)-2));Pa(y,f)}),
 lognorm_effet_log=rt(function(){y<-rnorm(45,0,1.5)+.6*(as.integer(f)-2);Pa(y,f)}),
 t2_H0=rt(function(){y<-rt(45,2);Pa(y,f)})),3))
f2<-factor(rep(1:3,c(5,5,40))); set.seed(43)
print(round(rbind(hetero_H0_petit_grand_sd=rt(function(){y<-rnorm(50,0,c(4,4,1)[f2]);Pa(y,f2)}),
  hetero_H0_inverse=rt(function(){y<-rnorm(50,0,c(1,1,4)[f2]);Pa(y,f2)})),3))

cat("\n=== 5. Sphericite violee : T=6, n=20, AR1(.95) ; H0 vraie ; BF intercept aleatoire\n")
set.seed(51); nT<-6;ns<-20;Sig<-.95^abs(outer(1:nT,1:nT,"-"));L<-chol(Sig)
z<-t(replicate(120,{Y<-matrix(rnorm(ns*nT),ns)%*%L; dl<-data.frame(y=c(Y),T=factor(rep(1:nT,each=ns)),id=factor(rep(1:ns,nT)))
  b<-extractBF(lmBF(y~T+id,dl,whichRandom="id",iterations=2000)/lmBF(y~id,dl,whichRandom="id",iterations=2000))$bf
  a<-summary(aov(y~T+Error(id/T),dl))[[1]][[1]]; c(P=b/(1+b),p=a[1,5])}))
cat(sprintf(" P>.5 : %.3f ; P>.95 : %.3f ; p<.05 non corrige : %.3f ; med P : %.3f (120 rep)\n",mean(z[,1]>.5),mean(z[,1]>.95),mean(z[,2]<.05),median(z[,1])))
cat("\nFIN\n")
