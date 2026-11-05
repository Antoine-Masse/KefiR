.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); options(BFprogress=FALSE)
B <- function(e) extractBF(e)$bf
cat("=== A. Croisement pur 2x2, n=10/cellule, A sans effet moyen, A:B fort\n")
set.seed(11); d<-expand.grid(A=gl(2,1),B=gl(2,1),r=1:10); s<-ifelse(d$A==1,1,-1)*ifelse(d$B==1,1,-1)
d$y<-1.0*s+rnorm(40); d$y<-d$y-ave(d$y,d$A)+mean(d$y)+0  # moyennes marginales de A strictement egales
it<-2e5
m<-list(A=B(lmBF(y~A,d,iterations=it)),B=B(lmBF(y~B,d,iterations=it)),AB=B(lmBF(y~A+B,d,iterations=it)),F=B(lmBF(y~A+B+A:B,d,iterations=it)))
bfA<-(m$A+m$AB)/(1+m$B); bfI<-m$F/m$AB
cat(sprintf(" modeles : A=%.3g B=%.3g A+B=%.3g plein=%.4g\n",m$A,m$B,m$AB,m$F))
cat(sprintf(" Apparie : P(A)=%.3f  P(A:B)=%.5f  => P(A:B)>P(A) ; 'Non pour A avec p=%.3f'\n",bfA/(1+bfA),bfI/(1+bfI),1/(1+bfA)))
tot<-1+m$A+m$B+m$AB+m$F; cat(sprintf(" Uniforme 5 modeles : P(A)=%.5f ; Chipman (1/4,1/4,1/4,1/8,1/8): P(A)=%.5f\n",(m$A+m$AB+m$F)/tot,(m$A/4+m$AB/8+m$F/8)/(1/4+m$A/4+m$B/4+m$AB/8+m$F/8)))
print(anova(lm(y~A*B,d))); print(sapply(split(d,d$B),function(x) t.test(y~A,x,var.equal=TRUE)$p.value))

cat("\n=== B. Apparie = conditionnel de la loi uniforme rejetee (identite)\n")
cat(sprintf(" bfA/(1+bfA)=%.6f ; [m(A)+m(A+B)]/[1+m(A)+m(B)+m(A+B)]=%.6f\n",bfA/(1+bfA),(m$A+m$AB)/(1+m$A+m$B+m$AB)))
cat(sprintf(" Type II exact A+B vs B : BF=%.3f ; apparie : %.3f (agent_02 : 6.4449/0.31457=%.2f contre 18.7)\n",m$AB/m$B,bfA,6.4448757/0.31457))

cat("\n=== C. Calibration du P(A) apparie sous la loi jointe (2x2, n=10/cell, quadrature)\n")
N<-40; u<-seq(-9,14,length.out=26); gs<-exp(u); r<-.5
wt<-(r^2/2)^.5/gamma(.5)*gs^(-1.5)*exp(-r^2/(2*gs))*gs*(u[2]-u[1]); wt<-wt/sum(wt)
nA<-N/2; nI<-N/4
wA<-nA*gs/(1+nA*gs); wI<-nI*gs/(1+nI*gs); pA<-wt*(1+nA*gs)^(-.5); pI<-wt*(1+nI*gs)^(-.5)
evid <- function(sA,sB,sI){ # sX = SS_X/SST
  e<- -(N-1)/2
  mA<-sum(pA*(1-wA*sA)^e); mB<-sum(pA*(1-wA*sB)^e)
  GA<-outer(wA,wA,"+"); PA<-outer(pA,pA)
  mAB<-sum(PA*(1-(outer(wA*sA,wA*sB,"+")))^e)
  M3<-0; for(k in seq_along(wI)) M3<-M3+pI[k]*sum(PA*(1-outer(wA*sA,wA*sB,"+")-wI[k]*sI)^e)
  c(n=1,A=mA,B=mB,AB=mAB,F=M3)}
# validation contre lmBF
ss<-function(y,xA,xB,xI){yc<-y-mean(y);S<-sum(yc^2);c(sum(yc*xA)^2/sum(xA^2)/S,sum(yc*xB)^2/sum(xB^2)/S,sum(yc*xI)^2/sum(xI^2)/S)}
x<-expand.grid(a=c(-1,1),b=c(-1,1),r=1:10); xA<-x$a/sqrt(2);xB<-x$b/sqrt(2);xI<-x$a*x$b/2
set.seed(12); yy<-.5*x$a+.3*x$a*x$b+rnorm(40); s0<-ss(yy,xA,xB,xI); ev<-evid(s0[1],s0[2],s0[3]); dd<-data.frame(y=yy,A=factor(x$a),B=factor(x$b))
cat(" validation quadrature / lmBF : ",round(c(ev[2],B(lmBF(y~A,dd)),ev[4],B(lmBF(y~A+B,dd,iterations=1e5)),ev[5],B(lmBF(y~A*B,dd,iterations=1e5))),3),"\n")
simC<-function(pm,R=3000){res<-matrix(NA,R,6);for(i in 1:R){k<-sample(5,1,prob=pm)
  inA<-k%in%c(2,4,5);inB<-k%in%c(3,4,5);inI<-k==5
  g<-r^2/rchisq(3,1);y<-rnorm(N)+inA*xA*rnorm(1,0,sqrt(g[1]))+inB*xB*rnorm(1,0,sqrt(g[2]))+inI*xI*rnorm(1,0,sqrt(g[3]))
  s<-ss(y,xA,xB,xI); e<-evid(s[1],s[2],s[3])
  bf<-(e[2]+e[4])/(e[1]+e[3]); pc<-sum((pm*e)[c(2,4,5)])/sum(pm*e)
  res[i,]<-c(bf/(1+bf),pc,inA,inI,e[5]/(e[4]+e[5]),k)};res}
cal<-function(p,t,lab){b<-cut(p,c(0,.1,.3,.5,.7,.9,.99,1),include.lowest=TRUE);cat(" ",lab,"\n");print(round(rbind(annonce=tapply(p,b,mean),realise=tapply(t,b,mean),eff=tapply(t,b,length)),3))
 cat("   moyenne p annonce=",round(mean(p),3)," frequence vraie=",round(mean(t),3),"\n")}
set.seed(13); t0<-Sys.time()
R<-simC(c(.2,.2,.2,.2,.2),1500); cat(" temps",format(Sys.time()-t0),"\n")
cal(R[,1],R[,3],"Generation uniforme 5 modeles, P(A) apparie, A dans le modele vrai")
cal(R[,2],R[,3],"Temoin : P(A) coherent uniforme")
sub<-R[,6]!=5; cal(R[sub,1],R[sub,3],"Idem P(A) apparie, restreint aux modeles sans interaction (la seule lecture correcte)")
R2<-simC(c(.25,.25,.25,.125,.125),1500); cal(R2[,1],R2[,3],"Generation Chipman (ce que veut dire 'pi=1/2 par terme'), P(A) apparie")
cat(" Dans ces simulations : P(A)<.5 et A vrai parmi donnees a A:B vrai :",round(mean(R[R[,6]==5,1]<.5),3),"\n")

cat("\n=== D. Multiplicite : 2^3 sous H0 globale, 7 termes, pi=1/2 apparie (MC a priori)\n")
mult<-function(rp,Rn=300,Mg=3000,th=c(.5,.75,.9,.95)){
 des<-expand.grid(a=c(-1,1),b=c(-1,1),c=c(-1,1),r=1:rp); Nn<-nrow(des)
 cols<-with(des,list(A=a,B=b,C=c,AB=a*b,AC=a*c,BC=b*c,ABC=a*b*c)); cn<-names(cols)
 ord<-c(1,1,1,2,2,2,3); nt<-Nn/2^ord
 Xs<-lapply(cols,function(z) z/sqrt(sum(z^2)))
 allm<-as.matrix(expand.grid(rep(list(0:1),7))); colnames(allm)<-cn
 ok<-apply(allm,1,function(z){ (z["AB"]<=z["A"]*z["B"])&(z["AC"]<=z["A"]*z["C"])&(z["BC"]<=z["B"]*z["C"])&(z["ABC"]<=z["AB"]*z["AC"]*z["BC"])})
 Ms<-allm[ok,]; stopifnot(nrow(Ms)==19)
 hi<-list(A=c("AB","AC","ABC"),B=c("AB","BC","ABC"),C=c("AC","BC","ABC"),AB="ABC",AC="ABC",BC="ABC",ABC=character())
 key<-function(z) paste(z,collapse="")
 keys<-apply(Ms,1,key)
 mp<-lapply(cn,function(t){plus<-which(Ms[,t]==1&rowSums(Ms[,hi[[t]],drop=FALSE])==0); minus<-sapply(plus,function(i){z<-Ms[i,];z[t]<-0;match(key(z),keys)}); list(p=plus,m=minus)}); names(mp)<-cn
 G<-matrix(.25/rchisq(7*Mg,1),Mg,7)
 W<-sweep(G,2,nt,function(g,n) n*g/(1+n*g)); P<-sweep(G,2,nt,function(g,n) -.5*log(1+n*g))
 pv<-matrix(NA,Rn,7); res<-matrix(NA,Rn,7)
 for(i in 1:Rn){y<-rnorm(Nn);yc<-y-mean(y);S<-sum(yc^2);ssv<-sapply(Xs,function(x) sum(yc*x)^2)/S
  ev<-apply(Ms,1,function(z){ if(sum(z)==0) return(1); id<-which(z==1); mean(exp(rowSums(P[,id,drop=FALSE])-(Nn-1)/2*log(pmax(1-(W[,id,drop=FALSE]%*%ssv[id]),1e-12))))})
  res[i,]<-sapply(cn,function(t){b<-sum(ev[mp[[t]]$p])/sum(ev[mp[[t]]$m]);b/(1+b)})
  sse<-S*(1-sum(ssv)); pv[i,]<-1-pf(ssv*S/(sse/(Nn-8)),1,Nn-8)}
 cat(sprintf(" N=%d (%d rep/cellule)\n",Nn,rp))
 for(t in th) cat(sprintf("   seuil P>%.2f : par terme %.3f ; au moins 1 sur 7 : %.3f\n",t,mean(res>t),mean(apply(res>t,1,any))))
 cat(sprintf("   p<.05 : par terme %.3f ; au moins 1 sur 7 : %.3f\n",mean(pv<.05),mean(apply(pv<.05,1,any))))
 cat(sprintf("   nombre moyen de termes 'Oui' (P>.5) : %.2f ; mediane du max des P sur 7 : %.3f ; P mediane par terme : %.3f\n",mean(rowSums(res>.5)),median(apply(res,1,max)),median(res)))}
set.seed(14); mult(2); mult(5)
cat("\nFIN\n")
