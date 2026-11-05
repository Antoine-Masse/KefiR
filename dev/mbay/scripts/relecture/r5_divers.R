.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); options(BFprogress=FALSE)
B<-function(e) extractBF(e)$bf
cat("=== 1. Interaction moyenne : P(A:B) sous 1/2 => A lu sans interaction ; retrait implicite\n")
set.seed(61); it<-1e5
d<-expand.grid(A=gl(2,1),B=gl(2,1),r=1:12); sa<-ifelse(d$A==1,-1,1); sb<-ifelse(d$B==1,-1,1)
best<-NULL
for(s in 1:200){set.seed(600+s); y<-.0*sa+.32*sa*sb+rnorm(48); pI<-summary(lm(y~sa*sb))$coef[4,4]; pA<-summary(lm(y~sa*sb))$coef[2,4]; if(pI>.03&pI<.08&pA>.3){best<-y;break}}
d$y<-best
m<-list(A=B(lmBF(y~A,d,iterations=it)),B=B(lmBF(y~B,d,iterations=it)),AB=B(lmBF(y~A+B,d,iterations=it)),F=B(lmBF(y~A*B,d,iterations=it)))
bfA<-(m$A+m$AB)/(1+m$B);bfI<-m$F/m$AB;tot<-1+m$A+m$B+m$AB+m$F
print(anova(lm(y~A*B,d)))
cat(sprintf(" P(A:B | A,B)=%.3f (donc pas d'alerte si <.5) ; P(A apparie)=%.3f ; P(A uniforme 5)=%.3f ; P(A Chipman)=%.3f\n",bfI/(1+bfI),bfA/(1+bfA),(m$A+m$AB+m$F)/tot,(m$A/4+m$AB/8+m$F/8)/(1/4+m$A/4+m$B/4+m$AB/8+m$F/8)))

cat("\n=== 2. Calibration 2 groupes, n=30/gr, p annonce JZS r=.707 pi=.5, sous d'autres mondes (20000 jeux)\n")
jzs<-function(t,N,nu,r=sqrt(2)/2){u<-seq(-12,14,length.out=260);g<-exp(u);dg<-(r^2/2)^.5/gamma(.5)*g^(-1.5)*exp(-r^2/(2*g))*g*(u[2]-u[1])
 sapply(t,function(tt) sum(dg*(1+N*g)^(-.5)*((1+tt^2/((1+N*g)*nu))/(1+tt^2/nu))^(-(nu+1)/2)))}
n<-30;N<-n/2;nu<-2*n-2
monde<-function(pi1,rdelta,R=20000){h<-rbinom(R,1,pi1);del<-h*rdelta();tt<-rt(R,nu,ncp=del*sqrt(N));bf<-jzs(tt,N,nu);P<-bf/(1+bf);b<-cut(P,c(0,.25,.5,.75,.9,.95,1),include.lowest=TRUE);round(rbind(annonce=tapply(P,b,mean),realise=tapply(h,b,mean),eff=tapply(h,b,length)),3)}
set.seed(62);R<-20000
cat(" Monde defaut (pi=.5, Cauchy .707) :\n");print(monde(.5,function() rcauchy(R,0,.707)))
cat(" pi=.5, delta~N(0,.2^2) (effets petits) :\n");print(monde(.5,function() rnorm(R,0,.2)))
cat(" pi=.1, delta~N(0,.5^2) (criblage, peu de vrais effets) :\n");print(monde(.1,function() rnorm(R,0,.5)))
cat(" pi=.9, delta~N(0,3^2) (procede, gros effets) :\n");print(monde(.9,function() rnorm(R,0,3)))

cat("\n=== 3. Apparie : ttestBF sur differences vs lmBF A+id vs id, meme echelle r\n")
set.seed(63);res<-t(replicate(6,{ni<-25;u<-rnorm(ni);y1<-u+rnorm(ni,0,.5);y2<-u+.25+rnorm(ni,0,.5);dp<-data.frame(y=c(y1,y2),t=gl(2,ni),id=factor(rep(1:ni,2)));c(tt=B(ttestBF(y2,y1,paired=TRUE)),lm=B(lmBF(y~t+id,dp,whichRandom="id",iterations=2e4)/lmBF(y~id,dp,whichRandom="id",iterations=2e4)),p=t.test(y2,y1,paired=TRUE)$p.value)}))
print(round(res,3))

cat("\n=== 4. Moyenne des P sur grille r x cotes pour un jeu limite : la grille conclut-elle ?\n")
set.seed(64);a<-gl(3,15);y<-rnorm(45)+.45*(as.integer(a)-2);p<-summary(aov(y~a))[[1]][1,5]
g<-sapply(c(.25,.5,1),function(r) B(anovaBF(y~a,data.frame(y,a),rscaleFixed=r)))
cat(" p aov",signif(p,3)," BF(r/2,r,2r)",round(g,3),"\n");print(round(outer(g,c(.25,1,4),function(b,o) o*b/(1+o*b)),3))
cat(" Seuil BF pour 'stable' sur cotes 1:4..4:1 : BF>4 ou BF<1/4 ; plancher BF a n=45 (F=0)=0.164 ; a n=15 (F=0, 5/gr)=.310 => jamais 'stable pour non' a 5/gr\n")

cat("\n=== 5. Covariable forte : meme effet brut, sigma change, ROPE et P\n")
set.seed(65);n<-40;x<-rnorm(2*n);tr<-gl(2,n);eff<-.15;y<-eff*(tr==2)+5*x+rnorm(2*n,0,.5);dc<-data.frame(y,tr,x)
s<-sd(resid(lm(y~tr+x)));cat(sprintf(" Sans x : sd=%.2f delta=%.3f P(T)=%.3f ; avec x : sd=%.2f delta=%.3f P(T|x)=%.3f\n",sd(resid(lm(y~tr))),eff/sd(resid(lm(y~tr))),{b<-B(lmBF(y~tr,dc));b/(1+b)},s,eff/s,{b<-B(lmBF(y~tr+x,dc)/lmBF(y~x,dc));b/(1+b)}))
cat("\nFIN\n")
