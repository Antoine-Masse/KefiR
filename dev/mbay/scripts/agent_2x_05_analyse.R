# Analyse de toutes les sorties. Ecrit resultats/agent_2x_tables.md, et une figure PNG.
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/agent_2x_commun.R")
OUT <- file.path(RES,"agent_2x_tables.md"); cat("", file=OUT)
md <- function(df, titre=NULL, dg=3) { if (!is.null(titre)) cat("\n**",titre,"**\n\n",sep="",file=OUT,append=TRUE)
  df2 <- as.data.frame(lapply(df,function(x) if(is.numeric(x)) ifelse(is.na(x),"NA",formatC(x,digits=dg,format="f")) else as.character(x)),check.names=FALSE)
  names(df2) <- names(df)
  cat(paste0("| ",paste(names(df2),collapse=" | ")," |\n|",paste(rep("---",ncol(df2)),collapse="|"),"|\n"),file=OUT,append=TRUE)
  for(i in seq_len(nrow(df2))) cat(paste0("| ",paste(unlist(df2[i,]),collapse=" | ")," |\n"),file=OUT,append=TRUE); cat("\n",file=OUT,append=TRUE) }
sec <- function(s) cat("\n### ",s,"\n",sep="",file=OUT,append=TRUE)
ci <- function(x,n) { w<-wilson(x,n); sprintf("%.3f [%.3f ; %.3f]",x/n,w[1],w[2]) }

# ---- E1
if (file.exists(file.path(RES,"agent_2x_E1.rds"))) { E1 <- readRDS(file.path(RES,"agent_2x_E1.rds")); sec("Etape 1 : sous l'a priori")
  gl1 <- list(); dec <- list()
  png(file.path(RES,"agent_2x_calibration_E1.png"),1400,720,res=110); par(mfrow=c(2,4),mar=c(4,4,2.5,1))
  for (s in E1) { P <- prob(s$l,s$pi); N <- length(P); m<-mean(P); se<-sd(P)/sqrt(N); z<-(m-s$pi)/se
    v<-P*(1-P); zc<-sum(s$H-P)/sqrt(sum(v))
    b<-cut(P,seq(0,1,.1),include.lowest=TRUE); rows<-list(); zs<-c()
    for(lv in levels(b)){ w<-which(b==lv); if(!length(w)) next; zb<-sum(s$H[w]-P[w])/sqrt(sum(v[w])); zs<-c(zs,zb)
      wl<-wilson(sum(s$H[w]),length(w)); rows[[length(rows)+1]]<-data.frame(bin=lv,n=length(w),moyP=mean(P[w]),nH1=sum(s$H[w]),freqH1=mean(s$H[w]),ICbas=wl[1],IChaut=wl[2],Z=zb,check.names=FALSE)}
    t<-do.call(rbind,rows); chi<-sum(zs^2); pc<-pchisq(chi,length(zs),lower.tail=FALSE)
    gl1[[length(gl1)+1]]<-data.frame(pi=s$pi,k=s$k,n=s$n,N=N,`moy P`=m,IC95=sprintf("[%.4f ; %.4f]",m-1.96*se,m+1.96*se),`z(moy=pi)`=z,`p`=2*pnorm(-abs(z)),`freq H1`=mean(s$H),`z(H-P)`=zc,`chi2 deciles`=chi,ddl=length(zs),`p chi2`=pc,`max|Z| bin`=max(abs(zs)),replis=sum(s$rp),`moy P sans replis`=mean(P[!s$rp]),sec=s$sec,check.names=FALSE)
    dec[[length(dec)+1]]<-cbind(pi=s$pi,k=s$k,n.grp=s$n,t)
    plot(c(0,1),c(0,1),type="n",xlab="P annoncee",ylab="frequence de H1",main=sprintf("pi=%.1f k=%d n=%d",s$pi,s$k,s$n)); abline(0,1,col="grey60")
    segments(t$moyP,t$ICbas,t$moyP,t$IChaut,col="#1f6f8b"); points(t$moyP,t$freqH1,pch=19,col="#1f6f8b",cex=.8)
  }; dev.off()
  md(do.call(rbind,gl1),"Test DAP et calibration globale (par scenario)",4); md(do.call(rbind,dec),"Deciles",3)
}
# ---- E2
op <- function(l,p) { P<-prob(l,.5); N<-length(P); data.frame(N=N,`moy P`=mean(P),`P>0,95`=ci(sum(P>.95),N),`P<0,05`=ci(sum(P<.05),N),`0,25<P<0,75`=ci(sum(P>.25&P<.75),N),`p<0,05 aov`=ci(sum(p<.05),N),check.names=FALSE) }
if (file.exists(file.path(RES,"agent_2x_E2.rds"))) { E2<-readRDS(file.path(RES,"agent_2x_E2.rds")); sec("Etape 2 : effets fixes")
  mk<-function(q) do.call(rbind,lapply(q,function(s) cbind(data.frame(k=s$k,n=s$n,delta=s$d,f=if(s$k==2)s$d/2 else s$d/sqrt(6),erreur=s$err,replis=s$rep),op(s$l,s$p))))
  md(mk(Filter(function(s)s$err=="normal",E2)),"Erreur normale"); nn<-Filter(function(s)s$err!="normal",E2); if(length(nn)) md(mk(nn),"Non-normalite, n=30 (reference normale : tableau precedent)") }
if (file.exists(file.path(RES,"agent_2x_EH.rds"))) { H<-readRDS(file.path(RES,"agent_2x_EH.rds"))
  md(do.call(rbind,lapply(H,function(s) cbind(data.frame(cas=s$nom,delta=s$d,`NA anovaBF`=sum(is.na(s$l)),sec=s$sec,check.names=FALSE),op(s$l[!is.na(s$l)],s$p[!is.na(s$l)])))),"Heteroscedasticite k=3, n=10,30,50 (anovaBF complet)") }
# ---- E3
if (file.exists(file.path(RES,"agent_2x_E3.rds"))) { E3<-readRDS(file.path(RES,"agent_2x_E3.rds")); sec("Etape 3 : sensibilite a r (memes jeux, pi=1/2, r sur l'echelle reduite ; k=2 : r_delta = r*sqrt2)")
  md(do.call(rbind,lapply(E3,function(s){P<-plogis(s$L);N<-nrow(P);C<-P>.5
    data.frame(k=s$k,n=s$n,delta=s$d,N=N,`moy P r=.25`=mean(P[,1]),`moy P r=.5`=mean(P[,2]),`moy P r=1`=mean(P[,3]),
      `P>1/2 r=.25`=mean(C[,1]),`P>1/2 r=.5`=mean(C[,2]),`P>1/2 r=1`=mean(C[,3]),
      `change .25 vs .5`=ci(sum(C[,1]!=C[,2]),N),`change 1 vs .5`=ci(sum(C[,3]!=C[,2]),N),`change qqc`=ci(sum(C[,1]!=C[,2]|C[,3]!=C[,2]),N),
      `med |P(.25)-P(1)|`=median(abs(P[,1]-P[,3])),`max |P(.25)-P(1)|`=max(abs(P[,1]-P[,3])),check.names=FALSE)})),NULL) }
# ---- E4
fA<-file.path(RES,"agent_2x_04_AxB.rds")
if (file.exists(fA)) { D<-as.data.frame(readRDS(fA)); sec("Etape 4 : A*B 2x3, n=10 par cellule")
  N_m<-as.numeric(table(factor(D$model,1:5))); ls2<-function(a,b) pmax(a,b)+log1p(exp(pmin(a,b)-pmax(a,b)))
  Lm<-as.matrix(D[,c("l0","lA","lB","lAB","lAxB")])
  Pm<-data.frame(A=plogis(ls2(D$lA,D$lAB)-ls2(D$l0,D$lB)),B=plogis(ls2(D$lB,D$lAB)-ls2(D$l0,D$lA)),AxB=plogis(D$lAxB-D$lAB))
  inc<-rbind(A=c(0,1,0,1,1),B=c(0,0,1,1,1),AxB=c(0,0,0,0,1))
  cat(sprintf("\nJeux : %d (par modele vrai : %s). Duree moyenne par jeu : %.3f s. Replis (A ou B seuls, integrate echoue) : %d. NA restants : %d. Erreur MC relative mediane/max de A + B : %.4f/%.4f ; A*B : %.4f/%.4f.\n",nrow(D),paste(N_m,collapse=", "),max(D$sec)/nrow(D),sum(D$nrep),sum(is.na(Lm)),median(D$eAB,na.rm=T),max(D$eAB,na.rm=T),median(D$eAxB,na.rm=T),max(D$eAxB,na.rm=T)),file=OUT,append=TRUE)
  ko<-!complete.cases(Pm); if(any(ko)) { D<-D[!ko,]; Pm<-Pm[!ko,]; Lm<-Lm[!ko,]; N_m<-as.numeric(table(factor(D$model,1:5))) }
  tab<-as.data.frame(t(sapply(1:5,function(m) c(modele=m,colMeans(Pm[D$model==m,])))));tab$modele<-c("nul","A","B","A+B","A*B"); md(tab,"Moyenne des P appariees par modele vrai")
  bins<-function(P) cut(P,seq(0,1,.1),include.lowest=TRUE)
  inc_T<-function(T) sapply(1:5,function(m) T[D$model==m][1])
  stat<-function(P,T,w,minn=10){ # P, verite T (par jeu), poids w par strate. O-E par decile avec variance stratifiee
    b<-bins(P); r<-list()
    for(lv in levels(b)){ I<-b==lv; U<-(T-P)*I; a<-sapply(1:5,function(m) mean(U[D$model==m])); v<-sapply(1:5,function(m) sum((P*(1-P))[I&D$model==m])/N_m[m]^2) # variance SOUS calibration, y|P~Bern(P)
      nn<-sum(I & w[D$model]>0); pb<-sum(w*sapply(1:5,function(m) mean(I[D$model==m]))); oe<-sum(w*a)
      mp<-sum(w*sapply(1:5,function(m) mean((P*I)[D$model==m])))/pb; fr<-sum(w*sapply(1:5,function(m) mean((T*I)[D$model==m])))/pb
      r[[lv]]<-data.frame(bin=lv,n=nn,masse=pb,moyP=mp,freqT=fr,Z=if(nn>=minn) oe/sqrt(sum(w^2*v)) else NA)}
    t<-do.call(rbind,r); a<-sapply(1:5,function(m) mean((T-P)[D$model==m])); v<-sapply(1:5,function(m) var((T-P)[D$model==m])/N_m[m])
    zs<-t$Z[!is.na(t$Z)&is.finite(t$Z)]
    list(t=t,mP=sum(w*sapply(1:5,function(m) mean(P[D$model==m]))),pT=sum(w*inc_T(T)),d=sum(w*a),se=sqrt(sum(w^2*v)),chi=sum(zs^2),df=length(zs)) }
  G<-list(`G1 Chipman (1/4,1/4,1/4,1/8,1/8)`=c(.25,.25,.25,.125,.125),`G2 uniforme 5 modeles (rejet)`=rep(.2,5),`G3 additif (1/4 x4, A*B exclu)`=c(.25,.25,.25,.25,0),`G4 A+B ou A*B (1/2,1/2)`=c(0,0,0,.5,.5),`G5 nul ou A*B (seul a inclusions 1/2)`=c(.5,0,0,0,.5))
  rows<-list(); det<-list()
  for(g in names(G)){w<-G[[g]]; pst<-function(Ti) { lw<-sweep(Lm,2,log(pmax(w,1e-300)),"+"); mx<-apply(lw,1,max); e<-exp(lw-mx); rowSums(sweep(e,2,Ti,"*"))/rowSums(e)}
    for(tn in c("A","B","AxB")){T<-inc[tn,D$model]; s<-stat(Pm[[tn]],T,w); pv<-pchisq(s$chi,s$df,lower.tail=FALSE)
      rows[[length(rows)+1]]<-data.frame(generateur=g,terme=tn,proba="appariee",`P(T) vraie`=s$pT,`moy P`=s$mP,`moy P - P(T)`=-s$d,se=s$se,z=-s$d/s$se,`chi2 dec`=s$chi,ddl=s$df,`p chi2`=pv,check.names=FALSE)
      det[[length(det)+1]]<-cbind(generateur=g,terme=tn,s$t)
      s2<-stat(pst(inc[tn,]),T,w); pv2<-pchisq(s2$chi,s2$df,lower.tail=FALSE)
      rows[[length(rows)+1]]<-data.frame(generateur=g,terme=tn,proba="coherente (controle)",`P(T) vraie`=s2$pT,`moy P`=s2$mP,`moy P - P(T)`=-s2$d,se=s2$se,z=-s2$d/s2$se,`chi2 dec`=s2$chi,ddl=s2$df,`p chi2`=pv2,check.names=FALSE)
    }}
  md(do.call(rbind,rows),"Calibration par generateur (strates reponderees)")
  dd<-do.call(rbind,det); md(dd[dd$generateur %in% names(G)[c(1,3,4)],],"Deciles des P appariees : G1, G3, G4 (Z=NA si moins de 10 jeux)")
  # recherche post hoc sur tout melange des 5 modeles : minimise le chi2 total des 3 termes
  pre<-lapply(c("A","B","AxB"),function(tn){T<-inc[tn,D$model];P<-Pm[[tn]];b<-bins(P); lapply(levels(b),function(lv){I<-b==lv;U<-(T-P)*I; list(a=sapply(1:5,function(m) mean(U[D$model==m])),v=sapply(1:5,function(m) sum((P*(1-P))[I&D$model==m])/N_m[m]^2),c=sapply(1:5,function(m) sum(I[D$model==m])))})})
  pl<-unlist(pre,recursive=FALSE)
  a<-as.matrix(do.call(rbind,lapply(pl,`[[`,"a"))); v<-as.matrix(do.call(rbind,lapply(pl,`[[`,"v"))); cc<-as.matrix(do.call(rbind,lapply(pl,`[[`,"c")))
  gr<-expand.grid(a=0:20,b=0:20,c=0:20,d=0:20); gr<-gr[rowSums(gr)<=20,]; gr$e<-20-rowSums(gr); W<-as.matrix(gr)/20
  chiW<-function(w){ k<-as.numeric(cc%*%(w>0))>=10 & as.numeric(v%*%w^2)>0; sum(((a%*%w)^2/(v%*%w^2))[k]) }; dfW<-function(w) sum(as.numeric(cc%*%(w>0))>=10 & as.numeric(v%*%w^2)>0)
  ok<-W[,5]>0 & W[,5]<1; Wk<-W[ok,]; ch<-apply(Wk,1,chiW); dfw<-apply(Wk,1,dfW); pv<-pchisq(ch,dfw,lower.tail=FALSE)
  o<-order(pv,decreasing=TRUE)[1:5]; res<-data.frame(Wk[o,,drop=FALSE],chi2=ch[o],ddl=dfw[o],p=pv[o]); names(res)[1:5]<-c("w_nul","w_A","w_B","w_A+B","w_A*B")
  md(res,sprintf("Recherche post hoc sur %d melanges (pas 0,05, 0<w_A*B<1) : 5 meilleurs (tests multiples, sans valeur inferentielle)",nrow(Wk)))
  cat(sprintf("\nMelanges avec p chi2 > 0,05 : %d sur %d.\n",sum(pv>.05),nrow(Wk)),file=OUT,append=TRUE)
}
cat("analyse OK\n")
