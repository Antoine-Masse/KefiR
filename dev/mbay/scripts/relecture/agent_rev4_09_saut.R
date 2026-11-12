source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev4_lib.R")
n<-20; g<-rep(c("a","b","c"),each=n)
set.seed(3); b<-round(c(rnorm(n,0,1),rnorm(n,.3,1),rnorm(n,-.2,1)),3)
cat("--- A. saut au seuil : 10 puis 11 ex-aequo (n=20) a la meme valeur ; 11 force a 4 lois (SS>0, integrale finie)\n")
for(val in c(0,1,3)) for(m in c(9,10,11,12)){ y<-b; y[g=="a"][1:m]<-val; st<-.bay_laws_stats(y,g); v<-m.bay(y,g,verbose=FALSE,sensitivity=FALSE)
 tab<-.bay_laws_table(st,.5,1.5,4000L,1L); a<-.bay_laws_average(tab,.5)
 cat(sprintf("val=%g m=%d SS_a=%.3f m.bay P=%.4f logBF=%.2f [%s] | 4 lois forcees P=%.4f logBF=%.2f err=%.1f%% lois %s ess %s\n",val,m,st$ss[1],v$prob_effect,v$effects$log_bf,paste(v$laws$law,collapse=","),a$p_h1,a$log_bf,a$error_pct,paste(sprintf("%.3f",a$p_law),collapse="/"),paste(round(tab$ess),collapse="/")));flush.console()}
cat("--- B. seuil global : 50 puis 51 % de l'echantillon (aucun groupe > 50 %)\n")
y<-b; for(gr in c("a","b","c")) y[g==gr][1:10]<-0; st<-.bay_laws_stats(y,g); print(.bay_laws_allowed(st,.bay_law_names())); try_it("30/60",m.bay(y,g,verbose=FALSE,sensitivity=FALSE))
cat("    par groupe tout a 50 % a valeurs distinctes, 30/60 mais aucun commun\n")
y<-b; y[g=="a"][1:10]<-0; y[g=="b"][1:10]<-1; y[g=="c"][1:10]<-2; try_it("3x50% distincts",m.bay(y,g,verbose=FALSE,sensitivity=FALSE))
cat("--- C. grille brute t_eq en sinh autour des valeurs repetees contre IS\n")
sinh_grid<-function(centers,lo,hi,A=16,dx=.2,base=1e-6){ x<-sinh(seq(-A,A,by=dx))*base; pts<-sort(unique(c(lo,hi,seq(lo,hi,length.out=60),as.vector(outer(centers,x,"+"))))); pts<-pts[pts>=lo&pts<=hi]; w<-c(diff(pts)/2,0)+c(0,diff(pts)/2); list(x=pts,w=w)}
brute<-function(st,centers,s=NULL,r=.5){ lo<-min(unlist(st$y))-1; hi<-max(unlist(st$y))+1; Gm<-sinh_grid(centers,lo,hi);
 Gt<-list(x=seq(-4,40,by=.25));Gt$w<-rep(.25,length(Gt$x)); Gv<-list(x=seq(-14,5,by=.5));Gv$w<-rep(.5,length(Gv$x))
 if(is.null(s)){ X<-as.matrix(expand.grid(Gm$x,Gt$x,Gv$x)); lw<-log(as.vector(outer(outer(Gm$w,Gt$w),Gv$w))); keep<-lw>-Inf; h1<-FALSE } else {
   ka<-sinh_grid(c(0,outer(centers,centers,"-")),lo-hi,hi-lo,A=14,dx=.5); Gm2<-sinh_grid(centers,lo,hi,A=14,dx=.5)
   Gt$x<-seq(-4,30,by=.5);Gt$w<-rep(.5,length(Gt$x));Gv$x<-seq(-10,4,by=1);Gv$w<-rep(1,length(Gv$x))
   X<-as.matrix(expand.grid(Gm2$x,ka$x,Gt$x,Gv$x)); lw<-log(as.vector(outer(outer(outer(Gm2$w,ka$w),Gt$w),Gv$w))); keep<-lw>-Inf; h1<-TRUE
   cat("  taille grille",nrow(X),"\n") }
 X<-X[keep,,drop=FALSE]; lw<-lw[keep]; out<-numeric(nrow(X))
 for(b in split(seq_len(nrow(X)),ceiling(seq_len(nrow(X))/2e5))){ pa<-.bay_t_parts(X[b,,drop=FALSE],st,h1,FALSE,1.5,s); out[b]<-pa$base+(if(h1).bay_t_prior_theta(pa,r,st$k) else 0)+lw[b]}
 ok<-is.finite(out); m<-max(out[ok]); w<-exp(out[ok]-m); top<-X[ok,,drop=FALSE][which.max(out[ok]),]
 cat("  mode de la grille :",round(top,3)," part de la plus forte cellule",round(max(w)/sum(w),4),"\n"); m+log(sum(w))}
run<-function(lab,y,centers_raw){ st<-.bay_laws_stats(y,g); un<-sapply(centers_raw,function(v) unique(st$y[[1]][y[g=="a"]==v])[1]); cat("##",lab,"\n")
 r0<-brute(st,un); is0<-sapply(1:3,function(sd){set.seed(sd);x<-.bay_t_ml(st,FALSE,FALSE,1.5,.5,4000L);c(x$log_ml,x$ess)})
 cat(sprintf(" H0 grille %.3f ; IS %s ; ess %s\n",r0,paste(round(is0[1,],3),collapse="/"),paste(round(is0[2,]),collapse="/")));flush.console()
 r1<-brute(st,un,s=1L); is1<-sapply(1:3,function(sd){set.seed(sd);x<-.bay_t_ml(st,TRUE,FALSE,1.5,.5,4000L,s=1L);c(x$log_ml,x$ess)})
 cat(sprintf(" s=1 grille %.3f ; IS %s ; ess %s\n",r1,paste(round(is1[1,],3),collapse="/"),paste(round(is1[2,]),collapse="/")));flush.console()}
y<-b; for(gr in c("a","b","c")) y[g==gr][1:10]<-0; run("3 groupes 50 % a 0",y,0)
y<-b; for(gr in c("a","b","c")) y[g==gr][1:9]<-0; run("3 groupes 45 % a 0",y,0)
y<-b; for(gr in c("a","b","c")){id<-which(g==gr); y[id[1:9]]<- -1; y[id[10:18]]<- 1}; run("bimodal 45+45",y,c(-1,1))
