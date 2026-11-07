source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev4_lib.R")
cat("--- A. une paire d'ex-aequo dans un petit groupe\n")
y<-c(10,10,11, 1,8,20, -5,3,30); g<-rep(c("a","b","c"),each=3); st<-.bay_laws_stats(y,g); print(.bay_laws_allowed(st,.bay_law_names()))
v<-m.bay(y,g,verbose=FALSE); cat("P m.bay",v$prob_effect,"\n",v$notes$fr,sep="\n")
tab<-.bay_laws_table(st,.5,1.5,4000L,1L); a<-.bay_laws_average(tab,.5); cat("P 4 lois forcees",a$p_h1,"ecartees",a$laws_dropped,"lois",round(a$p_law,3),"ess",tab$ess,"\n")
cat("--- B. vrai H0, variances inegales, entiers : part reduite a normal_eq et erreur\n")
set.seed(41); R<-150; res<-matrix(NA,R,3)
for(i in 1:R){ g<-rep(c("a","b","c"),each=6); yy<-round(c(rnorm(6,10,.5),rnorm(6,10,3),rnorm(6,10,6))); if(length(unique(yy))<3) next
 st<-.bay_laws_stats(yy,g); al<-.bay_laws_allowed(st,.bay_law_names()); res[i,1]<-length(al)
 const<-.bay_law_const(st);one<-matrix(1,1,3); h0<-const+.bay_law_kernel(st,one,0); h1<-.bay_lse2(h0+.bay_bf_oneway_zs(yy,g,.5),const+.bay_log_sum_exp(sapply(1:3,function(s) drop(.bay_law_gint(st,one,.5*sqrt(2),s=s))))-log(3))+log(.5)
 res[i,2]<-h1-h0
 if(!any(st$ss==0)){ if(i%%5==0){tb<-.bay_laws_table(st,.5,1.5,1000L,i,laws=c("normal_eq","normal_het")); res[i,3]<-.bay_laws_average(tb,.5)$log_bf}}
}
cat("lois retenues : ");print(table(res[,1])); sel<-res[,1]==1
cat(sprintf("reduit a normal_eq : %d/%d ; parmi eux P(H1)>.5 : %.3f, >.9 : %.3f ; hors reduction >.5 : %.3f\n",sum(sel,na.rm=T),sum(!is.na(res[,1])),mean(res[sel,2]>0,na.rm=T),mean(res[sel,2]>log(9),na.rm=T),mean(res[!sel,2]>0,na.rm=T)))
cat("  sur les jeux reduits avec SS>0, logBF normal_eq contre eq+het :\n"); w<-which(sel&!is.na(res[,3])); print(round(cbind(eq=res[w,2],eq_het=res[w,3]),2)[seq_len(min(12,length(w))),])
cat("--- C. ex-aequo non reconnus apres centrage (bruit d'arrondi pres de 0)\n")
e<-.Machine$double.eps*.25
v0<-.3+c(0,1,2,-1,-2,3)*e*1; print(sprintf("%.17g",v0))
set.seed(1); o<-rnorm(6); o<-o-mean(o); yy<-c(v0,.3+o,.3-o); g<-rep(c("a","b","c"),each=6); st<-.bay_laws_stats(yy,g)
cat("groupe a brut : share",max(table(yy[1:6]))/6," centre : share",max(table(st$y[[1]]))/6," ss",st$ss[1]," retenues",.bay_laws_allowed(st,.bay_law_names()),"\n")
for(sd in 1:2) try_it(paste("groupe a arrondi, graine",sd),m.bay(yy,g,verbose=FALSE,sensitivity=FALSE,seed=sd))
try_it("meme groupe a exactement constant",m.bay(c(rep(.3,6),.3+o,.3-o),g,verbose=FALSE,sensitivity=FALSE))
cat("--- D. laws reduites : het calcule quand meme ? temps et plantage\n")
g<-rep(c("a","b","c"),each=200); yy<-c(rep(1,200),rep(2,200),rnorm(200)); print(system.time(try_it("2 groupes constants n=200",m.bay(yy,g,verbose=FALSE,sensitivity=FALSE))))
g<-rep(sprintf("g%02d",1:12),each=30); set.seed(2); yy<-rnorm(360); yy[g=="g01"]<-1; print(system.time(try_it("k=12 un constant",m.bay(yy,g,verbose=FALSE,sensitivity=FALSE))))
print(system.time(try_it("k=12 sans ex-aequo laws=normal_eq",m.bay(rnorm(360),g,verbose=FALSE,sensitivity=FALSE,laws="normal_eq"))))
cat("--- E. discontinuite exact contre 1e-8\n")
set.seed(3); g<-rep(c("a","b","c"),each=20); b<-c(rnorm(20,0,1),rnorm(20,1,1),rnorm(20,-.5,1));
y1<-b; y1[g=="a"]<-3; y2<-y1; y2[g=="a"]<-3+1e-8*rnorm(20); y3<-y1; y3[g=="a"]<-3+1e-2*rnorm(20)
for(nm in c("y1","y2","y3")) try_it(nm,m.bay(get(nm),g,verbose=FALSE,sensitivity=FALSE))
