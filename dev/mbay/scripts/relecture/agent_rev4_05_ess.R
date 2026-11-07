source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev4_lib.R")
show <- function(lab,y,g,seeds=1:3,draws=4000L,u=1.5){ st<-.bay_laws_stats(y,g)
 for(sd in seeds){ t0<-Sys.time(); tab<-.bay_laws_table(st,c(.25,.5,1),u,draws,sd)
  a<-.bay_laws_average(tab,.5); z<-.bay_laws_average(tab,.5,h1_form="zs")
  sets<-sapply(c(.25,.5,1),function(r) paste(.bay_laws_average(tab,r)$laws_dropped,collapse="+"))
  cat(sprintf("%s graine=%d draws=%d (%.0fs)\n  P=%.4f err=%.1f%% lois=%s\n  ecartees r=.25/.5/1 : [%s] [%s] [%s] ; forme 1 seule ecarte [%s] P=%.4f\n  ess: %s\n  logBF/loi %s\n",lab,sd,draws,as.numeric(Sys.time()-t0,units="secs"),a$p_h1,a$error_pct,paste(sprintf("%s:%.3f",names(a$p_law),a$p_law),collapse=" "),sets[1],sets[2],sets[3],paste(z$laws_dropped,collapse="+"),z$p_h1,paste(sprintf("%s/%s:%.0f",tab$law[tab$r==.5],tab$hypothesis[tab$r==.5],tab$ess[tab$r==.5]),collapse=" "),paste(sprintf("%.2f",a$log_bf_law),collapse="/")));flush.console()}
 invisible(NULL)}
mk <- function(k,n) rep(sprintf("g%02d",1:k),each=n)
# 2. k=8, variances tres differentes
set.seed(32); k<-8;n<-6;g<-mk(k,n);sds<-10^seq(-2,2,length.out=k);y<-unlist(lapply(sds,function(s) rnorm(n,0,s)))
show("k=8 H0 sd 0,01..100",y,g,seeds=1:3)
show("k=8 H0 sd 0,01..100 draws max",y,g,seeds=1,draws=40000L)
# 3. groupe tres eloigne : normal_het ajuste sur le melange, la colonne forme 1 suit-elle ?
set.seed(33); k<-8;n<-8;g<-mk(k,n);y<-rnorm(k*n);y[g=="g08"]<-y[g=="g08"]+1000
show("k=8 d=1000",y,g,seeds=1:2)
set.seed(34); k<-5;n<-4;g<-mk(k,n);y<-rnorm(k*n);y[g=="g05"]<-y[g=="g05"]+50
show("k=5 n=4 d=50",y,g,seeds=1:2)
# 4. tres petits groupes, singleton
set.seed(35); g<-c("a",rep("b",3),rep("c",3),rep("d",3),"e"); y<-rnorm(length(g)); show("singletons",y,g,seeds=1:3)
# 5. petit u
set.seed(36); g<-mk(3,10); y<-rnorm(30); show("u=0,01",y,g,seeds=1:2,u=.01)
cat("--- 6. une seule piece t faible : ESS par forme, k=8, queues tres lourdes\n")
set.seed(37); g<-mk(8,5);y<-rt(40,1.5); st<-.bay_laws_stats(y,g)
for(het in c(FALSE,TRUE)){ set.seed(1);x<-.bay_t_ml(st,TRUE,het,1.5,.5,2000L);cat(sprintf("het=%s tous: %.2f ess %.0f\n",het,x$log_ml,x$ess));for(s in 1:8){x<-.bay_t_ml(st,TRUE,het,1.5,.5,2000L,s=s);cat(sprintf("  s=%d %.2f ess %.0f\n",s,x$log_ml,x$ess))}}
cat("--- 7. rscale petit, 4 lois : ecart annonce ou non ?\n")
set.seed(38); g<-mk(3,6);y<-rnorm(18)
for(r in c(.001,.01,.05,.1,.25,.5)) try_it(sprintf("rscale=%g",r),{v<-m.bay(y,g,rscale=r,verbose=FALSE,sensitivity=FALSE,draws=4000L);cat(paste(substr(v$notes$fr,1,60),collapse=" | "),"\n");v})
# 1. k=20, variances tres differentes (I4 non refait), le plus long, en dernier
set.seed(31); k<-20;n<-6;g<-mk(k,n);sds<-10^seq(-1.5,1.5,length.out=k);y<-unlist(lapply(sds,function(s) rnorm(n,0,s)))
show("k=20 H0 sd 0,03..30",y,g,seeds=1)
