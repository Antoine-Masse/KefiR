source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev4_lib.R")
cut_mix <- .bay_t_h1_mix
full_mix <- cut_mix; b <- deparse(body(cut_mix)); b <- sub("if \\(k > 8L\\)","if (FALSE)",b); body(full_mix) <- parse(text=b)[[1]]; environment(full_mix)<-globalenv()
run2 <- function(lab,y,g,laws=.bay_law_names(),seed=1L,draws=2000L){
  t0<-Sys.time()
  assign(".bay_t_h1_mix",cut_mix,envir=globalenv()); a<-m.bay(y,g,verbose=FALSE,sensitivity=FALSE,laws=laws,seed=seed,draws=draws)
  t1<-Sys.time()
  assign(".bay_t_h1_mix",full_mix,envir=globalenv()); b<-m.bay(y,g,verbose=FALSE,sensitivity=FALSE,laws=laws,seed=seed,draws=draws)
  t2<-Sys.time()
  assign(".bay_t_h1_mix",cut_mix,envir=globalenv())
  st<-.bay_laws_stats(y,g); ord<-order(st$m); mb<-st$m[ord][which(cumsum(st$n[ord])/st$N>=.5)[1]]; sel<-order(-abs(st$m-mb)*sqrt(st$n))[1:min(5,st$k)]
  cat(sprintf("%s k=%d groupes coupe=%s\n  coupe : P=%.4f logBF=%.2f lois %s (%.0fs)\n  entier: P=%.4f logBF=%.2f lois %s (%.0fs)\n  ecart logBF (coupe-entier)=%.2f ; log(k/5)=%.2f\n  logBF/loi coupe %s ; entier %s\n",lab,st$k,paste(sel,collapse=","),a$prob_effect,a$effects$log_bf,paste(round(a$laws$posterior,3),collapse="/"),as.numeric(t1-t0,units="secs"),b$prob_effect,b$effects$log_bf,paste(round(b$laws$posterior,3),collapse="/"),as.numeric(t2-t1,units="secs"),a$effects$log_bf-b$effects$log_bf,log(st$k/5),paste(round(a$laws$log_bf,2),collapse="/"),paste(round(b$laws$log_bf,2),collapse="/")))
  flush.console(); invisible(list(a=a,b=b))
}
mk <- function(k,n) rep(sprintf("g%02d",1:k),each=n)
set.seed(21); run2("H0 normale",rnorm(72),mk(9,8))
set.seed(22); run2("H0 normale",rnorm(96),mk(12,8))
set.seed(23); run2("H0 t3",rt(96,3),mk(12,8))
set.seed(24); y<-rt(96,3); y[mk(12,8)=="g12"]<-y[mk(12,8)=="g12"]+1.5; run2("un groupe decale 1,5 ; t3",y,mk(12,8))
# Aberrants : six groupes ont une valeur extreme qui gonfle leur moyenne ; le vrai decale (g12) a une moyenne moindre
set.seed(25); g<-mk(12,10); y<-rnorm(120,0,.5); for(i in 1:6) y[g==sprintf("g%02d",i)][1]<-20; y[g=="g12"]<-y[g=="g12"]+1.2
print(round(tapply(y,g,mean),2)); r<-run2("aberrants+g12 decale, lois t",y,g,laws=c("t_eq","t_het")); r<-run2("aberrants+g12 decale, 4 lois",y,g)
# Variance tres petite : ecart modeste en unite commune, enorme en unite propre
set.seed(26); g<-mk(12,8); y<-rnorm(96,0,3); y[g=="g12"]<-.8+rnorm(8,0,.02); print(round(tapply(y,g,sd),2)); r<-run2("g12 tres serre, ecart modeste",y,g)
set.seed(27); g<-mk(16,6); run2("H0 normale",rnorm(96),g,laws=c("normal_eq","t_eq"))
cat("--- continuite k=8/k=9 : un groupe de plus sous H0\n")
set.seed(28); g<-mk(9,8); y<-rnorm(72); run2("k=9",y,g,laws=c("t_eq"));a<-m.bay(y[g!="g09"],g[g!="g09"],verbose=FALSE,sensitivity=FALSE,laws="t_eq",draws=2000L);cat("k=8 sans g09 : logBF",a$effects$log_bf,"\n")
