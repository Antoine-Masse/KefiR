# Agent 25b : suite de agent_25 apres a) (replanifie pour le delai) : b) N=60 puis c) jusqu'a l'echeance. Memes formules de graines.
source("agent_24_commun.R")
FIN <- as.POSIXct(paste(Sys.Date(), "14:15:15"))
OUT <- readRDS(paste0(SORTIES, "agent_25_v3abc.rds")); sauve <- function() saveRDS(OUT, paste0(SORTIES, "agent_25_v3abc.rds"))
ZS1 <- function(y,g) tryCatch(suppressWarnings(as.numeric(ZS(y,g,.5))), error=function(e) NA_real_)
cells <- expand.grid(k=c(2,3,5), n=c(5,20,100)); br <- vector("list", nrow(cells)); OUT$b_cells<-cells
for (rep in 1:60) { for (ci in seq_len(nrow(cells))) { k<-cells$k[ci]; n<-cells$n[ci]; gj<-3000000L+ci*1000L+rep; set.seed(gj)
    g<-factor(rep(1:k,each=n)); y<-rnorm(k*n); r<-het_safe(y,g,.5,500L,gj)
    br[[ci]] <- rbind(br[[ci]], c(rep=rep,het=r$log_bf,zs=ZS1(y,g),err=r$error_pct,ess=r$ess,aov=anova(lm(y~g))[1,5],welch=oneway.test(y~g)$p.value)) } 
  if (rep %% 20==0) { OUT$b<-br; sauve(); cat("b",rep,format(Sys.time()),"\n"); flush.console() } }
OUT$b<-br; sauve()
cc <- expand.grid(k=c(2,3), d=c(.5,1), n=c(10,30,100), cond=c("egal","sd13_n_egal","sd13_grand_n_sur_sd1","sd13_grand_n_sur_sd3"), stringsAsFactors=FALSE)
plan_c <- function(k,n,cond) { sdv <- if (cond=="egal") rep(1,k) else c(rep(1,k-1),3)
  nv <- switch(cond, egal=rep(n,k), sd13_n_egal=rep(n,k), sd13_grand_n_sur_sd1=c(rep(3*n,k-1),n), sd13_grand_n_sur_sd3=c(rep(n,k-1),3*n)); list(sd=sdv,n=nv) }
cr <- vector("list", nrow(cc)); OUT$c_cells<-cc
for (rep in 1:100) { for (ci in seq_len(nrow(cc))) { if (Sys.time()>FIN) break
    p<-plan_c(cc$k[ci],cc$n[ci],cc$cond[ci]); k<-cc$k[ci]; tau<-sqrt(mean(p$sd^2)); gj<-5000000L+ci*1000L+rep; set.seed(gj)
    mus<-c(rep(0,k-1),cc$d[ci]*tau); g<-factor(rep(1:k,p$n)); y<-rnorm(sum(p$n),rep(mus,p$n),rep(p$sd,p$n)); r<-het_safe(y,g,.5,500L,gj)
    cr[[ci]] <- rbind(cr[[ci]], c(rep=rep,het=r$log_bf,zs=ZS1(y,g),welch=oneway.test(y~g)$p.value,aov=anova(lm(y~g))[1,5],err=r$error_pct)) }
  if (Sys.time()>FIN) break
  if (rep %% 5==0) { OUT$c<-cr; sauve(); cat("c",rep,format(Sys.time()),"\n"); flush.console() } }
OUT$c<-cr; sauve(); cat("FIN 25b", format(Sys.time()), "\n")
