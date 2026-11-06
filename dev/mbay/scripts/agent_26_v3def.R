# Agent 26 : V3 f) degeneres, d) calibration sous l'a priori, e) sensibilite a r.
# Lance apres agent_24_v1v2.R dans le meme processus logique (jamais plus de 2 R).
source("agent_24_commun.R")
FIN <- as.POSIXct(paste(Sys.Date(), "14:15:40"))
OUT <- list(graines = "f: 6100+cas ; d: 7000000+cell*100000+rep ; e: a: agent_21 ; c: 8000000+rep")
sauve <- function() saveRDS(OUT, paste0(SORTIES, "agent_26_v3def.rds"))
# ---- f) ----
mk <- function(n, sd, seed, f = NULL) { set.seed(seed); y <- rnorm(sum(n), 0, rep(sd, n)); g <- factor(rep(seq_along(n), n)); if (!is.null(f)) y <- f(y, g); list(y = y, g = g) }
CAS <- list(
  n2_par_groupe = mk(c(2,2), c(1,1), 6101), n2_k4 = mk(c(2,2,2,2), rep(1,4), 6102),
  n_2_2_50 = mk(c(2,2,50), c(1,1,1), 6103), var_enorme = mk(c(10,10,10), c(1,1,1e4), 6104),
  var_enorme2 = mk(c(10,10), c(1e-3,1e3), 6105), quasi_cst = mk(c(10,10,10), c(1,1,1e-8), 6106),
  quasi_cst_2 = mk(c(10,10), c(1,1e-12), 6107), aberrante_45 = mk(c(15,15,15), c(1,1,1), 6108, function(y,g){y[1]<-y[1]+30;y}),
  k20_n5 = mk(rep(5,20), rep(1,20), 6109), k20_n30 = mk(rep(30,20), rep(1,20), 6110), k20_het = mk(rep(5,20), rep(c(1,5),10), 6111))
ftab <- list()
for (nm in names(CAS)) { C <- CAS[[nm]]; v <- list(); war <- character(0)
  t0 <- Sys.time()
  for (s in 1:10) v[[s]] <- withCallingHandlers(tryCatch(HET(C$y, C$g, .5, 4000L, s), error=function(e) list(log_bf=NA,error_pct=NA,ess=NA,n_draws=NA,low_ess=NA,err=conditionMessage(e))), warning=function(w){ war <<- c(war, conditionMessage(w)); invokeRestart("muffleWarning") })
  tt <- as.numeric(difftime(Sys.time(), t0, units="secs"))/10
  lb <- vapply(v, `[[`, 0, "log_bf"); ep <- vapply(v, `[[`, 0, "error_pct"); es <- vapply(v, `[[`, 0, "ess"); nd <- vapply(v, function(z) as.numeric(z$n_draws), 0)
  zs <- tryCatch(suppressWarnings(as.numeric(ZS(C$y, C$g, .5))), error=function(e) NA_real_)
  ftab[[nm]] <- data.frame(cas=nm, log_bf_moy=mean(lb), sd_logBF_pct=100*sd(lb), err_pct_moy=mean(ep), ess_min=min(es), n_draws_max=max(nd), n_alertes=length(war), sec_appel=tt, zs=zs,
    message=if (is.na(lb[1])) substr(v[[1]]$err,1,60) else if (length(war)) substr(war[1],1,60) else "")
  print(ftab[[nm]], digits=4); OUT$f <- ftab; sauve(); flush.console() }
OUT$f <- do.call(rbind, ftab); sauve()
# ---- d) ----
sim_prior <- function(n, R, r, seed) { set.seed(seed); k <- length(n); s <- exp(runif(k, -log(R), log(R)))
  H1 <- runif(1) < .5; a <- rep(0, k)
  if (H1) { g <- r^2/rchisq(1, 1); e <- rnorm(k); a <- sqrt(g*mean(s^2))*(e-mean(e)) }
  list(y = rnorm(sum(n), rep(a, n), rep(s, n)), g = factor(rep(seq_len(k), n)), H1 = H1) }
dc <- expand.grid(R=c(20,3), des=c("10,10","5,5,20"), stringsAsFactors=FALSE); dres <- vector("list", 4)
FIN_d <- as.POSIXct(paste(Sys.Date(), "14:14:00"))
for (rep in 1:1000) { for (ci in 1:4) { if (Sys.time()>FIN_d) break
    n <- as.numeric(strsplit(dc$des[ci],",")[[1]]); gj <- 7000000L+ci*100000L+rep; S <- sim_prior(n, dc$R[ci], .5, gj); r <- het_safe(S$y,S$g,.5,500L,gj)
    dres[[ci]] <- rbind(dres[[ci]], c(rep=rep, H1=S$H1, log_bf=r$log_bf, err=r$error_pct, ess=r$ess)) }
  if (Sys.time()>FIN_d) break
  if (rep %% 100==0) { OUT$d_cells<-dc; OUT$d<-dres; sauve(); cat("d",rep,format(Sys.time()),"\n"); flush.console() } }
OUT$d_cells<-dc; OUT$d<-dres; sauve()
# ---- e) ----
r21 <- readRDS(paste0(SORTIES,"agent_21_heterosc.rds")); n3 <- c(5,5,40); g3 <- factor(rep(1:3,n3))
eres <- list()
for (i in 1:300) { if (Sys.time()>FIN) break
  for (sc in c("A","B")) { sdv <- r21$plan[[paste0("sd_",sc)]]; gj <- r21[[sc]]$graine_jeu[i]; set.seed(gj); y <- rnorm(50,0,rep(sdv,n3))
    v <- vapply(c(.25,.5,1), function(r) het_safe(y,g3,r,500L,gj)$log_bf, 0); eres[[length(eres)+1]] <- c(sc=if(sc=="A")1 else 2, i=i, v) }
  # scenario c) : k=3, d=.5, n=30, sd=1,1,3 grand n sur sd1 ; k=2, d=1, n=30, egal
  gj <- 8000000L+i; set.seed(gj); n<-c(90,90,30); s<-c(1,1,3); m<-c(0,0,.5*sqrt(mean(s^2))); y<-rnorm(210,rep(m,n),rep(s,n)); g<-factor(rep(1:3,n))
  v <- vapply(c(.25,.5,1), function(r) het_safe(y,g,r,500L,gj)$log_bf, 0); eres[[length(eres)+1]] <- c(sc=3,i=i,v)
  gj <- 8100000L+i; set.seed(gj); y<-rnorm(60,rep(c(0,1),each=30),1); g<-factor(rep(1:2,each=30))
  v <- vapply(c(.25,.5,1), function(r) het_safe(y,g,r,500L,gj)$log_bf, 0); eres[[length(eres)+1]] <- c(sc=4,i=i,v)
  if (i %% 50==0) { OUT$e<-do.call(rbind,eres); sauve(); cat("e",i,format(Sys.time()),"\n"); flush.console() } }
OUT$e<-do.call(rbind,eres); colnames(OUT$e)<-c("sc","i","r025","r05","r1"); sauve()
cat("FIN A", format(Sys.time()), "\n")
