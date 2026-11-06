# Agent 25 : V3 a) (appariement agent_21), effet de n_draws, b) homoscedastique, c) puissance
source("agent_24_commun.R")
FIN <- as.POSIXct(paste(Sys.Date(), "14:14:15"))
rds21 <- readRDS(paste0(SORTIES, "agent_21_heterosc.rds"))
n3 <- c(5,5,40); g3 <- factor(rep(1:3, n3)); ZS1 <- function(y,g) tryCatch(suppressWarnings(as.numeric(ZS(y,g,.5))), error=function(e) NA_real_)
OUT <- list(graines="a: graine_jeu de agent_21 ; hetero seed=graine_jeu ; b: 3000000+cell*1000+rep ; c: 5000000+cell*1000+rep")
sauve <- function() saveRDS(OUT, paste0(SORTIES, "agent_25_v3abc.rds"))
tgen <- function(sd, gj) { set.seed(gj); rnorm(50, 0, rep(sd, n3)) }
# Effet de n_draws
nd <- list()
for (sc in c("A","B")) { sdv <- rds21$plan[[paste0("sd_",sc)]]; d <- rds21[[sc]]
  for (i in 1:25) { y <- tgen(sdv, d$graine_jeu[i])
    v <- vapply(c(500L,2000L,4000L,20000L), function(m) het_safe(y,g3,.5,m,d$graine_jeu[i])$log_bf, 0)
    nd[[length(nd)+1]] <- c(sc=match(sc,c("A","B")), jeu=i, v) } }
OUT$ndraws <- do.call(rbind, nd); colnames(OUT$ndraws) <- c("sc","jeu","nd500","nd2000","nd4000","nd20000"); sauve()
cat("ndraws fini", format(Sys.time()), "\n")
# a)
for (sc in c("A","B")) { sdv <- rds21$plan[[paste0("sd_",sc)]]; d <- rds21[[sc]][1:1000,]
  res <- matrix(NA_real_, 1000, 10, dimnames=list(NULL,c("log_bf","err","ess","zs","p_aov","p_welch","p_kw","t","gj","low_ess")))
  for (i in 1:1000) { y <- tgen(sdv, d$graine_jeu[i]); t0 <- Sys.time()
    r <- het_safe(y,g3,.5,2000L,d$graine_jeu[i])
    res[i,] <- c(r$log_bf,r$error_pct,r$ess,ZS1(y,g3),anova(lm(y~g3))[1,5],oneway.test(y~g3)$p.value,kruskal.test(y~g3)$p.value,as.numeric(Sys.time()-t0),d$graine_jeu[i],as.numeric(isTRUE(r$low_ess)))
    if (i %% 200==0) { OUT[[paste0("a_",sc)]] <- res; sauve(); cat(sc,i,format(Sys.time()),"\n") } }
  OUT[[paste0("a_",sc)]] <- res; OUT[[paste0("a_",sc,"_ref")]] <- d; sauve() }
# b)
cells <- expand.grid(k=c(2,3,5), n=c(5,20,100)); br <- vector("list", nrow(cells))
for (rep in 1:100) { for (ci in seq_len(nrow(cells))) { k<-cells$k[ci]; n<-cells$n[ci]; gj<-3000000L+ci*1000L+rep; set.seed(gj)
    g<-factor(rep(1:k,each=n)); y<-rnorm(k*n); r<-het_safe(y,g,.5,500L,gj)
    br[[ci]] <- rbind(br[[ci]], c(rep=rep,het=r$log_bf,zs=ZS1(y,g),err=r$error_pct,ess=r$ess,aov=anova(lm(y~g))[1,5],welch=oneway.test(y~g)$p.value)) }
  if (rep %% 25==0) { OUT$b_cells<-cells; OUT$b<-br; sauve(); cat("b",rep,format(Sys.time()),"\n") } }
OUT$b_cells<-cells; OUT$b<-br; sauve()
# c)
cc <- expand.grid(k=c(2,3), d=c(.5,1), n=c(10,30,100), cond=c("egal","sd13_n_egal","sd13_grand_n_sur_sd1","sd13_grand_n_sur_sd3"), stringsAsFactors=FALSE)
plan_c <- function(k,n,cond) { sdv <- if (cond=="egal") rep(1,k) else c(rep(1,k-1),3)
  nv <- switch(cond, egal=rep(n,k), sd13_n_egal=rep(n,k), sd13_grand_n_sur_sd1=c(rep(3*n,k-1),n), sd13_grand_n_sur_sd3=c(rep(n,k-1),3*n)); list(sd=sdv,n=nv) }
cr <- vector("list", nrow(cc)); OUT$c_cells<-cc
for (rep in 1:100) { for (ci in seq_len(nrow(cc))) { if (Sys.time()>FIN) break
    p<-plan_c(cc$k[ci],cc$n[ci],cc$cond[ci]); k<-cc$k[ci]; tau<-sqrt(mean(p$sd^2)); gj<-5000000L+ci*1000L+rep; set.seed(gj)
    mus<-c(rep(0,k-1),cc$d[ci]*tau); g<-factor(rep(1:k,p$n)); y<-rnorm(sum(p$n),rep(mus,p$n),rep(p$sd,p$n)); r<-het_safe(y,g,.5,500L,gj)
    cr[[ci]] <- rbind(cr[[ci]], c(rep=rep,het=r$log_bf,zs=ZS1(y,g),welch=oneway.test(y~g)$p.value,aov=anova(lm(y~g))[1,5],err=r$error_pct)) }
  if (Sys.time()>FIN) break
  if (rep %% 10==0) { OUT$c<-cr; sauve(); cat("c",rep,format(Sys.time()),"\n") } }
OUT$c<-cr; sauve(); cat("FIN B", format(Sys.time()), "\n")
