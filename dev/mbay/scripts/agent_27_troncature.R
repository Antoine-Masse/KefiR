# Agent 27 : ecart lie a la troncature. Memes jeux que agent_26 d) (memes graines) ;
# BF de reference sous a priori tronque : grille deterministe de agent_24_v1v2.R, bornes = +/- log R.
if (Sys.getenv("AG27") != "1") quit("no")
source("agent_24_commun.R")
ex <- parse("agent_24_v1v2.R"); for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]])[1] %in% c("lse","stats_of","grille_t","route_A")) eval(e)
FIN <- as.POSIXct(paste(Sys.Date(), "14:15:50"))
sim_prior <- function(n, R, r, seed) { set.seed(seed); k <- length(n); s <- exp(runif(k, -log(R), log(R)))
  H1 <- runif(1) < .5; a <- rep(0, k)
  if (H1) { g <- r^2/rchisq(1, 1); e <- rnorm(k); a <- sqrt(g*mean(s^2))*(e-mean(e)) }
  list(y = rnorm(sum(n), rep(a, n), rep(s, n)), g = factor(rep(seq_len(k), n)), H1 = H1) }
dc <- expand.grid(R=c(20,3), des=c("10,10","5,5,20"), stringsAsFactors=FALSE); res <- vector("list",4)
d <- readRDS(paste0(SORTIES,"agent_26_v3def.rds"))$d
for (rep in 1:400) { for (ci in c(1,2,3,4)) { if (Sys.time()>FIN) break
    if (ci==3 && rep>60) next
    n <- as.numeric(strsplit(dc$des[ci],",")[[1]]); gj <- 7000000L+ci*100000L+rep; S <- sim_prior(n, dc$R[ci], .5, gj); R<-dc$R[ci]
    a <- tryCatch(route_A(S$y,S$g,.5,1,c(-log(R),log(R))), error=function(e) list(lbf=NA))
    res[[ci]] <- rbind(res[[ci]], c(rep=rep,H1=S$H1,lbf_tr=a$lbf,lbf_fn=d[[ci]][rep,"log_bf"])) }
  if (Sys.time()>FIN) break
  if (rep %% 20==0) { saveRDS(list(cells=dc,res=res), paste0(SORTIES,"agent_27_troncature.rds")); cat(rep, format(Sys.time()),"\n"); flush.console() } }
saveRDS(list(cells=dc,res=res), paste0(SORTIES,"agent_27_troncature.rds")); cat("FIN 27", format(Sys.time()), "\n")
