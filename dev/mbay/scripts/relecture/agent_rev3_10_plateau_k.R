# Plateau en k : un groupe de 3 mesures a 1000 sd, les autres identiques. Le BF passe-t-il sous 1 ?
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
for (k in c(12,20,30)) for (ns in c(3,5)) for (d in c(10,1e3)) {
  des <- c(rep(10,k-1),ns); g <- rep(seq_len(k),des)
  set.seed(7); y <- c(rnorm(sum(des[-k])), rnorm(ns,d,1))
  t0 <- Sys.time()
  a <- tryCatch(ans_of(y,g,draws=1000L,laws=c("normal_eq","normal_het")),error=function(e) {cat("ERREUR",conditionMessage(e),"\n");NULL})
  if(!is.null(a)) show(sprintf("k=%d n_s=%d d=%g Welch p=%.0e (%.0fs)",k,ns,d,oneway.test(y~factor(g))$p.value,as.numeric(difftime(Sys.time(),t0,units="secs"))),a)
  flush.console()
}
cat("--- Toutes les lois, k=12\n")
des <- c(rep(10,11),3); g <- rep(seq_len(12),des); set.seed(7); y <- c(rnorm(110), rnorm(3,1000,1))
t0 <- Sys.time(); a <- ans_of(y,g,draws=800L); show(sprintf("k=12 n_s=3 d=1000 (%.0fs)",as.numeric(difftime(Sys.time(),t0,units="secs"))),a)
print(a$models[,c("law","hypothesis","log_ml","posterior","rel_se","ess")],digits=4)
cat("--- Fixed by group : plan cible ou m.bay imprime\n")
r <- m.bay(y,g,draws=800L,sensitivity=FALSE)
