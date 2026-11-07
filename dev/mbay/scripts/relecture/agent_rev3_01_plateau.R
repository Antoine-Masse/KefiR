# Conjecture : BF moyenne ~ max(cte, d^(n_s + 2u - k)). Plateau des que k >= n_s + 2u.
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
plans <- list(c(rep(10,7),3), c(rep(10,9),5), c(rep(10,7),2), c(rep(10,5),1))
for (des in plans) {
  k <- length(des); ns <- des[k]
  cat("\nPlan", paste(des,collapse=","), " k=",k," n_s=",ns," exposant prevu ",ns+3-k,"\n")
  for (d in c(5,50,1e3,1e5)) {
    set.seed(7)
    y <- c(rnorm(sum(des[-k])), rnorm(ns,d,1))
    g <- rep(seq_len(k),des)
    t0 <- Sys.time()
    a <- tryCatch(ans_of(y,g,draws=1000L),error=function(e) {cat("ERREUR",conditionMessage(e),"\n");NULL})
    if (!is.null(a)) show(sprintf("d=%g (%.0fs, Welch p=%.1e, zs=%.0f)",d,as.numeric(difftime(Sys.time(),t0,units="secs")),oneway.test(y~factor(g))$p.value,.bay_bf_oneway_zs(y,g)),a)
    flush.console()
  }
}
