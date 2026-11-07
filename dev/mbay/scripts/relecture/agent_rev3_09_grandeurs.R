# Temps et comportement a N grand, k grand. Tirages reduits, sans sensibilite sauf mention.
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
run <- function(lab, y, g, draws=1000L, laws=.bay_law_names()) {
  t0 <- Sys.time()
  r <- tryCatch(setTimeLimit(elapsed=280,transient=TRUE),error=function(e) NULL)
  a <- tryCatch(ans_of(y,g,draws=draws,laws=laws),error=function(e) {cat(lab,"ERREUR",conditionMessage(e),"\n");NULL})
  setTimeLimit()
  el <- as.numeric(difftime(Sys.time(),t0,units="secs"))
  if(!is.null(a)) show(sprintf("%s (%.0fs)",lab,el),a) else cat(sprintf("%s arret apres %.0fs\n",lab,el))
  flush.console()
}
set.seed(1)
y <- c(rnorm(2e4),rnorm(2e4,.05)); g <- rep(1:2,each=2e4)
run("N=4e4 k=2, normales", y,g,laws=c("normal_eq","normal_het"))
run("N=4e4 k=2, toutes", y,g,draws=500L)
y <- rnorm(2e5); g <- rep(1:2,each=1e5)
run("N=2e5 k=2, normales", y,g,laws=c("normal_eq","normal_het"))
run("N=2e5 k=2, toutes", y,g,draws=500L)
set.seed(2); for (k in c(15,30,50)) { g <- rep(seq_len(k),each=6); y <- rnorm(length(g))
  run(sprintf("k=%d n=6 normales",k),y,g,laws=c("normal_eq","normal_het"))
  run(sprintf("k=%d n=6 toutes",k),y,g,draws=500L)
}
