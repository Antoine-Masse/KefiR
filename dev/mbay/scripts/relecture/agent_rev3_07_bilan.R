# Bilan : ce qu'il dit sur des cas ou la bonne reponse est evidente. Sortie brute.
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
sink_print <- function(lab, expr) {
  cat("\n========", lab, "\n")
  t0 <- Sys.time()
  r <- tryCatch(expr, error=function(e) {cat("ERREUR:",conditionMessage(e),"\n"); NULL})
  cat(sprintf("[%.0f s]\n", as.numeric(difftime(Sys.time(),t0,units="secs"))))
  if (!is.null(r)) { print(r$models[,c("law","hypothesis","log_ml","posterior","rel_se","ess")],digits=4)
    if (!is.null(r$sensitivity)) {print(r$sensitivity$grid,digits=3); cat("by_u",r$sensitivity$by_u,"by_law",r$sensitivity$by_law,"range",r$sensitivity$range_lead,"verdict",r$sensitivity$verdict,"\n")}}
  flush.console(); invisible(r)
}
# 1. Valeur aberrante du test : le gaussien seul perd l'effet, les lois t le gardent.
set.seed(11); y <- c(rnorm(30), rnorm(30, 1.2)); y[1L] <- 40; g <- rep(1:2, each=30)
sink_print("aberrante, effet reel", m.bay(y,g,draws=2000L))
# 2. Effet enorme, propre.
set.seed(2); y <- c(rnorm(20),rnorm(20,10)); g <- rep(1:2,each=20)
sink_print("effet enorme", m.bay(y,g,draws=600L))
# 3. N grand, effet infime.
set.seed(3); y <- c(rnorm(5000),rnorm(5000,.08)); g <- rep(1:2,each=5000)
sink_print("N=1e4, ecart .08 sd", m.bay(y,g,draws=500L,sensitivity=FALSE))
# 4. Rien : H0 vraie, n moyen.
set.seed(4); y <- rnorm(60); g <- rep(1:3,each=20)
sink_print("H0 vraie", m.bay(y,g,draws=600L))
# 5. Plateau k=8, sous forme bilan (compare a Welch p).
set.seed(7); y <- c(rnorm(70), rnorm(3,1000,1)); g <- rep(1:8,c(rep(10,7),3))
print(oneway.test(y~factor(g))$p.value)
sink_print("odds prior_h1=.2", m.bay(c(rnorm(20),rnorm(20,1)),rep(1:2,each=20),draws=300L,prior_h1=.2))
