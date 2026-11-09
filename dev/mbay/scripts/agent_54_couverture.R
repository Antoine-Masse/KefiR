# agent_54 : couverture des intervalles de credibilite du post-hoc (D11), toutes les paires.
# Les intervalles sont conditionnels au modele "tous les groupes peuvent differer" ; on mesure,
# par scenario, la part des jeux ou l'intervalle a 95 % de l'ecart (diff) et de delta contient
# la vraie valeur, separement pour les paires vraiment egales (vraie valeur 0) et differentes.
# Usage : Rscript agent_54_couverture.R N coeurs
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 300L
nc <- if (length(args) >= 2L) as.integer(args[2L]) else 4L
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
.libPaths(c(LIB, .libPaths()))
library(parallel)
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
cl <- makeCluster(nc)
invisible(clusterCall(cl, function(LIB, BASE) {
	.libPaths(c(LIB, .libPaths()))
	source(file.path(BASE, "dev/mbay/scripts/agent_load_bay.R"))
	.kefir_bay(BASE)
	suppressMessages({library(BayesFactor); library(bayestestR)})
	NULL
}, LIB, BASE))
SCEN <- list(
	list(nom = "k2_n10_d0.5", mu = c(0, 0.5), n = 10), list(nom = "k3_n10_H0", mu = c(0, 0, 0), n = 10),
	list(nom = "k3_n15_un_ecart", mu = c(0, 0, 1), n = 15), list(nom = "k4_n20_deux_blocs", mu = c(0, 0, 0.8, 0.8), n = 20),
	list(nom = "k4_n10_graduel", mu = c(0, 0.3, 0.6, 0.9), n = 10))
un <- function(i, s) {
	set.seed(540000L + 977L * nchar(s$nom) + 31L * s$n + i)
	k <- length(s$mu)
	g <- factor(rep(letters[seq_len(k)], each = s$n))
	y <- s$mu[as.integer(g)] + stats::rnorm(length(g))
	pr <- .bay_posthoc_oneway(y, g, 0.5, 0.5, 4000L, i)$pairs
	tv <- s$mu[match(pr$group2, letters)] - s$mu[match(pr$group1, letters)]
	c(cd = as.numeric(pr$diff_low <= tv & tv <= pr$diff_high), cs = as.numeric(pr$delta_low <= tv & tv <= pr$delta_high))
}
clusterExport(cl, "un")
out <- c(R.version.string, sprintf("agent_54 : couverture des intervalles a 95 %% du post-hoc, N = %d jeux par scenario, sigma = 1, r = 0.5", N),
	"scenario : paires egales -> couverture ecart / delta (nb de paires) || paires differentes -> couverture ecart / delta (nb)")
for (s in SCEN) {
	r <- do.call(rbind, parLapplyLB(cl, seq_len(N), un, s = s))
	k <- length(s$mu); np <- k * (k - 1) / 2; pc <- utils::combn(k, 2L)
	eq <- s$mu[pc[1L, ]] == s$mu[pc[2L, ]]
	cd <- r[, seq_len(np), drop = FALSE]; cs <- r[, np + seq_len(np), drop = FALSE]
	f <- function(m, w) if (any(w)) sprintf("%.3f", mean(m[, w])) else "-"
	out <- c(out, sprintf("%s : %s / %s (%d) || %s / %s (%d)", s$nom, f(cd, eq), f(cs, eq), sum(eq), f(cd, !eq), f(cs, !eq), sum(!eq)))
	writeLines(out, file.path(BASE, "dev/mbay/scripts/resultats/agent_54.txt"))
}
stopCluster(cl)
