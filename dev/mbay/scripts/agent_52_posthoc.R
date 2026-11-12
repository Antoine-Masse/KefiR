# agent_52 : post-hoc de m.bay a un facteur (D11), a cote de Tukey HSD.
# Par scenario (k groupes, n par groupe, moyennes vraies en ecarts-types, sigma = 1) :
#   - paires vraiment egales : part des jeux ou au moins une est dite "differente"
#     (P > .5, > .75, > .95) ; Tukey : au moins une p < .05 ;
#   - paires vraiment differentes : puissance moyenne P > .5, P > .75 ; Tukey p < .05 ;
#   - couverture de l'intervalle a 95 % de delta sur les paires vraiment differentes ;
#   - regroupement le plus probable = regroupement vrai.
# Calcul par .bay_posthoc_oneway de KefiR. Reprise possible : resultats/agent_52_cas/.
# Usage : Rscript agent_52_posthoc.R N coeurs
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 300L
nc <- if (length(args) >= 2L) as.integer(args[2L]) else 4L
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
.libPaths(c(LIB, .libPaths()))
library(parallel)
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
CAS <- file.path(BASE, "dev/mbay/scripts/resultats/agent_52_cas")
dir.create(CAS, showWarnings = FALSE, recursive = TRUE)
cl <- makeCluster(nc)
invisible(clusterCall(cl, function(LIB, BASE) {
	.libPaths(c(LIB, .libPaths()))
	source(file.path(BASE, "dev/mbay/scripts/agent_load_bay.R"))
	.kefir_bay(BASE)
	suppressMessages({library(BayesFactor); library(bayestestR)})
	NULL
}, LIB, BASE))

SCEN <- list(
	list(nom = "k3_n10_H0", mu = c(0, 0, 0), n = 10), list(nom = "k3_n20_H0", mu = c(0, 0, 0), n = 20),
	list(nom = "k3_n15_un_ecart", mu = c(0, 0, 1), n = 15), list(nom = "k3_n15_trois", mu = c(0, 0.6, 1.2), n = 15),
	list(nom = "k4_n10_H0", mu = c(0, 0, 0, 0), n = 10), list(nom = "k4_n15_un_ecart", mu = c(0, 0, 0, 1), n = 15),
	list(nom = "k4_n15_deux_blocs", mu = c(0, 0, 0.8, 0.8), n = 15), list(nom = "k5_n10_H0", mu = rep(0, 5), n = 10))

un <- function(i, s) {
	set.seed(520000L + 1000L * nchar(s$nom) + 17L * s$n + i)
	k <- length(s$mu)
	g <- factor(rep(letters[seq_len(k)], each = s$n))
	y <- s$mu[as.integer(g)] + stats::rnorm(length(g))
	ph <- .bay_posthoc_oneway(y, g, 0.5, 0.5, 2000L, i)
	pr <- ph$pairs
	tk <- stats::TukeyHSD(stats::aov(y ~ g))$g[, "p adj"]
	d1 <- s$mu[match(pr$group1, letters)]; d2 <- s$mu[match(pr$group2, letters)]
	true_delta <- d2 - d1
	code_true <- paste(as.integer(factor(s$mu, levels = unique(s$mu))), collapse = "")
	c(p = pr$p_differ, tk = unname(tk), cover = as.numeric(pr$delta_low <= true_delta & true_delta <= pr$delta_high),
		top_ok = as.numeric(identical(ph$partitions$code[1L], code_true)))
}
clusterExport(cl, "un")

lignes <- c(R.version.string, sprintf("agent_52 : post-hoc a un facteur, N = %d jeux par scenario, r = 0.5, prior_h1 = 0.5, Ewens alpha = 1", N),
	"scenario : paires egales -> un faux P>.5 / P>.75 / P>.95 | Tukey un faux p<.05 || paires differentes -> P>.5 / P>.75 | Tukey p<.05 | couverture HDI delta || regroupement vrai en tete")
for (s in SCEN) {
	r <- NULL
	for (d0 in seq(1L, N, by = 50L)) {
		fic <- file.path(CAS, sprintf("%s_%05d.rds", s$nom, d0))
		if (!file.exists(fic)) saveRDS(do.call(rbind, parLapplyLB(cl, d0:min(d0 + 49L, N), un, s = s)), fic)
		r <- rbind(r, readRDS(fic))
	}
	k <- length(s$mu); np <- k * (k - 1) / 2
	pc <- utils::combn(k, 2L)
	eq <- s$mu[pc[1L, ]] == s$mu[pc[2L, ]]
	P <- r[, seq_len(np), drop = FALSE]; TK <- r[, np + seq_len(np), drop = FALSE]; CV <- r[, 2 * np + seq_len(np), drop = FALSE]
	faux <- function(th) if (any(eq)) mean(apply(P[, eq, drop = FALSE] > th, 1L, any)) else NA
	lignes <- c(lignes, sprintf("%s (N=%d) : %s / %s / %s | %s || %s / %s | %s | %s || %.3f", s$nom, nrow(r),
		format(round(faux(.5), 3)), format(round(faux(.75), 3)), format(round(faux(.95), 3)),
		if (any(eq)) format(round(mean(apply(TK[, eq, drop = FALSE] < .05, 1L, any)), 3)) else "-",
		if (any(!eq)) format(round(mean(P[, !eq] > .5), 3)) else "-", if (any(!eq)) format(round(mean(P[, !eq] > .75), 3)) else "-",
		if (any(!eq)) format(round(mean(TK[, !eq] < .05), 3)) else "-",
		if (any(!eq)) format(round(mean(CV[, !eq]), 3)) else "-", mean(r[, "top_ok"])))
	writeLines(lignes, file.path(BASE, "dev/mbay/scripts/resultats/agent_52.txt"))
}
stopCluster(cl)
