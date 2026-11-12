# agent_56 : calibration des donnees appariees et des mesures repetees (D15).
# Par scenario (n sujets, k conditions, effet vrai en ecarts-types residuels, ecart-type
# sujet 2) : part des jeux ou m.bay dit "difference" (P > .5, > .75, > .95), a cote du t
# apparie (k = 2) ou de aov(y ~ temps + Error(id/temps)) (k >= 3) a p < .05 ; couverture de
# l'intervalle de dz (k = 2) ou de delta (k >= 3) ; k >= 3 : paires egales faussement dites
# differentes (P > .75), regroupement vrai en tete.
# Usage : Rscript agent_56_repete.R N coeurs
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 200L
nc <- if (length(args) >= 2L) as.integer(args[2L]) else 4L
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
.libPaths(c(LIB, .libPaths()))
library(parallel)
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
cl <- makeCluster(nc)
invisible(clusterCall(cl, function(LIB, BASE) {
	.libPaths(c(LIB, .libPaths()))
	suppressMessages(devtools::load_all(BASE, quiet = TRUE))
	NULL
}, LIB, BASE))
SCEN <- list(
	list(nom = "k2_n12_H0", mu = c(0, 0), n = 12), list(nom = "k2_n20_H0", mu = c(0, 0), n = 20),
	list(nom = "k2_n15_d0.7", mu = c(0, 0.7), n = 15),
	list(nom = "k3_n12_H0", mu = c(0, 0, 0), n = 12), list(nom = "k3_n15_un_ecart", mu = c(0, 0, 0.8), n = 15),
	list(nom = "k4_n12_graduel", mu = c(0, 0.3, 0.6, 0.9), n = 12))
un <- function(i, s) {
	set.seed(560000L + 101L * nchar(s$nom) + 7L * s$n + i)
	k <- length(s$mu); n <- s$n
	d <- data.frame(id = factor(rep(seq_len(n), k)), temps = factor(rep(paste0("t", seq_len(k)), each = n)))
	d$y <- rnorm(n, 0, 2)[d$id] + s$mu[d$temps] + rnorm(n * k)
	r <- suppressMessages(m.bay(y ~ temps | id, data = d, verbose = FALSE, draws = 2000L, sensitivity = FALSE, seed = i))
	pr <- r$posthoc$pairs
	if (k == 2L) {
		p_cl <- t.test(d$y[d$temps == "t2"], d$y[d$temps == "t1"], paired = TRUE)$p.value
		dd <- d$y[d$temps == "t2"] - d$y[d$temps == "t1"]
		true_dz <- (s$mu[2] - s$mu[1]) / sqrt(2)
		cover <- as.numeric(pr$delta_low <= true_dz & true_dz <= pr$delta_high)
		return(c(p = r$prob_effect, pcl = p_cl, cover = cover, eqfalse = NA, top = NA))
	}
	a <- summary(aov(y ~ temps + Error(id / temps), data = d))
	p_cl <- a[["Error: id:temps"]][[1]][["Pr(>F)"]][1]
	tv <- s$mu[as.integer(sub("t", "", pr$group2))] - s$mu[as.integer(sub("t", "", pr$group1))]
	cover <- mean(pr$delta_low <= tv & tv <= pr$delta_high)
	eq <- tv == 0
	eqfalse <- if (any(eq)) as.numeric(any(pr$p_differ[eq] > 0.75)) else NA
	code_true <- paste(as.integer(factor(s$mu, levels = unique(s$mu))), collapse = "")
	c(p = r$prob_effect, pcl = p_cl, cover = cover, eqfalse = eqfalse, top = as.numeric(r$posthoc$partitions$code[1L] == code_true))
}
clusterExport(cl, "un")
out <- c(R.version.string, sprintf("agent_56 : appariees et mesures repetees, N = %d jeux par scenario, sujet sd 2, residu sd 1, r = 0.5, prior_h1 = 0.5", N),
	"scenario : P(diff)>.5 / >.75 / >.95 | classique p<.05 || couverture delta (dz si k=2) || paire egale dite differente (P>.75) || regroupement vrai en tete")
for (s in SCEN) {
	r <- do.call(rbind, parLapplyLB(cl, seq_len(N), un, s = s))
	f <- function(x) if (all(is.na(x))) "-" else sprintf("%.3f", mean(x, na.rm = TRUE))
	out <- c(out, sprintf("%s : %s / %s / %s | %s || %s || %s || %s", s$nom, f(r[, "p"] > .5), f(r[, "p"] > .75), f(r[, "p"] > .95),
		f(r[, "pcl"] < .05), f(r[, "cover"]), f(r[, "eqfalse"]), f(r[, "top"])))
	writeLines(out, file.path(BASE, "dev/mbay/scripts/resultats/agent_56.txt"))
}
stopCluster(cl)
