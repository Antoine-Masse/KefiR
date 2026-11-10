# agent_60 : calibration des comptages (D17, voie 2). Par scenario (3 groupes de n, moyennes
# vraies mu, binomiale negative de parametre theta, ou Poisson si theta = Inf) : part des jeux
# ou m.bay dit "difference" (P > .5, > .75, > .95), a cote de Poisson seul (meme
# approximation BIC, sans la loi binomiale negative) et du test du rapport de vraisemblance
# de glm.nb (p < .05) ; part des jeux ou la binomiale negative est en tete.
# Usage : Rscript agent_60_comptages.R N
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 300L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
SCEN <- list(list(nom = "H0_poisson_mu5_n20", mu = c(5, 5, 5), th = Inf, n = 20),
	list(nom = "H0_surdisp_mu5_th1_n20", mu = c(5, 5, 5), th = 1, n = 20),
	list(nom = "H0_surdisp_mu3_th0.5_n15", mu = c(3, 3, 3), th = 0.5, n = 15),
	list(nom = "H1_surdisp_mu4-4-8_th2_n20", mu = c(4, 4, 8), th = 2, n = 20),
	list(nom = "H1_poisson_mu4-4-6_n20", mu = c(4, 4, 6), th = Inf, n = 20))
out <- c(R.version.string, sprintf("agent_60 : comptages, N = %d jeux par scenario, prior_h1 = 0.5", N),
	"scenario : m.bay P>.5 / >.75 / >.95 || Poisson seul P>.5 / >.75 / >.95 || glm.nb LRT p<.05 || BN en tete")
for (s in SCEN) {
	r <- t(vapply(seq_len(N), function(i) {
		set.seed(600000L + 7L * nchar(s$nom) + i)
		g <- factor(rep(1:3, each = s$n)); mu <- s$mu[g]
		y <- if (is.finite(s$th)) MASS::rnegbin(length(mu), mu = mu, theta = s$th) else rpois(length(mu), mu)
		f0 <- .bay_counts_fit(y, factor(rep(1L, length(y)))); f1 <- .bay_counts_fit(y, g)
		po <- .bay_counts_post(f0, f1, 0.5)
		pp <- 1 / (1 + exp(-(f1$lm_pois - f0$lm_pois)))
		lrt <- tryCatch({
			m1 <- suppressWarnings(MASS::glm.nb(y ~ g)); m0 <- suppressWarnings(MASS::glm.nb(y ~ 1))
			stats::pchisq(2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0))), 2, lower.tail = FALSE)
		}, error = function(e) NA_real_)
		c(p = po$p_h1, pp = pp, lrt = lrt, nb = po$p_nb)
	}, numeric(4)))
	f <- function(x) sprintf("%.3f", mean(x, na.rm = TRUE))
	out <- c(out, sprintf("%s : %s / %s / %s || %s / %s / %s || %s || %s", s$nom, f(r[, "p"] > .5), f(r[, "p"] > .75), f(r[, "p"] > .95),
		f(r[, "pp"] > .5), f(r[, "pp"] > .75), f(r[, "pp"] > .95), f(r[, "lrt"] < .05), f(r[, "nb"] > .5)))
	writeLines(out, file.path(BASE, "dev/mbay/scripts/resultats/agent_60.txt"))
}
