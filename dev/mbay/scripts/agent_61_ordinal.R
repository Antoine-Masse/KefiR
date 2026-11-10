# agent_61 : calibration de la voie ordinale (D17, voie 3). Likert a K niveaux par seuils
# d'une variable latente normale ; 3 groupes de n ; decalages latents d. Part des jeux ou
# m.bay dit "difference" (P > .5, > .75, > .95), a cote du test du rapport de vraisemblance
# de polr et de kruskal.test (p < .05). Un scenario a cotes non proportionnelles : meme
# moyenne latente, dispersion latente differente dans un groupe.
# Usage : Rscript agent_61_ordinal.R N
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 300L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
SCEN <- list(
	list(nom = "K5_n10_H0", K = 5, n = 10, d = c(0, 0, 0), s = c(1, 1, 1)),
	list(nom = "K5_n30_H0", K = 5, n = 30, d = c(0, 0, 0), s = c(1, 1, 1)),
	list(nom = "K7_n10_H0", K = 7, n = 10, d = c(0, 0, 0), s = c(1, 1, 1)),
	list(nom = "K7_n30_H0", K = 7, n = 30, d = c(0, 0, 0), s = c(1, 1, 1)),
	list(nom = "K5_n10_effet0.8", K = 5, n = 10, d = c(0, 0, 0.8), s = c(1, 1, 1)),
	list(nom = "K5_n30_effet0.5", K = 5, n = 30, d = c(0, 0, 0.5), s = c(1, 1, 1)),
	list(nom = "K7_n30_effet0.5", K = 7, n = 30, d = c(0, 0, 0.5), s = c(1, 1, 1)),
	list(nom = "K5_n30_non_prop", K = 5, n = 30, d = c(0, 0, 0), s = c(1, 1, 2.5)))
out <- c(R.version.string, sprintf("agent_61 : voie ordinale, N = %d jeux par scenario, 3 groupes, prior_h1 = 0.5", N),
	"scenario : m.bay P>.5 / >.75 / >.95 || polr LRT p<.05 || kruskal p<.05 || voie choisie ordinale")
for (s in SCEN) {
	cuts <- stats::qnorm(seq_len(s$K - 1L) / s$K)
	r <- t(vapply(seq_len(N), function(i) {
		set.seed(610000L + 13L * nchar(s$nom) + i)
		g <- factor(rep(c("a", "b", "c"), each = s$n))
		z <- s$d[as.integer(g)] + s$s[as.integer(g)] * stats::rnorm(length(g))
		y <- findInterval(z, cuts) + 1L
		rr <- suppressWarnings(m.bay(y, g, verbose = FALSE, sensitivity = FALSE, posthoc = FALSE, response = "ordinal"))
		auto <- identical(.bay_int_route(y, g), "ordinal")
		yo <- factor(y, ordered = TRUE)
		lrt <- tryCatch({
			m1 <- suppressWarnings(MASS::polr(yo ~ g)); m0 <- suppressWarnings(MASS::polr(yo ~ 1))
			stats::pchisq(2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0))), 2, lower.tail = FALSE)
		}, error = function(e) NA_real_)
		c(p = rr$prob_effect, lrt = lrt, kw = stats::kruskal.test(y, g)$p.value, auto = auto)
	}, numeric(4)))
	f <- function(x) sprintf("%.3f", mean(x, na.rm = TRUE))
	out <- c(out, sprintf("%s : %s / %s / %s || %s || %s || %s", s$nom, f(r[, "p"] > .5), f(r[, "p"] > .75), f(r[, "p"] > .95),
		f(r[, "lrt"] < .05), f(r[, "kw"] < .05), f(r[, "auto"])))
	writeLines(out, file.path(BASE, "dev/mbay/scripts/resultats/agent_61.txt"))
}
