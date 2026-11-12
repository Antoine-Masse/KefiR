# agent_59 : calibration du tableau de contingence (D17, voie 1), a cote du khi-deux.
# Par scenario (k groupes de n, proportions vraies p) : part des jeux ou m.bay dit "les
# proportions different" (P > .5, > .75, > .95), khi-deux p < .05 (simule si effectifs
# attendus < 5), couverture de l'intervalle a 95 % de la difference de proportions.
# Usage : Rscript agent_59_contingence.R N
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 500L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
SCEN <- list(list(nom = "k2_n20_H0_p.3", p = c(.3, .3), n = 20), list(nom = "k2_n50_H0_p.5", p = c(.5, .5), n = 50),
	list(nom = "k3_n30_H0_p.2", p = c(.2, .2, .2), n = 30), list(nom = "k2_n40_ecart_.2", p = c(.3, .5), n = 40),
	list(nom = "k3_n30_un_ecart", p = c(.5, .5, .8), n = 30), list(nom = "k2_n10_H0_petit", p = c(.2, .2), n = 10))
out <- c(R.version.string, sprintf("agent_59 : contingence, N = %d jeux par scenario, a = 1, prior_h1 = 0.5", N),
	"scenario : P(diff)>.5 / >.75 / >.95 | khi-deux p<.05 || couverture diff. de proportions (paires)")
for (s in SCEN) {
	r <- t(vapply(seq_len(N), function(i) {
		set.seed(590000L + 13L * nchar(s$nom) + i)
		k <- length(s$p)
		g <- factor(rep(seq_len(k), each = s$n)); y <- factor(rbinom(k * s$n, 1, s$p[g]), levels = 0:1)
		tab <- table(g, y)
		lb <- .bay_cont_bf_pkg(tab, 1); p <- 1 / (1 + exp(-lb))
		ex <- suppressWarnings(stats::chisq.test(tab)$expected)
		pc <- if (any(ex < 5)) stats::chisq.test(tab, simulate.p.value = TRUE, B = 2000)$p.value else suppressWarnings(stats::chisq.test(tab, correct = FALSE)$p.value)
		ph <- .bay_cont_posthoc(tab, 1, 0.5, 2000L, i)$pairs
		tv <- s$p[as.integer(ph$group2)] - s$p[as.integer(ph$group1)]
		c(p = p, pc = pc, cover = mean(ph$diff_low <= tv & tv <= ph$diff_high))
	}, numeric(3)))
	f <- function(x) sprintf("%.3f", mean(x, na.rm = TRUE))
	out <- c(out, sprintf("%s : %s / %s / %s | %s || %s", s$nom, f(r[, "p"] > .5), f(r[, "p"] > .75), f(r[, "p"] > .95), f(r[, "pc"] < .05), f(r[, "cover"])))
	writeLines(out, file.path(BASE, "dev/mbay/scripts/resultats/agent_59.txt"))
}
