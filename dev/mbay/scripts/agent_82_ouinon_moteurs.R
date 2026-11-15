# agent_82 : oui/non a un facteur (I4 de l'expert, passe 3). Sur les memes tables, la forme close
# a marges neutres (a = 2, voie a un facteur, commune avec g.bay) et le facteur de Bayes de
# deviance a g de Zellner-Siow (r = sqrt(2)/4, moteur des plans a plusieurs facteurs, D31).
# P(difference) a cotes egales ; 1000 tables par scenario.
#   t22_H0 / t22_eff : 2 groupes de 20, p = 0,4 ; puis 0,3 contre 0,6
#   t32_H0 / t32_eff : 3 groupes de 15, p = 0,5 ; puis 0,5 / 0,5 / 0,8
# Mesures : ecart |P1 - P2| (mediane, 95e centile, maximum), part des tables ou les deux tombent
# dans deux classes differentes (non < .25, non tranche, oui > .75), BF > 3 et BF > 20 de chacun.
args <- commandArgs(TRUE)
sc <- args[1L]
N <- if (length(args) >= 2L) as.integer(args[2L]) else 1000L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
r <- t(vapply(seq_len(N), function(i) {
	set.seed(820000L + i + 10000L * match(sc, c("t22_H0", "t22_eff", "t32_H0", "t32_eff")))
	k <- if (startsWith(sc, "t22")) 2L else 3L; n <- if (k == 2L) 20L else 15L
	p <- switch(sc, t22_H0 = c(.4, .4), t22_eff = c(.3, .6), t32_H0 = rep(.5, 3), t32_eff = c(.5, .5, .8))
	g <- factor(rep(seq_len(k), each = n)); y <- stats::rbinom(k * n, 1, p[g])
	tab <- table(g, factor(y, levels = 0:1))
	l1 <- .bay_cont_lbf(tab, 2)
	z <- stats::glm(y ~ 1, family = stats::binomial())$deviance - stats::glm(y ~ g, family = stats::binomial())$deviance
	l2 <- .bay_dev_lbf(z, k - 1, k * n)
	p1 <- 1 / (1 + exp(-l1)); p2 <- 1 / (1 + exp(-l2))
	cl <- function(p) findInterval(p, c(0.25, 0.75))
	c(abs(p1 - p2), cl(p1) != cl(p2), l1 > log(3), l1 > log(20), l2 > log(3), l2 > log(20))
}, numeric(6)))
out <- c(R.version.string, sprintf("agent_82 : oui/non a un facteur, scenario %s, N = %d", sc, N),
	sprintf("|P marges neutres - P deviance| : mediane %.3f, 95e centile %.3f, max %.3f | classes differentes : %.3f",
		stats::median(r[, 1]), stats::quantile(r[, 1], 0.95), max(r[, 1]), mean(r[, 2])),
	sprintf("marges neutres BF>3 %.3f, BF>20 %.3f | deviance BF>3 %.3f, BF>20 %.3f", mean(r[, 3]), mean(r[, 4]), mean(r[, 5]), mean(r[, 6])))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_82_%s.txt", sc))); cat(out, sep = "\n")
