# agent_70 : le defaut des queues lourdes a plusieurs facteurs, mesure avant toute correction.
# 2 x 2, n = 12 par case, bruit de Student a 3 ddl (ecart-type 1). Scenarios : H0, et effet de
# A de 0,8. Part des jeux ou P(A) > .75 / .95 pour m.bay a plusieurs facteurs (normale seule,
# D19 compris), p < .05 pour aov et pour aov sur les rangs, et P(A) de la voie a un facteur
# (4 lois, t comprises) appliquee a A seul, B ignore. 3e argument : suffixe du fichier de
# resultat (« apres » : rejoue apres D24, loi t a plusieurs facteurs).
args <- commandArgs(TRUE)
sc <- if (length(args) >= 1L) args[1L] else "H0"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 150L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
r <- t(vapply(seq_len(N), function(i) {
	set.seed(700000L + i + if (sc == "H0") 0L else 9000L)
	d <- expand.grid(rep = 1:12, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
	d$y <- stats::rt(nrow(d), 3) / sqrt(3) + if (sc == "H0") 0 else 0.8 * (d$A == "a2")
	f <- suppressMessages(m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE, draws = 2000))
	pA <- f$effects$p_h1[1L]
	o <- suppressMessages(m.bay(y ~ A, data = d, verbose = FALSE, sensitivity = FALSE, draws = 2000, posthoc = FALSE))
	pa <- summary(stats::aov(y ~ A * B, data = d))[[1]][1L, "Pr(>F)"]
	d$rk <- rank(d$y)
	pr <- summary(stats::aov(rk ~ A * B, data = d))[[1]][1L, "Pr(>F)"]
	c(f75 = pA > 0.75, f95 = pA > 0.95, o75 = o$prob_effect > 0.75, o95 = o$prob_effect > 0.95,
		heavy = isTRUE(o$heavy), aov = pa < 0.05, rang = pr < 0.05)
}, numeric(7)))
out <- c(R.version.string, sprintf("agent_70 : 2 x 2, bruit t3, scenario %s, N = %d", sc, N),
	sprintf("m.bay facteurs (normale) P(A)>.75 : %.3f | >.95 : %.3f || voie un facteur (4 lois) sur A P>.75 : %.3f | >.95 : %.3f | queues lourdes dites : %.3f || aov p<.05 : %.3f | aov rangs p<.05 : %.3f",
		mean(r[, 1]), mean(r[, 2]), mean(r[, 3]), mean(r[, 4]), mean(r[, 5]), mean(r[, 6]), mean(r[, 7])))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_70_%s%s.txt", sc, if (length(args) >= 3L) paste0("_", args[3L]) else ""))); cat(out, sep = "\n")
