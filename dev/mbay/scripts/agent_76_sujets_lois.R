# agent_76 : queues lourdes et variances inegales avec des sujets (mesures repetees, D15),
# mesurees avant toute correction. 15 sujets x 3 conditions ; effet sujet d'ecart-type 1.
# Scenarios :
#   t3_H0 / t3_A   : residus de Student a 3 ddl (ecart-type 1) ; effet de 0,6 sur la condition 3.
#   het_H0 / het_A : residus normaux d'ecart-types (0,5 ; 0,5 ; 2) selon la condition ; effet 0,6
#                    sur la condition 1 (la plus precise).
# Part des jeux ou P(cond) > .75 / .95 (m.bay, sujet aleatoire, loi normale), p < .05 pour
# l'ANOVA a mesures repetees (non corrigee) et pour Friedman. Temoins norm_H0 / norm_A (residus
# normaux, variances egales). 3e argument : suffixe (« D28 » : rejoue avec la prudence provisoire).
args <- commandArgs(TRUE)
sc <- if (length(args) >= 1L) args[1L] else "t3_H0"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 200L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
n <- 15L; k <- 3L
r <- t(vapply(seq_len(N), function(i) {
	set.seed(760000L + i + match(sc, c("t3_H0", "t3_A", "het_H0", "het_A", "norm_H0", "norm_A")) * 10000L)
	u <- stats::rnorm(n)
	E <- if (startsWith(sc, "t3")) matrix(stats::rt(n * k, 3) / sqrt(3), n, k) else if (startsWith(sc, "norm")) matrix(stats::rnorm(n * k), n, k) else sweep(matrix(stats::rnorm(n * k), n, k), 2L, c(0.5, 0.5, 2), "*")
	mu <- switch(sc, t3_A = c(0, 0, 0.6), het_A = c(0.6, 0, 0), norm_A = c(0, 0, 0.6), c(0, 0, 0))
	Y <- u + E + rep(mu, each = n)
	d <- data.frame(y = as.vector(Y), cond = factor(rep(paste0("c", 1:k), each = n)), id = factor(rep(seq_len(n), k)))
	rr <- suppressMessages(m.bay(y ~ cond | id, data = d, verbose = FALSE, sensitivity = FALSE, posthoc = FALSE, draws = 2000))
	p <- rr$effects$p_h1[1L]
	pa <- summary(stats::aov(y ~ cond + Error(id / cond), data = d))[["Error: id:cond"]][[1]][1L, "Pr(>F)"]
	pf <- stats::friedman.test(Y)$p.value
	cap <- isTRUE(rr$het_cells)
	c(p75 = p > 0.75, p95 = p > 0.95, aov = pa < 0.05, fried = pf < 0.05, cap = cap, shown95 = p > 0.95 && !cap)
}, numeric(6)))
out <- c(R.version.string, sprintf("agent_76 : sujets, 15 x 3, scenario %s, N = %d", sc, N),
	sprintf("m.bay P(cond)>.75 : %.3f | >.95 : %.3f || ANOVA MR p<.05 : %.3f | Friedman p<.05 : %.3f || prudence D28 declenchee : %.3f | oui fort affiche (P>.95 sans plafond) : %.3f",
		mean(r[, 1]), mean(r[, 2]), mean(r[, 3]), mean(r[, 4]), mean(r[, 5]), mean(r[, 6])))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_76_%s%s.txt", sc, if (length(args) >= 3L) paste0("_", args[3L]) else ""))); cat(out, sep = "\n")
