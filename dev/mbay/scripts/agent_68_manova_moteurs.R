# agent_68 : choix du moteur de la MANOVA bayesienne (D22). Trois candidats, face a Pillai :
#   A. BFpack::BF sur aov multivarie (facteur de Bayes fractionnaire, Mulder et Gu) ;
#   B. BIC multivarie (loi normale a covariance libre, information unitaire, Kass et Raftery 1995) ;
#   C. lmBF sur les donnees empilees, reponses standardisees : y ~ g * vd + sujet aleatoire
#      (« au moins un ecart » = g et g:vd, moyenne sur les modeles), BayesFactor.
# Scenarios : H0 (d = 3 et d = 10), effet sur une VD (d = 3, d = 10), effet reparti (d = 3),
# et correlations tres inegales (0,8 / 0 / 0) sous H0 et avec effet sur une VD.
# 3 groupes de n = 15. N jeux par scenario. Mesures : P(effet) > .75 et > .95, Pillai p < .05.
args <- commandArgs(TRUE)
sc <- if (length(args) >= 1L) args[1L] else "H0_d3"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 100L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages({library(BFpack); library(BayesFactor)})
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
cfg <- list(H0_d3 = list(d = 3, eff = c(0, 0, 0)), H0_d10 = list(d = 10, eff = rep(0, 10)),
	une_d3 = list(d = 3, eff = c(0.8, 0, 0)), une_d10 = list(d = 10, eff = c(0.8, rep(0, 9))),
	rep_d3 = list(d = 3, eff = c(0.5, 0.5, 0.5)),
	H0_cor = list(d = 3, eff = c(0, 0, 0), R = matrix(c(1, 0.8, 0, 0.8, 1, 0, 0, 0, 1), 3)),
	une_cor = list(d = 3, eff = c(0.8, 0, 0), R = matrix(c(1, 0.8, 0, 0.8, 1, 0, 0, 0, 1), 3)))[[sc]]
n <- 15L
res <- t(vapply(seq_len(N), function(i) {
	set.seed(680000L + i)
	g <- factor(rep(c("a", "b", "c"), each = n))
	R <- matrix(0.3, cfg$d, cfg$d); diag(R) <- 1
	if (!is.null(cfg$R)) R <- cfg$R  # correlations far from equal: the stacked model assumes them equal
	Y <- matrix(stats::rnorm(3 * n * cfg$d), ncol = cfg$d) %*% chol(R)
	Y[g == "c", ] <- sweep(Y[g == "c", , drop = FALSE], 2L, cfg$eff, "+")
	colnames(Y) <- paste0("y", seq_len(cfg$d))
	pil <- summary(stats::manova(Y ~ g), test = "Pillai")$stats[1L, "Pr(>F)"]
	# A
	pa <- tryCatch({ b <- BFpack::BF(stats::aov(Y ~ g)); b$PHP_main[1L, 2L] }, error = function(e) NA_real_)
	# B
	f1 <- stats::lm(Y ~ g); f0 <- stats::lm(Y ~ 1)
	ll <- function(f) { E <- stats::residuals(f); S <- crossprod(E) / nrow(E); -nrow(E) / 2 * (as.numeric(determinant(S)$modulus) + ncol(E) * (1 + log(2 * pi))) }
	lbf <- (ll(f1) - ll(f0)) - 0.5 * (2 * cfg$d) * log(3 * n)
	pb <- 1 / (1 + exp(-lbf))
	# C
	Z <- scale(Y)
	long <- data.frame(y = as.vector(Z), g = rep(g, cfg$d), vd = factor(rep(colnames(Y), each = 3 * n)), id = factor(rep(seq_len(3 * n), cfg$d)))
	bf <- tryCatch(BayesFactor::generalTestBF(y ~ g * vd + id, data = long, whichRandom = "id", neverExclude = c("id", "^vd$"), progress = FALSE), error = function(e) NULL)
	pc <- NA_real_
	if (!is.null(bf)) {
		lb <- BayesFactor::extractBF(bf, logbf = TRUE)$bf; nm <- rownames(BayesFactor::extractBF(bf))
		base <- lb[nm == "vd + id"]
		has_g <- grepl("g", gsub("vd", "", gsub("id", "", nm)))
		w <- exp(c(base, lb[has_g]) - max(c(base, lb[has_g])))
		pc <- 1 - w[1L] / sum(w) # equal prior on the models with g, against the base
		pc <- 1 - 1 / (1 + sum(exp(lb[has_g] - base)) / sum(has_g))
	}
	c(pil = pil < 0.05, a75 = pa > 0.75, a95 = pa > 0.95, b75 = pb > 0.75, b95 = pb > 0.95, c75 = pc > 0.75, c95 = pc > 0.95)
}, numeric(7)))
out <- c(R.version.string, sprintf("agent_68 : MANOVA, scenario %s, 3 groupes x %d, N = %d", sc, n, N),
	sprintf("Pillai p<.05 : %.3f || A BFpack P>.75 %.3f P>.95 %.3f || B BIC P>.75 %.3f P>.95 %.3f || C lmBF empile P>.75 %.3f P>.95 %.3f",
		mean(res[, 1]), mean(res[, 2], na.rm = TRUE), mean(res[, 3], na.rm = TRUE), mean(res[, 4]), mean(res[, 5]), mean(res[, 6], na.rm = TRUE), mean(res[, 7], na.rm = TRUE)))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_68_%s.txt", sc))); cat(out, sep = "\n")
