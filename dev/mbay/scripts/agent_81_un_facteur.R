# agent_81 : voies a un facteur des notes entieres et des comptages (I1 de l'expert, passe 3),
# BIC (ancien) contre facteur de Bayes de deviance a g de Zellner-Siow, r = sqrt(2)/4 (nouveau),
# sur les memes jeux. 2 groupes de 20, 1000 jeux par scenario, sans sensibilite ni post-hoc.
#   notes_H0 / notes_06 : notes de 0 a 10, normale (6 ; 1,8) arrondie (le jeu de l'expert),
#                         ecart de 0,6 ecart-type pour le second groupe ; routage automatique
#   cnt_H0 / cnt_eff    : comptages, binomiale negative taille 2, moyenne 4, puis 4 -> 6,
#                         response = "counts"
# Mesures : « oui » P > .75, « non » P < .25 (a tort sous effet), P > .95 ; voie choisie ;
# memes jeux en voie continue (response = "continuous", loi normale et lois du bruit) pour les
# notes ; test t ou rapport de vraisemblance a p < .05.
args <- commandArgs(TRUE)
sc <- args[1L]
N <- if (length(args) >= 2L) as.integer(args[2L]) else 1000L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
r <- t(vapply(seq_len(N), function(i) {
	set.seed(810000L + i + 10000L * match(sc, c("notes_H0", "notes_06", "cnt_H0", "cnt_eff")))
	g <- factor(rep(c("a", "b"), each = 20))
	q <- function(...) suppressMessages(suppressWarnings(m.bay(..., verbose = FALSE, sensitivity = FALSE, posthoc = FALSE)))
	if (startsWith(sc, "notes")) {
		y <- pmin(10, pmax(0, round(stats::rnorm(40, 6 + (if (sc == "notes_06") 0.6 * 1.8 else 0) * (g == "b"), 1.8))))
		rr <- q(y ~ g)
		yo <- factor(y, levels = sort(unique(y)), ordered = TRUE)
		f0 <- .bay_polr_lm(yo, factor(rep(1L, 40))); f1 <- .bay_polr_lm(yo, g)
		old <- f1$lm - f0$lm
		cont <- q(y ~ g, response = "continuous", draws = 1000)$prob_effect
		pv <- stats::t.test(y ~ g)$p.value
		route <- if (identical(rr$plan, "ordinal")) 1 else 0
	} else {
		y <- stats::rnbinom(40, mu = ifelse(g == "b" & sc == "cnt_eff", 6, 4), size = 2)
		rr <- q(y ~ g, response = "counts")
		f0 <- .bay_counts_fit(y, factor(rep(1L, 40))); f1 <- .bay_counts_fit(y, g)
		lm <- c(f0$lm_pois, f1$lm_pois, f0$lm_nb, f1$lm_nb); ok <- is.finite(lm)
		w <- exp(lm[ok] - max(lm[ok])); nm <- c("p0", "p1", "n0", "n1")[ok]
		p1 <- sum(w[nm %in% c("p1", "n1")]) / sum(w)
		old <- log(p1) - log1p(-p1)
		cont <- NA_real_
		g1 <- MASS::glm.nb(y ~ g); g0 <- MASS::glm.nb(y ~ 1)
		pv <- stats::pchisq(2 * (as.numeric(stats::logLik(g1)) - as.numeric(stats::logLik(g0))), 1, lower.tail = FALSE)
		route <- 1
	}
	p <- rr$prob_effect; po <- 1 / (1 + exp(-old))
	c(p > 0.75, p < 0.25, p > 0.95, po > 0.75, po < 0.25, po > 0.95, cont > 0.75, cont < 0.25, pv < 0.05, route)
}, numeric(10)))
wil <- function(k, n) { z <- 1.96; ph <- k / n
	c((ph + z^2 / (2 * n) - z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n), (ph + z^2 / (2 * n) + z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)) }
f <- function(j) { v <- r[, j]; v <- v[!is.na(v)]; if (!length(v)) return("-"); k <- sum(v); n <- length(v); w <- wil(k, n); sprintf("%.3f [%.3f ; %.3f]", k / n, w[1], w[2]) }
out <- c(R.version.string, sprintf("agent_81 : un facteur, scenario %s, N = %d (voie ordinale ou comptages choisie : %.3f)", sc, N, mean(r[, 10])),
	paste0("nouveau (deviance ZS) : oui ", f(1), " | non ", f(2), " | P>.95 ", f(3)),
	paste0("ancien (BIC)          : oui ", f(4), " | non ", f(5), " | P>.95 ", f(6)),
	paste0("voie continue         : oui ", f(7), " | non ", f(8)),
	paste0("test t ou rapport de vraisemblance p<.05 : ", f(9)))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_81_%s.txt", sc))); cat(out, sep = "\n")
