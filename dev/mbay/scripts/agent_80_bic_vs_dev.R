# agent_80 : moteur du discret a plusieurs facteurs (D31). Sur les memes jeux, -BIC/2 (D25)
# contre le facteur de Bayes de deviance a g de Zellner-Siow, r = sqrt(2)/4 (Held et al. 2015),
# avec l'a priori HOP sur les modeles. 2 x 2, n = 25 par case, 1000 jeux par scenario.
#   bin_H0 / bin_or2 : oui/non, p = 0,4 ; rapport de cotes 2 sur A
#   cnt_H0 / cnt_rr15: comptages binomiale negative taille 2, moyenne 4 ; rapport 1,5 sur A
# Mesures : un terme a P > .95 et > .75 (sous H0) ; sur A sous effet : « oui » P(A) > .75,
# « non » BF(A) <= 1/3, « non, soutien fort » BF(A) <= 1/20 ; rapport de vraisemblance p < .05.
args <- commandArgs(TRUE)
sc <- args[1L]
N <- if (length(args) >= 2L) as.integer(args[2L]) else 1000L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
sp <- .bay_fact_space(c("A", "B", "A:B"))
lp <- .bay_fact_prior(sp, "hop")
ans <- function(lml, kind) {
	lm_m <- if (kind == "binary") lml[, 1L] else .bay_lse2(lml[, "poisson"], lml[, "negbin"]) - log(2)
	lb <- lm_m - lm_m[1L]
	e <- .bay_fact_terms(sp, lb, lp)
	e[match(.bay_term_key(c("A", "B", "A:B")), .bay_term_key(e$term)), c("p_h1", "log_bf")]
}
r <- t(vapply(seq_len(N), function(i) {
	set.seed(800000L + i + 10000L * match(sc, c("bin_H0", "bin_or2", "cnt_H0", "cnt_rr15")))
	d <- expand.grid(rep = 1:25, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
	kind <- if (startsWith(sc, "bin")) "binary" else "counts"
	if (kind == "binary") {
		p <- if (sc == "bin_or2") ifelse(d$A == "a2", 0.8 / 1.4, 0.4) else 0.4
		y <- stats::rbinom(nrow(d), 1, p)
	} else {
		y <- stats::rnbinom(nrow(d), mu = if (sc == "cnt_rr15") ifelse(d$A == "a2", 6, 4) else 4, size = 2)
	}
	fits <- suppressWarnings(.bay_df_fits(y, d[c("A", "B")], sp, kind))
	a1 <- ans(.bay_df_lml_engine(fits, "bic"), kind)
	a2 <- ans(.bay_df_lml_engine(fits, "dev"), kind)
	d$y <- y
	g1 <- if (kind == "binary") stats::glm(y ~ A + B, data = d, family = stats::binomial()) else MASS::glm.nb(y ~ A + B, data = d)
	g0 <- if (kind == "binary") stats::glm(y ~ B, data = d, family = stats::binomial()) else MASS::glm.nb(y ~ B, data = d)
	lrt <- stats::pchisq(2 * (as.numeric(stats::logLik(g1)) - as.numeric(stats::logLik(g0))), 1, lower.tail = FALSE)
	st <- function(a) c(any(a$p_h1 > 0.95), any(a$p_h1 > 0.75), a$p_h1[1L] > 0.75, a$log_bf[1L] <= -log(3), a$log_bf[1L] <= -log(20))
	c(st(a1), st(a2), lrt < 0.05)
}, numeric(11)))
wil <- function(k, n) { z <- 1.96; ph <- k / n
	c((ph + z^2 / (2 * n) - z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n), (ph + z^2 / (2 * n) + z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)) }
f <- function(j) { k <- sum(r[, j], na.rm = TRUE); w <- wil(k, N); sprintf("%.3f [%.3f ; %.3f]", k / N, w[1], w[2]) }
lab <- c("un terme P>.95", "un terme P>.75", "A oui (P>.75)", "A non (BF<=1/3)", "A non fort (BF<=1/20)")
out <- c(R.version.string, sprintf("agent_80 : BIC contre deviance ZS sqrt(2)/4, scenario %s, N = %d", sc, N),
	paste0("BIC  : ", paste0(lab, " ", vapply(1:5, f, ""), collapse = " | ")),
	paste0("ZS   : ", paste0(lab, " ", vapply(6:10, f, ""), collapse = " | ")),
	paste0("rapport de vraisemblance p<.05 : ", f(11)))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_80_%s.txt", sc))); cat(out, sep = "\n")
