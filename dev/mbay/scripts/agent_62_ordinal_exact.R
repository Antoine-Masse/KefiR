# agent_62 : approximation BIC de la voie ordinale contre une approximation de Laplace sous
# l'a priori propre d'information unitaire (Kass et Wasserman 1995), que le BIC approche.
# Effets beta ~ N(0, n V) (V : covariance de l'estimateur de polr), seuils communs aux deux
# modeles. brms (family = cumulative) + bridgesampling devait servir de reference exacte :
# la compilation Stan echoue sur ce poste (R 4.6 attend Rtools 4.5, seul Rtools 4.4 est
# installe ; journal dans resultats/agent_62.txt).
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
cuts <- stats::qnorm(1:4 / 5)
cases <- list(list(d = c(0, 0, 0), n = 20), list(d = c(0, 0, 0.6), n = 20), list(d = c(0, 0.5, 1), n = 15),
	list(d = c(0, 0, 0.4), n = 40))
out <- c("agent_62 : brms + bridgesampling impossible sur ce poste (compilation Stan : R 4.6 attend Rtools 4.5, absent).",
	"Reference de remplacement : Laplace sous a priori d'information unitaire propre, beta ~ N(0, n V).")
for (i in seq_along(cases)) {
	cs <- cases[[i]]
	set.seed(620 + i)
	g <- factor(rep(c("a", "b", "c"), each = cs$n))
	yo <- factor(findInterval(cs$d[as.integer(g)] + stats::rnorm(length(g)), cuts) + 1L, ordered = TRUE)
	m1 <- MASS::polr(yo ~ g, Hess = TRUE); m0 <- MASS::polr(yo ~ 1, Hess = TRUE)
	bic <- (-BIC(m1) + BIC(m0)) / 2
	b <- stats::coef(m1); V <- stats::vcov(m1)[names(b), names(b)]
	n <- length(yo); p <- length(b)
	S0 <- n * V                                   # unit-information prior covariance
	post_cov <- solve(solve(V) + solve(S0))
	lp <- -0.5 * (p * log(2 * pi) + as.numeric(determinant(S0)$modulus) + drop(t(b) %*% solve(S0) %*% b))
	lap <- as.numeric(logLik(m1)) + lp + 0.5 * (p * log(2 * pi) + as.numeric(determinant(post_cov)$modulus)) - as.numeric(logLik(m0))
	out <- c(out, sprintf("cas %d (d = %s, n = %d par groupe) : log BF  BIC %.2f | Laplace information unitaire %.2f | ecart %.2f",
		i, paste(cs$d, collapse = ","), cs$n, bic, lap, bic - lap))
}
writeLines(out, file.path(BASE, "dev/mbay/scripts/resultats/agent_62.txt"))
cat(out, sep = "\n")
