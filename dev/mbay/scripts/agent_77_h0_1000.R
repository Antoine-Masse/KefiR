# agent_77 : recalibration sous H0 sur au moins 1000 jeux, avec intervalle de Wilson a 95 %
# (audit de l'expert du 11 novembre : 60 a 200 jeux ne permettent pas de comparer a 0,7 %).
# Memes plans que les calibrations d'origine :
#   d23  : MANOVA 2 x 2, 3 reponses (correlations 0,3), n = 12 par case (agent_69)
#   d25b : oui/non 2 x 2, p = 0,4, n = 25 par case (agent_73)
#   d25c : comptages 2 x 2, binomiale negative taille 2, moyenne 4, response = "counts" (agent_73)
#   d26  : notes 1 a 5, 2 x 2, n = 25 par case, response = "ordinal" (agent_75)
#   d27  : mesures repetees, 15 sujets x 4 temps, AR(1) rho = 0,8 (agent_74)
# Mesure : part des jeux ou au moins un terme est a P > .95 (et > .75), sans sensibilite.
# Arguments : scenario, N, decalage des graines (tranches paralleles).
args <- commandArgs(TRUE)
sc <- args[1L]
N <- as.integer(args[2L])
off <- if (length(args) >= 3L) as.integer(args[3L]) else 0L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
one <- function(i) {
	set.seed(770000L + 100000L * match(sc, c("d23", "d25b", "d25c", "d26", "d27")) + off + i)
	q <- function(...) suppressMessages(suppressWarnings(m.bay(..., verbose = FALSE, sensitivity = FALSE)))
	p <- switch(sc,
		d23 = {
			d <- expand.grid(rep = 1:12, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
			R <- matrix(0.3, 3, 3); diag(R) <- 1
			Y <- matrix(stats::rnorm(nrow(d) * 3), ncol = 3) %*% chol(R)
			d$y1 <- Y[, 1]; d$y2 <- Y[, 2]; d$y3 <- Y[, 3]
			q(cbind(y1, y2, y3) ~ A * B, data = d, draws = 2000)$effects$p_h1
		},
		d25b = {
			d <- expand.grid(rep = 1:25, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
			d$y <- factor(ifelse(stats::rbinom(nrow(d), 1, 0.4) == 1, "oui", "non"), levels = c("non", "oui"))
			q(y ~ A * B, data = d)$effects$p_h1
		},
		d25c = {
			d <- expand.grid(rep = 1:25, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
			d$y <- stats::rnbinom(nrow(d), mu = 4, size = 2)
			q(y ~ A * B, data = d, response = "counts")$effects$p_h1
		},
		d26 = {
			d <- expand.grid(rep = 1:25, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
			d$y <- pmin(5, pmax(1, round(3 + stats::rnorm(nrow(d)))))
			q(y ~ A * B, data = d, response = "ordinal")$effects$p_h1
		},
		d27 = {
			k <- 4L; n <- 15L
			Y <- matrix(stats::rnorm(n * k), n, k) %*% chol(0.8^abs(outer(1:k, 1:k, "-")))
			d <- data.frame(y = as.vector(Y), t = factor(rep(paste0("T", 1:k), each = n)), id = factor(rep(1:n, k)))
			q(y ~ t | id, data = d, posthoc = FALSE, draws = 2000)$effects$p_h1
		})
	c(any(p > 0.95), any(p > 0.75))
}
r <- t(vapply(seq_len(N), one, numeric(2)))
wil <- function(k, n) {
	z <- 1.96; ph <- k / n
	c(ph, (ph + z^2 / (2 * n) - z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
		(ph + z^2 / (2 * n) + z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n))
}
k95 <- sum(r[, 1]); k75 <- sum(r[, 2])
w95 <- wil(k95, N); w75 <- wil(k75, N)
out <- c(R.version.string, sprintf("agent_77 : H0, scenario %s, N = %d, graines +%d", sc, N, off),
	sprintf("un terme P>.95 : %d/%d = %.4f [Wilson 95 %% : %.4f ; %.4f] | P>.75 : %d/%d = %.4f [%.4f ; %.4f]",
		k95, N, w95[1], w95[2], w95[3], k75, N, w75[1], w75[2], w75[3]))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_77_%s%s.txt", sc, if (off) paste0("_", off) else ""))); cat(out, sep = "\n")
