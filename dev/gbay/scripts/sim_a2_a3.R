# G12 : a = 2 ou a = 3 ? (desaccord avec l'expert, passe 2). Sous H0, P(lien) > 0,95 (BF > 19) et
# BF > 20 pour a = 2 et 3 ; scenarios demandes : 2 x 6, N = 60, marges tres inegales (0,8 / 0,2 ;
# 0,4 / 0,3 / 0,1 / 0,1 / 0,05 / 0,05, ceux de l'expert) et 4 x 5, N = 20 ; plus 2 x 6 a N = 30 et
# 120 et 4 x 5, N = 20 a marges inegales. Puis mediane a priori du V de Cramer sous H1.
#   Rscript dev/gbay/scripts/sim_a2_a3.R      (un seul processus, sans paquet)
set.seed(20261114)
R <- as.integer(Sys.getenv("GBAY_REPS", "1000"))
ldm <- function(n, a) lgamma(sum(a)) - lgamma(sum(a) + sum(n)) + sum(lgamma(a + n) - lgamma(a))
lbf <- function(t, a) sum(apply(t, 1L, ldm, a = rep(a, ncol(t)))) - ldm(colSums(t), rep(nrow(t) * a, ncol(t)))
wilson <- function(k, n) { z <- 1.96; p <- k / n; c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
	(p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)) }
sc <- list(
	list(lab = "2x6, 0,8/0,2, N = 60", P = outer(c(.8, .2), c(.4, .3, .1, .1, .05, .05)), N = 60),
	list(lab = "2x6, 0,8/0,2, N = 30", P = outer(c(.8, .2), c(.4, .3, .1, .1, .05, .05)), N = 30),
	list(lab = "2x6, 0,8/0,2, N = 120", P = outer(c(.8, .2), c(.4, .3, .1, .1, .05, .05)), N = 120),
	list(lab = "4x5 egales, N = 20", P = matrix(1 / 20, 4, 5), N = 20),
	list(lab = "4x5 inegales, N = 20", P = outer((1:4) / 10, (1:5) / 15), N = 20))
out <- NULL
for (s in sc) {
	X <- stats::rmultinom(R, s$N, as.numeric(s$P))
	m <- vapply(seq_len(R), function(i) {
		t <- matrix(X[, i], nrow(s$P)); t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]
		if (nrow(t) < 2L || ncol(t) < 2L) return(c(NA, NA))
		c(lbf(t, 2), lbf(t, 3))
	}, numeric(2))
	m <- m[, !is.na(m[1L, ]), drop = FALSE]
	for (k in 1:2) {
		n95 <- sum(m[k, ] > log(19)); w <- wilson(n95, ncol(m))
		out <- rbind(out, data.frame(scenario = s$lab, a = k + 1, n_ok = ncol(m), p95 = n95 / ncol(m), wilson_low = w[1], wilson_high = w[2],
			bf20 = mean(m[k, ] > log(20)), p75 = mean(m[k, ] > log(3))))
	}
}
cat("== Sous H0 : P(lien) > 0,95 (BF > 19), intervalle de Wilson, BF > 20, P > 0,75 ==\n"); print(out, digits = 3, width = 200)
V <- function(P) { E <- outer(rowSums(P), colSums(P)); sqrt(sum((P - E)^2 / E) / (min(dim(P)) - 1)) }
pv <- NULL
for (d in list(c(2, 2), c(2, 3), c(3, 3), c(3, 4), c(4, 5), c(2, 6))) for (a in c(2, 3)) {
	v <- replicate(4000, { p <- stats::rgamma(prod(d), a); V(matrix(p / sum(p), d[1], d[2])) })
	q <- stats::quantile(v, c(.5, .1, .9))
	pv <- rbind(pv, data.frame(table = paste0(d[1], "x", d[2]), a = a, median = q[1], q10 = q[2], q90 = q[3]))
}
rownames(pv) <- NULL
cat("\n== V de Cramer a priori sous H1 (cases Dirichlet(a)) ==\n"); print(pv, digits = 2)
saveRDS(list(h0 = out, prior_V = pv, R = R), "dev/gbay/scripts/sim_a2_a3.rds")
