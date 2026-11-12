# G12 : forme close a marges neutres, concentration a = 1, 2, 3, et effet du garde-fou de coherence
# (cases dont l'intervalle a 95 % exclut 1) sur ce que dit le bilan.
#   Rscript dev/gbay/scripts/sim_neutral.R
# « Fort pour l'independance » = P(independance) >= 0,95 (soutien fort ou tres fort) ; avec le
# garde-fou, il ne compte que si aucune case n'est signalee.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
set.seed(20261113)
R <- as.integer(Sys.getenv("GBAY_REPS", "1000"))
S <- 1000L
t0 <- Sys.time()
ldm <- function(n, a) lgamma(sum(a)) - lgamma(sum(a) + sum(n)) + sum(lgamma(a + n) - lgamma(a))
lb_neutral <- function(t, a = 1) {
	r <- nrow(t); cc <- ncol(t)
	ldm(as.vector(t), rep(a, r * cc)) - ldm(rowSums(t), rep(cc * a, r)) - ldm(colSums(t), rep(r * a, cc))
}
# Cells flagged by g.bay's rule (G6): Dirichlet(n + 1) posterior, P(theta > 1) >= 0.975 or <= 0.025.
flagged <- function(t) {
	r <- nrow(t); cc <- ncol(t); k <- r * cc
	g <- matrix(stats::rgamma(S * k, shape = rep(as.vector(t) + 1, each = S)), S, k)
	P <- g / rowSums(g)
	pr <- vapply(seq_len(r), function(i) rowSums(P[, (seq_len(cc) - 1) * r + i, drop = FALSE]), numeric(S))
	pc <- vapply(seq_len(cc), function(j) rowSums(P[, (j - 1) * r + seq_len(r), drop = FALSE]), numeric(S))
	th <- P / (pr[, rep(seq_len(r), cc)] * pc[, rep(seq_len(cc), each = r)])
	pm <- colMeans(th > 1)
	any(pm >= 0.975 | pm <= 0.025)
}
pr <- function(l) 1 / (1 + exp(-l))
cv <- function(P) { E <- outer(rowSums(P), colSums(P)); sqrt(sum((P - E)^2 / E) / (min(dim(P)) - 1)) }
mk <- function(r, cc, V) {
	P0 <- outer(r, cc); D <- matrix(-1 / (length(cc) - 1), length(r), length(cc))
	for (i in seq_len(min(length(r), length(cc)))) D[i, i] <- 1
	f <- function(d) { P <- P0 * (1 + d * D); if (any(P <= 0)) return(1e6); cv(P / sum(P)) - V }
	d <- stats::uniroot(f, c(0, 0.999 * (length(cc) - 1)))$root
	P <- P0 * (1 + d * D); P / sum(P)
}
as_ <- c(1, 2, 3)
run <- function(P, N) {
	r <- nrow(P)
	m <- vapply(seq_len(R), function(i) {
		t <- matrix(stats::rmultinom(1, N, as.numeric(P)), r)
		t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]
		if (nrow(t) < 2L || ncol(t) < 2L) return(rep(NA_real_, 4))
		c(vapply(as_, function(a) lb_neutral(t, a), numeric(1)), flagged(t))
	}, numeric(4))
	m[, !is.na(m[1L, ]), drop = FALSE]
}
summ <- function(m, lab) {
	out <- lab
	for (k in seq_along(as_)) {
		p <- pr(m[k, ])
		out[[paste0("a", as_[k], "_p95")]] <- mean(p > 0.95)
		out[[paste0("a", as_[k], "_h0strong")]] <- mean(p < 0.05)
		out[[paste0("a", as_[k], "_h0strong_guard")]] <- mean(p < 0.05 & m[4L, ] == 0)
		out[[paste0("a", as_[k], "_p75")]] <- mean(p > 0.75)
	}
	out$flag <- mean(m[4L, ] == 1)
	out
}
h0 <- NULL
for (dm in list(c(2, 2), c(3, 3), c(3, 4), c(4, 5))) for (mg in c("egales", "inegales")) for (N in c(30, 60, 200)) {
	rr <- if (mg == "egales") rep(1, dm[1]) else seq(1, 3, length.out = dm[1])
	cc <- if (mg == "egales") rep(1, dm[2]) else seq(1, 4, length.out = dm[2])
	h0 <- rbind(h0, summ(run(outer(rr / sum(rr), cc / sum(cc)), N), data.frame(table = paste0(dm[1], "x", dm[2]), margins = mg, N = N)))
}
cat("== Sous H0 ==\n"); print(h0, digits = 3, width = 250)
h1 <- NULL
for (dm in list(c(2, 2), c(3, 3), c(3, 4), c(4, 5))) for (V in c(0.2, 0.3)) for (N in c(50, 100, 200)) {
	h1 <- rbind(h1, summ(run(mk(rep(1 / dm[1], dm[1]), rep(1 / dm[2], dm[2]), V), N), data.frame(table = paste0(dm[1], "x", dm[2]), V = V, N = N)))
}
cat("\n== Sous un lien ==\n"); print(h1, digits = 3, width = 250)
saveRDS(list(h0 = h0, h1 = h1, R = R), "dev/gbay/scripts/sim_neutral.rds")
cat("\nDuree :", format(Sys.time() - t0), "\n")
