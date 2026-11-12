# G12 : ce que dit le bilan avec a = 2 et sa grille a = 1, 2, 4, ou (GBAY_A=3) a = 3 et sa grille
# a = 2, 3, 4 (chances de depart
# 1:2 a 2:1), garde-fou compris, sur les cas les plus durs :
#   sous H0, 2 x 6 a marges tres inegales (N = 60, 120) : « lien tres probable » (P >= 0,95, ordre stable) ;
#   sous un lien modere (V = 0,2, N = 100 et 200 ; 3 x 3, 3 x 4, 4 x 5) : « independance forte »
#   (P(indep) >= 0,95, ordre stable, aucune case signalee).
#   Rscript dev/gbay/scripts/sim_a2_bilan.R ; GBAY_A=3 Rscript dev/gbay/scripts/sim_a2_bilan.R
set.seed(20261115)
R <- as.integer(Sys.getenv("GBAY_REPS", "1000"))
S <- 1000L
ldm <- function(n, a) lgamma(sum(a)) - lgamma(sum(a) + sum(n)) + sum(lgamma(a + n) - lgamma(a))
lbf <- function(t, a) sum(apply(t, 1L, ldm, a = rep(a, ncol(t)))) - ldm(colSums(t), rep(nrow(t) * a, ncol(t)))
pr <- function(l, p = 0.5) 1 / (1 + exp(-l) * (1 - p) / p)
flagged <- function(t) {
	r <- nrow(t); cc <- ncol(t); k <- r * cc
	g <- matrix(stats::rgamma(S * k, shape = rep(as.vector(t) + 1, each = S)), S, k)
	P <- g / rowSums(g)
	prr <- vapply(seq_len(r), function(i) rowSums(P[, (seq_len(cc) - 1) * r + i, drop = FALSE]), numeric(S))
	pc <- vapply(seq_len(cc), function(j) rowSums(P[, (j - 1) * r + seq_len(r), drop = FALSE]), numeric(S))
	pm <- colMeans(P / (prr[, rep(seq_len(r), cc)] * pc[, rep(seq_len(cc), each = r)]) > 1)
	any(pm >= 0.975 | pm <= 0.025)
}
# Verdict du bilan : +2 lien tres probable, +1 probable, -2 independance forte (apres garde-fou), 0 sinon.
A0 <- as.numeric(Sys.getenv("GBAY_A", "2"))
GRID <- if (A0 == 2) c(1, 2, 4) else c(2, 3, 4)
said <- function(t, a0 = A0, grid = GRID) {
	lb <- lbf(t, a0); p <- pr(lb)
	g <- outer(vapply(grid, function(a) lbf(t, a), numeric(1)), c(1 / 3, 1 / 2, 2 / 3), function(l, q) pr(l, q))
	lead <- p >= 0.5
	open <- if (lead) min(g) <= 0.5 else max(g) >= 0.5
	if (open || abs(lb) < log(3)) return(0L)
	if (lead) return(if (p >= 0.95) 2L else if (p >= 0.75) 1L else 0L)
	if (1 - p >= 0.95 && !flagged(t)) return(-2L)
	0L
}
cv <- function(P) { E <- outer(rowSums(P), colSums(P)); sqrt(sum((P - E)^2 / E) / (min(dim(P)) - 1)) }
mk <- function(r, cc, V) {
	P0 <- outer(r, cc); D <- matrix(-1 / (length(cc) - 1), length(r), length(cc))
	for (i in seq_len(min(length(r), length(cc)))) D[i, i] <- 1
	d <- stats::uniroot(function(d) { P <- P0 * (1 + d * D); if (any(P <= 0)) return(1e6); cv(P / sum(P)) - V }, c(0, 0.999 * (length(cc) - 1)))$root
	P <- P0 * (1 + d * D); P / sum(P)
}
run <- function(P, N) {
	vapply(seq_len(R), function(i) {
		t <- matrix(stats::rmultinom(1, N, as.numeric(P)), nrow(P)); t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]
		if (nrow(t) < 2L || ncol(t) < 2L) return(NA_integer_)
		said(t)
	}, integer(1))
}
out <- NULL
Pu <- outer(c(.8, .2), c(.4, .3, .1, .1, .05, .05))
for (N in c(60, 120)) { s <- run(Pu, N); s <- s[!is.na(s)]
	out <- rbind(out, data.frame(scenario = paste0("H0 2x6 0,8/0,2, N = ", N), yes_very = mean(s == 2), yes_any = mean(s >= 1), indep_strong = mean(s == -2))) }
for (d in list(c(3, 3), c(3, 4), c(4, 5))) for (N in c(100, 200)) {
	s <- run(mk(rep(1 / d[1], d[1]), rep(1 / d[2], d[2]), 0.2), N); s <- s[!is.na(s)]
	out <- rbind(out, data.frame(scenario = paste0("H1 V = 0,2 ", d[1], "x", d[2], ", N = ", N), yes_very = mean(s == 2), yes_any = mean(s >= 1), indep_strong = mean(s == -2)))
}
print(out, digits = 3)
saveRDS(out, paste0("dev/gbay/scripts/sim_a", A0, "_bilan.rds"))
