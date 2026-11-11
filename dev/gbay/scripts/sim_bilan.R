# Taux de « oui » dit par le bilan de g.bay sous H0 (indépendance), sensibilité comprise :
# le bilan ne dit « un lien existe probablement » que si P >= 0,75, le facteur de Bayes >= 3 et
# l'ordre ne change pas quand a (1, 2, 4), les chances de départ (1:2 à 2:1) et, en mode
# automatique, le plan d'échantillonnage varient (DECISIONS, G3, G5, G10).
#   Rscript dev/gbay/scripts/sim_bilan.R
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(pkgload::load_all(".", quiet = TRUE))
set.seed(11112026)
R <- as.integer(Sys.getenv("GBAY_REPS", "2000"))
J <- function(t, a = 1) BayesFactor:::contingencyJointMultinomial(t, a)
Pn <- function(t, a = 1) BayesFactor:::contingencyPoisson(t, a)
Ir <- function(t, a = 1) BayesFactor:::contingencyIndepMultinomial(t, a)
Ic <- function(t, a = 1) BayesFactor:::contingencyIndepMultinomial(t(t), a)
for (i in 1:10) {
	tt <- matrix(stats::rpois(6, 5) + 1, 2)
	stopifnot(abs(Pn(tt) - .gbay_lbf_ct(tt, "poisson")) < 1e-10, abs(Ic(tt) - .gbay_lbf_ct(tt, "indepMulti", "cols")) < 1e-10)
}
pr <- function(l, p = 0.5) 1 / (1 + exp(-l) * (1 - p) / p)
phs <- c(1 / 3, 1 / 2, 2 / 3)

# Ce que dit le bilan : 0 (pas de oui), 1 (probablement), 2 (très probablement).
said <- function(t, plan = c("joint", "rows")) {
	plan <- match.arg(plan)
	f <- if (plan == "joint") J else Ir
	lbs <- c(f(t, 1), f(t, 2), f(t, 4))
	lb <- lbs[1L]; p <- pr(lb)
	if (p < 0.5) return(0L)
	pl <- c(J(t), Pn(t), Ir(t), Ic(t))
	vals <- c(outer(lbs, phs, pr), outer(pl, phs, pr))
	open <- min(vals) <= 0.5
	if (open || lb < log(3) || p < 0.75) return(0L)
	if (p >= 0.95) 2L else 1L
}
# Regle de Cochran (lecture courante) : plus de 20 % des cases d'attendu < 5, ou une case < 1.
sparse <- function(t) { E <- outer(rowSums(t), colSums(t)) / sum(t); mean(E < 5) > 0.2 || any(E < 1) }

out <- NULL
designs <- list(list(r = c(.5, .5), c = c(.5, .5)), list(r = c(.5, .5), c = c(.5, .3, .2)),
	list(r = c(.4, .35, .25), c = c(.1, .3, .3, .3)), list(r = rep(.25, 4), c = rep(.2, 5)))
for (ds in designs) {
	for (N in c(10, 20, 50, 100, 200)) {
		P <- outer(ds$r, ds$c)
		X <- stats::rmultinom(R, N, as.numeric(P))
		s <- vapply(seq_len(R), function(i) {
			t <- matrix(X[, i], length(ds$r)); t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]
			if (nrow(t) < 2L || ncol(t) < 2L) return(c(NA, NA, NA, NA))
			c(said(t, "joint"), pr(J(t)), mean(outer(rowSums(t), colSums(t)) / sum(t)), sparse(t))
		}, numeric(4))
		s <- s[, !is.na(s[1L, ]), drop = FALSE]
		nr <- length(ds$r); ng <- max(2L, round(N / nr))
		s2 <- vapply(seq_len(R), function(i) {
			t <- t(stats::rmultinom(nr, ng, ds$c)); t <- t[, colSums(t) > 0, drop = FALSE]
			if (ncol(t) < 2L) return(c(NA_real_, NA_real_))
			c(said(t, "rows"), sparse(t))
		}, numeric(2))
		s2 <- s2[, !is.na(s2[1L, ]), drop = FALSE]
		out <- rbind(out, data.frame(table = paste0(length(ds$r), "x", length(ds$c)), N = N, mean_expected = round(N / length(P), 1),
			joint_raw_p75 = mean(s[2L, ] > 0.75), joint_bilan_yes = mean(s[1L, ] >= 1), joint_bilan_very = mean(s[1L, ] == 2),
			rows_bilan_yes = mean(s2[1L, ] >= 1), rows_bilan_very = mean(s2[1L, ] == 2),
			sparse = mean(s[4L, ] == 1), joint_capped_yes = mean(s[1L, ] >= 1 & s[4L, ] == 0), rows_capped_yes = mean(s2[1L, ] >= 1 & s2[2L, ] == 0)))
	}
}
print(out, digits = 3, width = 200)

# Puissance du bilan sous H1, avec et sans la prudence des tables creuses.
pw <- NULL
for (cs in list(list(P = matrix(c(.25, .15, .1, .1, .15, .25), 2), lab = "2x3 fort"), list(P = matrix(c(.3, .2, .2, .3), 2), lab = "2x2 moyen"),
	list(P = outer(c(.4, .35, .25), c(.1, .3, .3, .3)) * (1 + 0.6 * outer(c(1, 0, -1), c(1, 1, -1, -1))), lab = "3x4 moyen"))) {
	P <- cs$P / sum(cs$P)
	for (N in c(20, 50, 100, 200)) {
		X <- stats::rmultinom(R / 2, N, as.numeric(P))
		s <- vapply(seq_len(R / 2), function(i) {
			t <- matrix(X[, i], nrow(P)); t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]
			if (nrow(t) < 2L || ncol(t) < 2L) return(c(NA, NA))
			c(said(t, "joint"), sparse(t))
		}, numeric(2))
		s <- s[, !is.na(s[1L, ]), drop = FALSE]
		pw <- rbind(pw, data.frame(case = cs$lab, N = N, yes = mean(s[1L, ] >= 1), yes_capped = mean(s[1L, ] >= 1 & s[2L, ] == 0), sparse = mean(s[2L, ] == 1)))
	}
}
print(pw, digits = 3)
saveRDS(list(h0 = out, h1 = pw), "dev/gbay/scripts/sim_bilan.rds")
