# Recalibration de l'independance apres G11 (test log-lineaire de l'interaction, BIC avec n = N).
#   Rscript dev/gbay/scripts/sim_loglin.R
# Grille demandee par la session parente : 2 x 2, 3 x 4, 4 x 5 ; N = 30, 60, 200 ; marges egales
# et inegales ; 2000 tables par scenario. Comparaison avec l'ancien moteur (contingencyTableBF,
# trois plans) et avec le test G. Puis puissance sous un lien de taille w.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(pkgload::load_all(".", quiet = TRUE))
set.seed(20261112)
R <- as.integer(Sys.getenv("GBAY_REPS", "2000"))
t0 <- Sys.time()

# Voie rapide : G vectorise ; egalite avec g.bay() (glm) et GTest verifiee ici.
.G <- function(t) {
	E <- outer(rowSums(t), colSums(t)) / sum(t)
	2 * sum(ifelse(t > 0, t * log(t / E), 0))
}
.lb_new <- function(t) (.G(t) - (nrow(t) - 1) * (ncol(t) - 1) * log(sum(t))) / 2
for (i in 1:20) {
	tt <- matrix(stats::rpois(12, 6), 3); tt <- tt[rowSums(tt) > 0, colSums(tt) > 0, drop = FALSE]
	stopifnot(abs(.lb_new(tt) - .gbay_lbf_ct(tt)) < 1e-8,
		abs(.G(tt) - unname(suppressWarnings(DescTools::GTest(tt)$statistic))) < 1e-8)
}
J <- function(t) BayesFactor:::contingencyJointMultinomial(t, 1)
Pn <- function(t) BayesFactor:::contingencyPoisson(t, 1)
Ir <- function(t) BayesFactor:::contingencyIndepMultinomial(t, 1)
pr <- function(l, p = 0.5) 1 / (1 + exp(-l) * (1 - p) / p)
# Ce que dit le bilan avec G11 : P >= 0,75, BF >= 3, et l'ordre ne change pas avec les chances
# de depart 1:2 a 2:1 (BF >= 2) ; « tres probable » a P >= 0,95.
said <- function(lb) ifelse(lb >= log(3) & pr(lb, 1 / 3) > 0.5, ifelse(pr(lb) >= 0.95, 2L, 1L), 0L)

clean <- function(x, r) { t <- matrix(x, r); t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE] }
h0 <- NULL
for (dm in list(c(2, 2), c(3, 4), c(4, 5))) for (mg in c("egales", "inegales")) for (N in c(30, 60, 200)) {
	rr <- if (mg == "egales") rep(1, dm[1]) else seq(1, 3, length.out = dm[1])
	cc <- if (mg == "egales") rep(1, dm[2]) else seq(1, 4, length.out = dm[2])
	P <- outer(rr / sum(rr), cc / sum(cc))
	X <- stats::rmultinom(R, N, as.numeric(P))
	m <- vapply(seq_len(R), function(i) {
		t <- clean(X[, i], dm[1])
		if (nrow(t) < 2L || ncol(t) < 2L) return(rep(NA_real_, 6))
		df <- (nrow(t) - 1) * (ncol(t) - 1)
		c(.lb_new(t), J(t), Pn(t), Ir(t), stats::pchisq(.G(t), df, lower.tail = FALSE), said(.lb_new(t)))
	}, numeric(6))
	m <- m[, !is.na(m[1L, ]), drop = FALSE]
	h0 <- rbind(h0, data.frame(table = paste0(dm[1], "x", dm[2]), margins = mg, N = N, n_ok = ncol(m),
		new_bf3 = mean(m[1L, ] > log(3)), new_bf20 = mean(m[1L, ] > log(20)), new_bilan_yes = mean(m[6L, ] >= 1), new_bilan_very = mean(m[6L, ] == 2),
		old_joint_bf20 = mean(m[2L, ] > log(20)), old_poisson_bf20 = mean(m[3L, ] > log(20)), old_rows_bf20 = mean(m[4L, ] > log(20)),
		gtest05 = mean(m[5L, ] < 0.05)))
}
cat("== Independance sous H0, G11 contre l'ancien moteur ==\n"); print(h0, digits = 3, width = 220)

# Puissance : 2 x 3 et 3 x 4, ecart de taille w (w de Cohen), plan joint.
pw <- NULL
mk <- function(r, c, D, w) { P0 <- outer(r, c); P0 + w / sqrt(sum(D^2 / P0)) * D }
cases <- list(list(lab = "2x2", r = c(.5, .5), c = c(.5, .5), D = outer(c(1, -1), c(1, -1))),
	list(lab = "2x3", r = c(.5, .5), c = c(.4, .35, .25), D = outer(c(1, -1), c(1, 0, -1))),
	list(lab = "3x4", r = c(.4, .35, .25), c = c(.1, .3, .3, .3), D = outer(c(1, 0, -1), c(0, 1, 0, -1))))
for (cs in cases) for (w in c(0.1, 0.3, 0.5)) for (N in c(30, 60, 200)) {
	P <- mk(cs$r, cs$c, cs$D, w)
	if (any(P <= 0)) next
	X <- stats::rmultinom(R / 2, N, as.numeric(P))
	m <- vapply(seq_len(R / 2), function(i) {
		t <- clean(X[, i], length(cs$r))
		if (nrow(t) < 2L || ncol(t) < 2L) return(rep(NA_real_, 4))
		df <- (nrow(t) - 1) * (ncol(t) - 1)
		c(.lb_new(t), J(t), stats::pchisq(.G(t), df, lower.tail = FALSE), said(.lb_new(t)))
	}, numeric(4))
	m <- m[, !is.na(m[1L, ]), drop = FALSE]
	pw <- rbind(pw, data.frame(table = cs$lab, w = w, N = N, new_p75 = mean(pr(m[1L, ]) > 0.75), new_p95 = mean(pr(m[1L, ]) > 0.95),
		new_bilan_yes = mean(m[4L, ] >= 1), old_joint_p75 = mean(pr(m[2L, ]) > 0.75), gtest05 = mean(m[3L, ] < 0.05)))
}
cat("\n== Puissance (taux de oui sous un lien) ==\n"); print(pw, digits = 3, width = 220)
saveRDS(list(h0 = h0, power = pw, R = R, seed = 20261112), "dev/gbay/scripts/sim_loglin.rds")
cat("\nDuree :", format(Sys.time() - t0), "\n")
