# Validation de g.bay() par simulation (dev/gbay/DECISIONS.md, G9).
# Lancer depuis la racine du paquet :
#   Rscript dev/gbay/scripts/sim_validation.R
# Sorties : dev/gbay/scripts/sim_validation.rds et un resume imprime.
# Les facteurs de Bayes sont ceux de g.bay (fonctions internes appelees directement pour la
# vitesse ; leur egalite avec g.bay() est verifiee en tete).

.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(pkgload::load_all(".", quiet = TRUE))
set.seed(20261111)
R <- as.integer(Sys.getenv("GBAY_REPS", "2000"))
t0 <- Sys.time()

# Controle : la voie rapide = g.bay().
tt <- matrix(c(12, 5, 7, 14, 3, 9), 2)
stopifnot(abs(.gbay_lbf_ct(tt, "jointMulti") - g.bay(tt, verbose = FALSE, sensitivity = FALSE)$log_bf) < 1e-12)
stopifnot(abs(.gbay_lbf_gof(c(5, 9, 2), rep(1 / 3, 3)) - g.bay(c(5, 9, 2), verbose = FALSE, sensitivity = FALSE)$log_bf) < 1e-12)

# Voie rapide : les fonctions internes que contingencyTableBF appelle (22 ms par appel pour
# contingencyTableBF, 0,15 ms pour elles) et une statistique G vectorisee (celle de GTest, sans
# correction). Leur egalite avec contingencyTableBF et GTest est verifiee ici.
.lbf_joint <- function(t) BayesFactor:::contingencyJointMultinomial(t, 1)
.lbf_rows <- function(t) BayesFactor:::contingencyIndepMultinomial(t, 1)
.g_p <- function(x, p = NULL) {
	if (is.null(p)) {
		E <- outer(rowSums(x), colSums(x)) / sum(x); df <- (nrow(x) - 1) * (ncol(x) - 1)
	} else {
		E <- sum(x) * p; df <- length(x) - 1
	}
	g <- 2 * sum(ifelse(x > 0, x * log(x / E), 0))
	stats::pchisq(g, df, lower.tail = FALSE)
}
for (i in 1:20) {
	tt <- matrix(stats::rpois(6, 8) + 1, 2)
	stopifnot(abs(.lbf_joint(tt) - .gbay_lbf_ct(tt, "jointMulti")) < 1e-10, abs(.lbf_rows(tt) - .gbay_lbf_ct(tt, "indepMulti", "rows")) < 1e-10,
		abs(.g_p(tt) - suppressWarnings(DescTools::GTest(tt)$p.value)) < 1e-10)
	xx <- stats::rpois(4, 5); pp <- c(.4, .3, .2, .1)
	stopifnot(abs(.g_p(xx, pp) - suppressWarnings(DescTools::GTest(xx, p = pp)$p.value)) < 1e-10)
}

pr <- function(lb) 1 / (1 + exp(-lb))
rates <- function(lb, pv) c(p75 = mean(pr(lb) > 0.75), p95 = mean(pr(lb) > 0.95), bf3 = mean(lb > log(3)),
	gtest05 = mean(pv < 0.05), p_h0_75 = mean(pr(lb) < 0.25))

gtest_p <- function(x, p = NULL) {
	if (is.null(p)) suppressWarnings(DescTools::GTest(x)$p.value) else suppressWarnings(DescTools::GTest(x, p = p)$p.value)
}

# --- 1. Ajustement sous H0 ----------------------------------------------------------------
gof_h0 <- NULL
for (pp in list(c(.5, .5), c(.75, .25), rep(1 / 3, 3), c(.6, .3, .1), rep(.2, 5), c(.4, .25, .15, .1, .1))) {
	for (N in c(10, 20, 50, 200, 1000)) {
		X <- stats::rmultinom(R, N, pp)
		lb <- apply(X, 2L, .gbay_lbf_gof, p = pp, a = 1)
		pv <- apply(X, 2L, .g_p, p = pp)
		gof_h0 <- rbind(gof_h0, data.frame(k = length(pp), p = paste(round(pp, 2), collapse = "/"), N = N, t(rates(lb, pv))))
	}
}
cat("\n== Ajustement sous H0 (taux de faux oui) ==\n"); print(gof_h0, digits = 3)

# --- 2. Independance sous H0 ----------------------------------------------------------------
# Plan joint (N fixe) analyse par le plan auto (joint) ; plan a lignes fixees analyse par
# y ~ g (independant, lignes) et par x, y (joint, plan mal choisi).
ind_h0 <- NULL
designs <- list(list(r = c(.5, .5), c = c(.5, .5)), list(r = c(.7, .3), c = c(.2, .8)),
	list(r = c(.5, .5), c = c(.5, .3, .2)), list(r = c(.4, .35, .25), c = c(.1, .3, .3, .3)))
for (ds in designs) {
	for (N in c(10, 20, 50, 200, 1000)) {
		P <- outer(ds$r, ds$c)
		Rk <- if (N >= 1000) R / 2 else R
		X <- stats::rmultinom(Rk, N, as.numeric(P))
		res <- lapply(seq_len(Rk), function(i) {
			t <- matrix(X[, i], length(ds$r))
			t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]
			if (nrow(t) < 2L || ncol(t) < 2L) return(c(NA, NA))
			c(.lbf_joint(t), .g_p(t))
		})
		m <- do.call(rbind, res); m <- m[!is.na(m[, 1L]), , drop = FALSE]
		ind_h0 <- rbind(ind_h0, data.frame(table = paste0(length(ds$r), "x", length(ds$c)), margins = if (max(P) / min(P) > 3) "uneven" else "even",
			sampling = "joint/joint", N = N, n_ok = nrow(m), t(rates(m[, 1L], m[, 2L]))))
		# Lignes fixees (groupes de meme taille), analyse independante et analyse jointe.
		nr <- length(ds$r); ng <- max(2L, round(N / nr))
		res <- lapply(seq_len(Rk), function(i) {
			t <- t(stats::rmultinom(nr, ng, ds$c))
			t <- t[, colSums(t) > 0, drop = FALSE]
			if (ncol(t) < 2L) return(c(NA, NA, NA))
			c(.lbf_rows(t), .lbf_joint(t), .g_p(t))
		})
		m <- do.call(rbind, res); m <- m[!is.na(m[, 1L]), , drop = FALSE]
		ind_h0 <- rbind(ind_h0, data.frame(table = paste0(nr, "x", length(ds$c)), margins = if (max(P) / min(P) > 3) "uneven" else "even",
			sampling = "rows/indep", N = ng * nr, n_ok = nrow(m), t(rates(m[, 1L], m[, 3L]))),
			data.frame(table = paste0(nr, "x", length(ds$c)), margins = if (max(P) / min(P) > 3) "uneven" else "even",
				sampling = "rows/joint", N = ng * nr, n_ok = nrow(m), t(rates(m[, 2L], m[, 3L]))))
	}
}
cat("\n== Independance sous H0 (taux de faux oui) ==\n"); print(ind_h0, digits = 3)

# --- 3. Puissance, comparee au test G -----------------------------------------------------
pw <- NULL
# Independance 2 x 3, plan joint ; ecart de taille w (w de Cohen = V * sqrt(min(r,c)-1)).
base_r <- c(.5, .5); base_c <- c(.4, .35, .25)
D <- outer(c(1, -1), c(1, 0, -1))
for (w in c(0.1, 0.2, 0.3)) {
	P0 <- outer(base_r, base_c)
	# delta tel que sum((P-P0)^2/P0) = w^2
	dlt <- w / sqrt(sum(D^2 / P0))
	P <- P0 + dlt * D
	for (N in c(50, 100, 200, 400)) {
		X <- stats::rmultinom(R / 2, N, as.numeric(P))
		m <- t(apply(X, 2L, function(x) {
			t <- matrix(x, 2); t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]
			if (nrow(t) < 2L || ncol(t) < 2L) return(c(NA, NA))
			c(.lbf_joint(t), .g_p(t))
		}))
		m <- m[!is.na(m[, 1L]), , drop = FALSE]
		pw <- rbind(pw, data.frame(case = "indep 2x3", w = w, N = N, t(rates(m[, 1L], m[, 2L]))))
	}
}
# Ajustement a 4 modalites, p egales ; ecart de taille w.
p0 <- rep(.25, 4); Dg <- c(1, 1, -1, -1)
for (w in c(0.1, 0.2, 0.3)) {
	p1 <- p0 + w / sqrt(sum(Dg^2 / p0)) * Dg
	for (N in c(50, 100, 200, 400)) {
		X <- stats::rmultinom(R / 2, N, p1)
		lb <- apply(X, 2L, .gbay_lbf_gof, p = p0, a = 1)
		pv <- apply(X, 2L, .g_p, p = p0)
		pw <- rbind(pw, data.frame(case = "gof k=4", w = w, N = N, t(rates(lb, pv))))
	}
}
cat("\n== Puissance (taux de oui sous H1) ==\n"); print(pw, digits = 3)

# --- 4. Cases signalees sous H0 (G6) ------------------------------------------------------
Rc <- as.integer(Sys.getenv("GBAY_CELL_REPS", "300"))
cells_h0 <- NULL
for (ds in list(list(r = c(.5, .5), c = c(.5, .5)), list(r = c(.5, .5), c = c(.5, .3, .2)), list(r = c(.4, .35, .25), c = c(.1, .3, .3, .3)))) {
	for (N in c(50, 200)) {
		P <- outer(ds$r, ds$c)
		X <- stats::rmultinom(Rc, N, as.numeric(P))
		fl <- vapply(seq_len(Rc), function(i) {
			t <- matrix(X[, i], length(ds$r), dimnames = list(paste0("r", seq_along(ds$r)), paste0("c", seq_along(ds$c))))
			t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]
			if (nrow(t) < 2L || ncol(t) < 2L) return(c(NA, NA))
			s <- .gbay_summ_ct(.gbay_draws_ct(t, "jointMulti", NULL, 2000L, i), t)
			c(any(s$cells$reading != "-"), mean(s$cells$reading != "-"))
		}, numeric(2))
		fl <- fl[, !is.na(fl[1L, ]), drop = FALSE]
		cells_h0 <- rbind(cells_h0, data.frame(table = paste0(length(ds$r), "x", length(ds$c)), N = N,
			tables_with_flag = mean(fl[1L, ]), cells_flagged = mean(fl[2L, ])))
	}
}
cat("\n== Cases signalees sous H0 ==\n"); print(cells_h0, digits = 3)

# --- 5. Couverture de l'intervalle du rapport des cotes (2 x 2, sous H1) --------------------
cov_or <- NULL
for (N in c(40, 100, 400)) {
	P <- matrix(c(.3, .2, .15, .35), 2); true_lor <- log(P[1, 1] * P[2, 2] / (P[1, 2] * P[2, 1]))
	X <- stats::rmultinom(Rc, N, as.numeric(P))
	hit <- vapply(seq_len(Rc), function(i) {
		t <- matrix(X[, i], 2, dimnames = list(c("a", "b"), c("u", "v")))
		if (any(rowSums(t) == 0) || any(colSums(t) == 0)) return(NA)
		o <- .gbay_summ_ct(.gbay_draws_ct(t, "jointMulti", NULL, 2000L, i), t)$odds_ratio
		log(o[["low"]]) <= true_lor && true_lor <= log(o[["high"]])
	}, logical(1))
	cov_or <- rbind(cov_or, data.frame(N = N, coverage = mean(hit, na.rm = TRUE), n_ok = sum(!is.na(hit))))
}
cat("\n== Couverture de l'intervalle a 95 % du rapport des cotes ==\n"); print(cov_or, digits = 3)

saveRDS(list(gof_h0 = gof_h0, ind_h0 = ind_h0, power = pw, cells_h0 = cells_h0, cov_or = cov_or, R = R, Rc = Rc,
	seed = 20261111, date = Sys.time()), "dev/gbay/scripts/sim_validation.rds")
cat("\nDuree :", format(Sys.time() - t0), "\n")
