# sys_bay_fact_laws.R
# Variances per cell with several factors between groups (m.bay, decision D19). Two laws:
# normal with a common variance (lmBF, sys_bay_factorial.R) and normal with one variance per
# cell. Each hereditary model gets a marginal likelihood under both laws, on the measure of
# sys_bay_laws.R; the Bayes factor of each model against the null is then averaged over the
# laws, which leaves the model space and its prior unchanged. Nothing is exported.
#
# Parametrisation (the cells play the part of the groups of sys_bay_laws.R):
#   y_ij ~ N(mu + x_c' beta, sigma_c^2), sigma_c^2 = 1 / (w_c taubar), w_c = K vartheta_c
#   p(mu, taubar) prop 1 / taubar, vartheta ~ Dirichlet(u, ..., u) on the K cells
#   beta_t | g_t ~ N(0, g_t / taubar I), g_t ~ InvGamma(1/2, r^2 / 2), one g per term
# x_c: the coding of BayesFactor (orthonormal contrasts per factor, products of them for an
# interaction), so that w = 1 gives back the Zellner-Siow prior of lmBF (Rouder et al. 2012);
# the result does not depend on the choice of the orthonormal basis, the prior being
# isotropic in it. 1 / taubar is the harmonic mean of the cell variances: the effects keep a
# common scale, which an additive model needs (one effect of A spans cells of different
# spreads). The price is defect B2 of the one-factor case: a very dispersed cell changes
# that scale a little. It is reported in DECISIONS, D19.
#
# Integrating mu, taubar and beta, given (w, g):
#   log m = C + sum (n_c - 1)/2 log w_c - 1/2 log|V| - 1/2 log(1'V^-1 1) - (N - 1)/2 log Q
#   V = diag(1 / (n_c w_c)) + X D X',  Q = sum w_c SS_c + m'V^-1 m - (1'V^-1 m)^2 / 1'V^-1 1
# with C that of .bay_law_const, so that both laws stand on one measure. The integral on
# (vartheta, g) is done by importance sampling around the mode (multivariate t, df 4) with a
# defensive part (the proposal of sys_bay_laws.R for vartheta, the prior for g).
#
# References: Dablander et al. (2024) for the Dirichlet on the proportions of precision;
# Rouder et al. (2012) for the prior on the effects (see sys_bay_laws.R).


# Orthonormal basis of the contrasts of a factor with a levels (a x (a - 1)).
.bay_fl_contr <- function(a) {
	e <- eigen(diag(a) - 1 / a, symmetric = TRUE)$vectors[, seq_len(a - 1L), drop = FALSE]
	e
}

# Design of the cells: one matrix per term (K x p_t), K the non-empty cells.
.bay_fl_design <- function(factors, keys) {
	# Same order of the cells as .bay_laws_stats, which sorts the labels.
	cl <- as.character(interaction(factors, drop = TRUE))
	first <- match(levels(factor(cl)), cl)
	lev <- lapply(factors, function(f) as.integer(f[first]))
	lapply(strsplit(keys, ":", fixed = TRUE), function(fs) {
		X <- NULL
		for (v in fs) {
			P <- .bay_fl_contr(nlevels(factors[[v]]))[lev[[v]], , drop = FALSE]
			X <- if (is.null(X)) P else {
				do.call(cbind, lapply(seq_len(ncol(P)), function(j) X * P[, j]))
			}
		}
		X
	})
}

# Kernel for S draws: W (S x K) weights, G (S x T) g per term of the model (T may be 0).
.bay_fl_kernel <- function(st, Xs, W, G) {
	S <- nrow(W)
	K <- st$k
	n <- st$n
	m <- st$m
	one <- rep(1, K)
	out <- numeric(S)
	XX <- lapply(Xs, tcrossprod)
	lw <- log(W)
	base <- drop(lw %*% ((n - 1) / 2))
	wss <- drop(W %*% st$ss)
	for (s in seq_len(S)) {
		V <- diag(1 / (n * W[s, ]), K)
		for (t in seq_along(XX)) V <- V + G[s, t] * XX[[t]]
		R <- tryCatch(chol(V), error = function(e) NULL)
		if (is.null(R)) {
			out[s] <- -Inf
			next
		}
		a <- backsolve(R, forwardsolve(t(R), cbind(one, m)))
		v11 <- sum(one * a[, 1L])
		v1m <- sum(m * a[, 1L])
		vmm <- sum(m * a[, 2L])
		Q <- wss[s] + vmm - v1m^2 / v11
		out[s] <- base[s] - sum(log(diag(R))) - 0.5 * log(v11) - (st$N - 1) / 2 * log(Q)
	}
	out
}

.bay_fl_ligp <- function(g, r) 0.5 * log(r^2 / 2) - lgamma(0.5) - 1.5 * log(g) - r^2 / (2 * g)

# Log marginal likelihood of one model under the law with one variance per cell, for each r.
# Xs: the design of the terms of the model. Returns log_ml (one per r), ess, rel_se.
.bay_fl_het_ml <- function(st, Xs, rs, u, draws) {
	K <- st$k
	nt <- length(Xs)
	r0 <- rs[1L]
	ldir <- function(th) lgamma(K * u) - K * lgamma(u) + (u - 1) * rowSums(log(th))
	split_z <- function(Z) {
		th <- .bay_alr_inv(Z[, seq_len(K - 1L), drop = FALSE])
		G <- if (nt) exp(Z[, K - 1L + seq_len(nt), drop = FALSE]) else matrix(0, nrow(Z), 0L)
		list(th = th, G = G)
	}
	# log target in z (alr of vartheta, log g), at the central r; per-r terms added later.
	lt_parts <- function(Z) {
		p <- split_z(Z)
		lk <- .bay_fl_kernel(st, Xs, K * p$th, p$G)
		lk + ldir(p$th) + rowSums(log(p$th))
	}
	lg_prior <- function(Z, r) {
		if (!nt) return(rep(0, nrow(Z)))
		lg <- Z[, K - 1L + seq_len(nt), drop = FALSE]
		rowSums(.bay_fl_ligp(exp(lg), r) + lg)
	}
	# Mode: start from the own variances of the cells and g = r^2.
	th_start <- (1 / pmax(st$ss / pmax(st$n - 1, 1), 1e-3))
	th_start <- th_start / sum(th_start)
	z0 <- c(log(th_start[-K]) - log(th_start[K]), rep(log(r0^2), nt))
	f <- function(z) {
		v <- lt_parts(matrix(z, 1L)) + lg_prior(matrix(z, 1L), r0)
		if (!is.finite(v)) -1e10 else v
	}
	opt <- stats::optim(z0, f, method = "BFGS", control = list(fnscale = -1, maxit = 200))
	H <- tryCatch(stats::optimHess(opt$par, f), error = function(e) NULL)
	d <- length(z0)
	ok <- !is.null(H) && all(is.finite(H))
	if (ok) {
		ev <- eigen(-H, symmetric = TRUE)
		ok <- all(ev$values > 0)
	}
	if (ok) {
		fit <- .bay_mvt_make(opt$par, 1 / pmax(ev$values, 1e-8), ev$vectors)
	} else {
		fit <- .bay_mvt_make(opt$par, rep(1, d), diag(d))
	}
	pr <- .bay_law_proposal(st, u)
	mix <- c(0.15, 0.85)
	ldef <- function(Z) {
		p <- split_z(Z)
		.bay_law_lq_def(pr, p$th) + rowSums(log(p$th)) + lg_prior(Z, r0)
	}
	draw <- function(S) {
		ns <- as.numeric(stats::rmultinom(1L, S, mix))
		Zd <- NULL
		if (ns[1L] > 0L) {
			th <- .bay_law_draw_def(pr, ns[1L])
			lg <- if (nt) log(matrix(r0^2 / 2 / stats::rgamma(ns[1L] * nt, 0.5), ns[1L], nt)) else NULL
			Zd <- cbind(.bay_alr(th), lg)
		}
		Zt <- if (ns[2L] > 0L) .bay_mvt_draw(ns[2L], fit, 4, 1.2) else NULL
		rbind(Zd, Zt)
	}
	lq <- function(Z) .bay_row_lse(cbind(log(mix[1L]) + ldef(Z), log(mix[2L]) + .bay_mvt_ld(Z, fit, 4, 1.2)))
	Z <- draw(draws)
	# Refit once on the weighted draws when the mode-based fit is poor.
	base <- lt_parts(Z)
	lw0 <- base + lg_prior(Z, r0) - lq(Z)
	if (!isTRUE(.bay_law_is(lw0)$ess >= draws / 5)) {
		fit <- .bay_mvt_fit(Z, lw0)
		Z <- draw(draws)
		base <- lt_parts(Z)
	}
	lqZ <- lq(Z)
	res <- lapply(rs, function(r) .bay_law_is(base + lg_prior(Z, r) - lqZ))
	list(log_ml = .bay_law_const(st) + vapply(res, `[[`, numeric(1), "log_mean"),
		ess = min(vapply(res, `[[`, numeric(1), "ess")),
		rel_se = vapply(res, `[[`, numeric(1), "rel_se"))
}

# The same with the common variance (w = 1), for the check against lmBF only: integral on g.
.bay_fl_eq_ml <- function(st, Xs, rs, draws) {
	K <- st$k
	nt <- length(Xs)
	W1 <- matrix(1, 1L, K)
	if (!nt) return(.bay_law_const(st) + .bay_fl_kernel(st, Xs, W1, matrix(0, 1L, 0L)))
	vapply(rs, function(r) {
		lg <- matrix(log(r^2 / 2 / stats::rgamma(draws * nt, 0.5)), draws, nt)
		lk <- .bay_fl_kernel(st, Xs, matrix(1, draws, K), exp(lg))
		.bay_law_const(st) + .bay_law_is(lk)$log_mean
	}, numeric(1))
}

# Averaged log BF of every model against the null: log_bf_eq (models x r) from lmBF, eq0 the
# log marginal of the null with the common variance, l_het (models x r) with one variance per
# cell. Equal prior on the two laws, as with one factor. The model space and its prior are
# those of the common variance: only the Bayes factors change.
.bay_fl_combine <- function(log_bf_eq, eq0, l_het) {
	l_eq <- log_bf_eq + eq0
	avg <- l_eq
	avg[] <- .bay_lse2(as.vector(l_eq), as.vector(l_het))
	list(log_bf = sweep(avg, 2L, avg[1L, ], "-"), l_eq = l_eq)
}

# P(one variance per cell | data), the models weighted by their prior.
.bay_fl_p_het <- function(l_eq, l_het, log_prior) {
	a <- .bay_log_sum_exp(log_prior + l_het)
	b <- .bay_log_sum_exp(log_prior + l_eq)
	1 / (1 + exp(b - a))
}

# Marginal likelihoods with one variance per cell, for every model and every r.
.bay_fl_het_all <- function(y, factors, sp, rs, draws, seed, u = 1.5) {
	cl <- interaction(factors, drop = TRUE)
	st <- .bay_laws_stats(y, cl)
	Xall <- .bay_fl_design(factors, sp$keys)
	het <- .bay_with_seed(seed, lapply(sp$models, function(m) {
		.bay_fl_het_ml(st, Xall[m], rs, u, draws)
	}))
	l_het <- matrix(t(vapply(het, `[[`, numeric(length(rs)), "log_ml")), ncol = length(rs))
	ess <- vapply(het, `[[`, numeric(1), "ess")
	eq0 <- .bay_law_const(st) + .bay_fl_kernel(st, list(), matrix(1, 1L, st$k), matrix(0, 1L, 0L))
	failed <- anyNA(l_het) || any(!is.finite(l_het)) || any(!(ess >= .bay_law_min_ess()))
	list(l_het = l_het, eq0 = eq0, ess = ess, failed = failed, stats = st)
}
