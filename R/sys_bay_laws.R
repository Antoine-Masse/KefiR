# sys_bay_laws.R
# Comparison of sampling laws for the k-group plan (m.bay, decision D6). Eight
# models: {normal, Student t} x {common variance, one variance per group} x
# {no effect, effect}. Each gets a marginal likelihood under the same measure on
# the shared parameters; answers are averaged on their posterior probabilities.
# No preliminary test and no p-value here. Nothing is exported.
#
# Parametrisation (dev/mbay/METHODOLOGIE.md, section 7):
#   y_ij ~ law(theta_i, sigma_i^2), sigma_i^2 = c_i / taubar, c_i = 1 / (k vartheta_i)
#   p(mu, taubar) prop 1 / taubar         shared by all the models
#   vartheta ~ Dirichlet(u, ..., u)       variances per group; c_i = 1 otherwise
#   H1, an even mixture of two forms:
#     all groups: theta_i = mu + a_i, a_i | g ~ N(0, g sigma_i^2), g ~ InvGamma(1/2, r^2/2)
#     one group, s uniform: theta_s = mu + a_s, a_s | g ~ N(0, g sigma_s^2),
#                g ~ InvGamma(1/2, (r sqrt(2))^2 / 2), the others at mu
#   H0: theta_i = mu
#   t: scale s_i with s_i^2 nu / (nu - 2) = sigma_i^2, p(nu) = exp(-(nu - 2)), nu > 2
# With all c_i = 1 and the normal law, the first form is exactly the Zellner-Siow prior
# of Rouder et al. (2012), hence .bay_bf_oneway_zs, and each piece of the second is that
# same Bayes factor on the factor "s against the rest". At k = 2 the mixture is not used:
# both forms coincide when the variances are equal.
#
# Why the second form: with a common g, one far group costs d^(-k) under the first form,
# against d^(-n_s - 2u) for "no effect, that group is just more dispersed". For
# n_s + 2u <= k the evidence then never grows with d, and H0 can lead with 83 % at
# k = 20, Welch p = 1e-111 (dev/mbay/RELECTURE_ADVERSE_2026-11-07.md, B1). The second form
# costs d^(-2) whatever k. Specific to KefiR; the idea of comparing partitions of the
# means is that of Gopalan & Berry (1998).
#
# Why a_i | g is scaled by its own sigma_i, and not by a common scale: this is
# defect B2 of RELECTURE_ADVERSE_2026-11-06.md. An unrelated dispersed group must
# not change what a difference between two other groups means. Specific to KefiR
# for k >= 3: no source consulted combines these pieces.
#
# References (comments only, never in messages: dev/bp.log, 7.4 and BP-008).
#   Dablander, F., van den Bergh, D., Wagenmakers, E.-J., & Ly, A. (2024). Default
#     Bayes factors for testing the (in)equality of several population variances.
#     Bayesian Analysis, 19(3), 699-723. https://doi.org/10.1214/23-BA1369
#     (Dirichlet on the proportions of precision, nesting with the right Haar
#     prior on (mu, taubar); section 2.2 read on the arXiv PDF.)
#   Maier, M., Bartos, F., Quintana, D. S., Dablander, F., van den Bergh, D.,
#     Marsman, M., Ly, A., & Wagenmakers, E.-J. (2025). Model-averaged Bayesian t
#     tests. Psychonomic Bulletin & Review, 32(3), 1007-1031.
#     https://doi.org/10.3758/s13423-024-02590-5
#     (Eight models, Beta(1.5, 1.5) on the proportion of precision, shifted
#     exponential prior on nu; two groups only.)
#   Rouder, J. N., Morey, R. D., Speckman, P. L., & Province, J. M. (2012).
#     Default Bayes factors for ANOVA designs. Journal of Mathematical Psychology,
#     56(5), 356-374. https://doi.org/10.1016/j.jmp.2012.08.001


# ---------------------------------------------------------------------------
# Preparation
# ---------------------------------------------------------------------------

# Standardises y (grand mean 0, pooled within-group variance 1). Every model is
# computed on the same standardised data with the same improper measure, so the
# ratios do not depend on this choice.
.bay_laws_stats <- function(y, g) {
	keep <- !(is.na(y) | is.na(g))
	y <- as.numeric(y[keep])
	g <- factor(as.character(g[keep]))
	st0 <- .bay_oneway_prep(y, g)
	if (st0$N <= st0$k) {
		.exit("No residual degree of freedom (N <= k).",
			"Aucun degr\u00e9 de libert\u00e9 r\u00e9siduel (N <= k).")
	}
	if (isTRUE(sum(st0$ss) == 0)) {
		if (isTRUE(sum(st0$n * st0$m^2) <= st0$N * .Machine$double.eps^2)) {
			.exit("y is constant: the probabilities are undetermined.",
				"y est constant : les probabilit\u00e9s sont ind\u00e9termin\u00e9es.")
		}
		.exit("No variation within any group: the normal and t laws do not apply.",
			"Aucune variation dans aucun groupe : les lois normale et t ne s'appliquent pas.")
	}
	# Scaled first: squares of y * 1e200 or y * 1e-200 would overflow or vanish.
	y1 <- y / max(abs(y))
	by_group <- split(y1, g)
	n <- as.numeric(lengths(by_group))
	N <- sum(n)
	unit <- sqrt(sum(vapply(by_group, function(v) sum((v - mean(v))^2), numeric(1))) /
		(N - length(n)))
	y_std <- (y1 - mean(y1)) / unit
	by_group <- split(y_std, g)
	list(k = length(n), n = n, N = N, levels = levels(g), y = by_group,
		y_raw = y, g_raw = g,
		m = vapply(by_group, mean, numeric(1), USE.NAMES = FALSE),
		ss = vapply(by_group, function(v) sum((v - mean(v))^2), numeric(1), USE.NAMES = FALSE),
		unit = unit)
}


# ---------------------------------------------------------------------------
# Normal laws: closed form given (c, g)
# ---------------------------------------------------------------------------

# Integrating theta, mu and taubar gives, for C the constant below,
#   log m(y | c, g) = C - 1/2 sum (n_i - 1) log c_i - 1/2 sum log v_i
#                     - 1/2 log W - (N - 1)/2 log Q
#   v_i = c_i (1/n_i + g), W = sum 1/v_i,
#   Q = sum SS_i / c_i + sum (m_i - mbar_w)^2 / v_i.
# g = 0 is H0. The constant keeps normal and t on one measure.
.bay_law_const <- function(st) {
	nu <- st$N - 1
	-nu / 2 * log(2 * pi) - 0.5 * sum(log(st$n)) + lgamma(nu / 2) + nu / 2 * log(2)
}

.bay_law_kernel <- function(st, cc, g) {
	n <- st$n
	S <- nrow(cc)
	k <- st$k
	lc <- log(cc)
	v <- sweep(cc, 2L, 1 / n + g, "*")
	w <- 1 / v
	W <- rowSums(w)
	M <- matrix(rep(st$m, each = S), S, k)
	mbar <- rowSums(w * M) / W
	Q <- rowSums(sweep(1 / cc, 2L, st$ss, "*")) + rowSums(w * (M - mbar)^2)
	-0.5 * rowSums(sweep(lc, 2L, n - 1, "*")) - 0.5 * rowSums(log(v)) - 0.5 * log(W) -
		(st$N - 1) / 2 * log(Q)
}

# The same on a vector of g: one column per g. No loop on g, matrix products only.
.bay_law_kernel_grid <- function(st, cc, gs, s = NULL) {
	n <- st$n
	k <- st$k
	S <- nrow(cc)
	lc <- log(cc)
	w0 <- 1 / cc
	# Means relative to the group of largest weight on each row: S2 - Wm^2 / W then
	# stays a difference of comparable numbers, even at extreme c.
	M <- matrix(rep(st$m, each = S), S, k) - st$m[max.col(w0, ties.method = "first")]
	if (is.null(s)) {
		a <- 1 / outer(1 / n, gs, "+")
	} else {
		# Only group s has an effect: v_i = c_i / n_i elsewhere.
		a <- matrix(n, k, length(gs))
		a[s, ] <- 1 / (1 / n[s] + gs)
	}
	W <- w0 %*% a
	Wm <- (w0 * M) %*% a
	S2 <- (w0 * M^2) %*% a
	Q <- drop(w0 %*% st$ss) + pmax(S2 - Wm^2 / W, 0)
	base <- -0.5 * drop(lc %*% (n - 1)) - 0.5 * rowSums(lc)
	sweep(-0.5 * log(W) - (st$N - 1) / 2 * log(Q), 2L, 0.5 * colSums(log(a)), "+") + base
}

# Integral over g for each row of cc and each r, trapezoids in t = log g. Same
# argument as for .bay_bf_oneway_zs: analytic integrand, doubly exponential decay
# on the left, rate k/2 on the right. Left bound r^2 / (100 k): the prior peak is
# near r^2 / k (RELECTURE_ADVERSE_2026-11-06.md, I3). Right bound: the peak is near
# g n_i ~ N * (between part) / (within part), at most over the rows of cc.
.bay_law_gint <- function(st, cc, rs, h = min(0.25, 0.5 / sqrt(st$k)), s = NULL) {
	k <- st$k
	span <- diff(range(st$m))
	t_lo <- log(min(rs)^2 / (100 * k))
	ratio <- k^2 * span^2 * max(rowSums(1 / cc)) / min(drop((1 / cc) %*% st$ss) + 1e-300)
	t_hi <- max(log(max(rs)^2), log(k^2 * span^2 + 1e-8), log(st$N * ratio + 1e-8)) + 30
	S <- nrow(cc)
	# Pilot on a coarse grid and a few rows: keeps the part of the range where the
	# integrand is within exp(-40) of its peak on those rows, plus a margin.
	if (S > 50L) {
		pilot <- cc[round(seq(1, S, length.out = 50L)), , drop = FALSE]
		t_coarse <- seq(t_lo, t_hi, by = 1)
		Lp <- .bay_law_kernel_grid(st, pilot, exp(t_coarse), s)
		seen <- logical(length(t_coarse))
		for (r in range(rs)) {
			lp <- -0.5 * t_coarse - r^2 / (2 * exp(t_coarse))
			A <- sweep(Lp, 2L, lp, "+")
			seen <- seen | colSums(sweep(A, 1L, apply(A, 1L, max), "-") > -40) > 0L
		}
		t_lo <- max(t_lo, min(t_coarse[seen]) - 3)
		t_hi <- min(t_hi, max(t_coarse[seen]) + 3)
	}
	t_grid <- seq(t_lo, t_hi, by = h)
	L <- .bay_law_kernel_grid(st, cc, exp(t_grid), s)
	res <- matrix(NA_real_, S, length(rs))
	for (i in seq_along(rs)) {
		r <- rs[i]
		lp <- 0.5 * log(r^2 / 2) - 0.5 * log(pi) - 0.5 * t_grid - r^2 / (2 * exp(t_grid))
		A <- sweep(L, 2L, lp, "+")
		top <- A[cbind(seq_len(S), max.col(A, ties.method = "first"))]
		res[, i] <- top + log(rowSums(exp(A - top))) + log(h)
	}
	res
}


# log(exp(a) + exp(b)) by element.
.bay_lse2 <- function(a, b) {
	top <- pmax(a, b)
	top + log(exp(a - top) + exp(b - top))
}

# The even mixture of the two forms of H1, for each row of cc and each r. Returns the
# mixture in the first length(rs) columns, the "all groups" form alone in the others.
.bay_law_gint_mix <- function(st, cc, rs) {
	zs <- .bay_law_gint(st, cc, rs)
	if (st$k == 2L) return(cbind(zs, zs))
	sp <- NULL
	for (s in seq_len(st$k)) {
		x <- .bay_law_gint(st, cc, rs * sqrt(2), s = s)
		sp <- if (is.null(sp)) x else .bay_lse2(sp, x)
	}
	cbind(.bay_lse2(zs, sp - log(st$k)) + log(0.5), zs)
}

# ---------------------------------------------------------------------------
# Normal laws: importance sampling on the proportions of precision
# ---------------------------------------------------------------------------

# Log density of th = x / sum(x), x_i ~ Gamma(a_i, rate b_i), independent.
.bay_law_ldgd <- function(th, a, b) {
	lgamma(sum(a)) - sum(lgamma(a)) + sum(a * log(b)) + drop(log(th) %*% (a - 1)) -
		sum(a) * log(drop(th %*% b))
}

.bay_law_rgd <- function(S, a, b) {
	x <- vapply(seq_along(a), function(i) stats::rgamma(S, a[i], rate = b[i]), numeric(S))
	if (S == 1L) x <- matrix(x, nrow = 1L)
	x / rowSums(x)
}

# Defensive mixture. 1: the prior. 2: each group's own variance. 3: the same with
# the distance of its mean to the weighted median of the means, which is what H0
# needs when a group is far away (a plain mean would be pulled by that group, and
# would make the proposal wrong for all the others). 4: component 2 flattened.
.bay_law_proposal <- function(st, u) {
	k <- st$k
	nu <- st$n - 1
	s2 <- sum(st$ss) / (st$N - k)
	own_b <- ifelse(nu > 0 & st$ss > 0, st$ss / 2, (u + nu / 2) * s2)
	ord <- order(st$m)
	mb <- st$m[ord][which(cumsum(st$n[ord]) / st$N >= 0.5)[1L]]
	a2 <- u + nu / 2
	a4 <- pmax(a2 / 4, 0.1)
	comps <- list(
		list(a = rep(u, k), b = rep(1, k)),
		list(a = a2, b = own_b),
		list(a = u + st$n / 2, b = own_b + st$n * (st$m - mb)^2 / 2),
		list(a = a4, b = own_b * a4 / a2)
	)
	list(comps = comps, weights = c(0.1, 0.5, 0.2, 0.2))
}

.bay_law_draw_def <- function(pr, S) {
	ns <- as.numeric(stats::rmultinom(1L, S, pr$weights))
	th <- do.call(rbind, lapply(seq_along(pr$comps), function(j) {
		if (ns[j] == 0L) return(NULL)
		.bay_law_rgd(ns[j], pr$comps[[j]]$a, pr$comps[[j]]$b)
	}))
	th <- pmax(th, 1e-300)
	th / rowSums(th)
}

.bay_law_lq_def <- function(pr, th) {
	lq <- vapply(seq_along(pr$comps), function(j) {
		log(pr$weights[j]) + .bay_law_ldgd(th, pr$comps[[j]]$a, pr$comps[[j]]$b)
	}, numeric(nrow(th)))
	if (nrow(th) == 1L) lq <- matrix(lq, nrow = 1L)
	.bay_row_lse(lq)
}

# log(sum(exp(.))) per row.
.bay_row_lse <- function(A) {
	top <- A[cbind(seq_len(nrow(A)), max.col(A, ties.method = "first"))]
	top + log(rowSums(exp(A - top)))
}

# log mean of exp(lw), ESS and relative standard error.
.bay_law_is <- function(lw) {
	if (any(!is.finite(lw) & !(lw == -Inf), na.rm = TRUE) || anyNA(lw)) {
		return(list(log_mean = NA_real_, ess = NA_real_, rel_se = NA_real_))
	}
	top <- max(lw)
	if (!is.finite(top)) return(list(log_mean = NA_real_, ess = NA_real_, rel_se = NA_real_))
	w <- exp(lw - top)
	mw <- mean(w)
	list(log_mean = top + log(mw), ess = sum(w)^2 / sum(w^2),
		rel_se = stats::sd(w) / (mw * sqrt(length(w))))
}

# Absolute ceiling on the importance draws per model, whatever the growth.
.bay_law_max_draws <- function() 40000L

# An estimate on too few effective draws is not an estimate: its stated error saturates
# (rel_se tends to 1 when one weight dominates). Such a model counts as failed, and its
# law is dropped with a message, instead of taking all the weight by chance.
.bay_law_min_ess <- function() 30

# Multivariate t components. A fit is a centre and a square root of the covariance.
.bay_mvt_make <- function(centre, cov_vals, cov_vecs) {
	d <- length(centre)
	list(centre = centre, d = d,
		root = cov_vecs %*% diag(sqrt(cov_vals), d),
		root_inv = cov_vecs %*% diag(1 / sqrt(cov_vals), d),
		log_det_root = 0.5 * sum(log(cov_vals)))
}

# Weighted fit on a pilot. Weights are flattened when the pilot has few effective
# draws, so that the fit is not one or two points.
.bay_mvt_fit <- function(Z, lw) {
	d <- ncol(Z)
	w <- exp(lw - max(lw))
	for (i in 1:3) {
		if (sum(w)^2 / sum(w^2) >= 5 * (d + 1)) break
		w <- sqrt(w)
	}
	w <- w / sum(w)
	centre <- colSums(Z * w)
	D <- sweep(Z, 2L, centre, "-")
	eig <- eigen(crossprod(D * sqrt(w)), symmetric = TRUE)
	vals <- pmax(eig$values, 1e-8 * max(eig$values), 1e-10)
	.bay_mvt_make(centre, vals, eig$vectors)
}

.bay_mvt_ld <- function(Z, fit, df, f) {
	d <- fit$d
	z2 <- rowSums((sweep(Z, 2L, fit$centre, "-") %*% fit$root_inv)^2)
	lgamma((df + d) / 2) - lgamma(df / 2) - d / 2 * log(df * pi) - d * log(f) -
		fit$log_det_root - (df + d) / 2 * log1p(z2 / (f^2 * df))
}

.bay_mvt_draw <- function(S, fit, df, f) {
	Z <- matrix(stats::rnorm(S * fit$d), S, fit$d)
	scale <- f * sqrt(df / stats::rchisq(S, df))
	sweep((Z %*% t(fit$root)) * scale, 2L, fit$centre, "+")
}

# Mixture of such components: list(list(fit, df, f, weight), ...).
.bay_mix_draw <- function(S, comps) {
	ns <- as.numeric(stats::rmultinom(1L, S, vapply(comps, `[[`, numeric(1), "weight")))
	do.call(rbind, lapply(seq_along(comps), function(j) {
		if (ns[j] == 0L) return(NULL)
		.bay_mvt_draw(ns[j], comps[[j]]$fit, comps[[j]]$df, comps[[j]]$f)
	}))
}

.bay_mix_ld <- function(Z, comps) {
	A <- vapply(comps, function(cp) {
		log(cp$weight) + .bay_mvt_ld(Z, cp$fit, cp$df, cp$f)
	}, numeric(nrow(Z)))
	if (nrow(Z) == 1L) A <- matrix(A, nrow = 1L)
	.bay_row_lse(A)
}

# Additive log ratios of a point of the simplex, and back. The density on the simplex
# is the density of z divided by prod(th).
.bay_alr <- function(th) log(th[, -ncol(th), drop = FALSE]) - log(th[, ncol(th)])

.bay_alr_inv <- function(z) {
	top <- pmax(0, apply(z, 1L, max))
	den <- exp(-top) + rowSums(exp(z - top))
	th <- cbind(exp(z - top), exp(-top)) / den
	th <- pmax(th, 1e-300)
	th / rowSums(th)
}

# Adaptive importance sampling on the proportions of precision, for one target: a
# function of th giving a matrix (draws x number of r). First a pilot from the
# defensive mixture, then the final draws from that mixture and a t fitted on the
# pilot. Only the final draws enter the estimate. The defensive part keeps the tails
# and any second mode covered. Column ic drives the pilot and the growth.
.bay_law_het_is <- function(st, u, draws, target, ic, max_grow = 2L) {
	k <- st$k
	pr <- .bay_law_proposal(st, u)
	ldir <- function(th) lgamma(k * u) - k * lgamma(u) + (u - 1) * rowSums(log(th))
	S0 <- max(300L, draws %/% 4L)
	th0 <- .bay_law_draw_def(pr, S0)
	lw0 <- target(th0)[, ic] + ldir(th0) - .bay_law_lq_def(pr, th0)
	fit <- .bay_mvt_fit(.bay_alr(th0), lw0)
	mix <- c(0.3, 0.5, 0.2)
	lw <- NULL
	extra <- draws
	grown <- 0L
	repeat {
		ns <- as.numeric(stats::rmultinom(1L, extra, mix))
		z_new <- rbind(if (ns[2L] > 0L) .bay_mvt_draw(ns[2L], fit, 5, 1),
			if (ns[3L] > 0L) .bay_mvt_draw(ns[3L], fit, 3, 2))
		th <- rbind(if (ns[1L] > 0L) .bay_law_draw_def(pr, ns[1L]),
			if (nrow(z_new) > 0L) .bay_alr_inv(z_new))
		z <- .bay_alr(th)
		ljac <- rowSums(log(th))
		lq <- .bay_row_lse(cbind(log(mix[1L]) + .bay_law_lq_def(pr, th),
			log(mix[2L]) + .bay_mvt_ld(z, fit, 5, 1) - ljac,
			log(mix[3L]) + .bay_mvt_ld(z, fit, 3, 2) - ljac))
		lw <- rbind(lw, target(th) + ldir(th) - lq)
		ess <- .bay_law_is(lw[, ic])$ess
		if (isTRUE(ess >= min(500, draws / 4)) || grown >= max_grow || nrow(lw) >= .bay_law_max_draws()) break
		extra <- min(3L * nrow(lw), .bay_law_max_draws() - nrow(lw))
		grown <- grown + 1L
	}
	lapply(seq_len(ncol(lw)), function(i) .bay_law_is(lw[, i]))
}

# The four-model table of the normal law. rs: values of r; H1 gives one value per r.
.bay_law_normal <- function(st, rs, u, draws, seed, max_grow = 2L) {
	k <- st$k
	const <- .bay_law_const(st)
	one <- matrix(1, 1L, k)
	h0_eq <- const + .bay_law_kernel(st, one, 0)
	h1_eq_zs <- vapply(rs, function(r) {
		h0_eq + .bay_bf_oneway_zs(st$y_raw, st$g_raw, r)
	}, numeric(1))
	h1_eq <- h1_eq_zs
	if (k > 2L) {
		# Same grid as for the other laws, on the row c = 1.
		sp <- NULL
		for (s in seq_len(k)) {
			x <- const + drop(.bay_law_gint(st, one, rs * sqrt(2), s = s))
			sp <- if (is.null(sp)) x else .bay_lse2(sp, x)
		}
		h1_eq <- .bay_lse2(h1_eq_zs, sp - log(k)) + log(0.5)
	}
	ic <- ceiling(length(rs) / 2)
	het <- .bay_with_seed(seed, list(
		h0 = .bay_law_het_is(st, u, draws, function(th) {
			matrix(.bay_law_kernel(st, 1 / (k * th), 0), ncol = 1L)
		}, 1L, max_grow),
		h1 = .bay_law_het_is(st, u, draws, function(th) {
			.bay_law_gint_mix(st, 1 / (k * th), rs)
		}, ic, max_grow)
	))
	s0 <- het$h0[[1L]]
	s1 <- het$h1
	floor_ess <- function(x) {
		if (!isTRUE(x$ess >= .bay_law_min_ess())) x$log_mean <- NA_real_
		x
	}
	s0 <- floor_ess(s0)
	s1 <- lapply(s1, floor_ess)
	nr <- length(rs)
	s1_zs <- s1[nr + seq_len(nr)]
	s1 <- s1[seq_len(nr)]
	list(
		eq_h0 = list(log_ml = h0_eq, rel_se = 0, ess = NA_real_),
		eq_h1 = list(log_ml = h1_eq, log_ml_zs = h1_eq_zs, rel_se = rep(0, length(rs)), ess = NA_real_),
		het_h0 = list(log_ml = const + s0$log_mean, rel_se = s0$rel_se, ess = s0$ess),
		het_h1 = list(log_ml = const + vapply(s1, `[[`, numeric(1), "log_mean"),
			log_ml_zs = const + vapply(s1_zs, `[[`, numeric(1), "log_mean"),
			rel_se = vapply(s1, `[[`, numeric(1), "rel_se"),
			ess = min(vapply(s1, `[[`, numeric(1), "ess")))
	)
}


# ---------------------------------------------------------------------------
# Student t laws: importance sampling around the mode
# ---------------------------------------------------------------------------

# Parameters of one row: theta_1..k (H1) or mu (H0), log taubar, the additive log
# ratios of vartheta if the variances differ, log(nu - 2). Returns the log posterior
# kernel without the prior on theta, and the quantities that prior needs, so that
# the same draws serve every r.
.bay_t_parts <- function(X, st, h1, het, u, s = NULL) {
	k <- st$k
	S <- nrow(X)
	p <- if (!h1) 1L else if (!is.null(s)) 2L else k
	th_loc <- X[, seq_len(p), drop = FALSE]
	lt <- X[, p + 1L]
	if (het) {
		z <- X[, p + 1L + seq_len(k - 1L), drop = FALSE]
		den <- 1 + rowSums(exp(z))
		vth <- cbind(exp(z), 1) / den
		lc <- -log(k) - log(vth)
		lp_dir <- lgamma(k * u) - k * lgamma(u) + u * rowSums(log(vth))
	} else {
		lc <- matrix(0, S, k)
		lp_dir <- 0
	}
	lv <- X[, ncol(X)]
	nu <- 2 + exp(lv)
	lp_nu <- -exp(lv) + lv
	norm_t <- lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(nu * pi)
	ll <- numeric(S)
	for (i in seq_len(k)) {
		# log of the squared scale: s^2 = sigma^2 (nu - 2) / nu
		ls2 <- lc[, i] - lt + lv - log(nu)
		loc <- if (!h1) {
			th_loc[, 1L]
		} else if (!is.null(s)) {
			if (i == s) th_loc[, 1L] + th_loc[, 2L] else th_loc[, 1L]
		} else {
			th_loc[, i]
		}
		yi <- st$y[[i]]
		acc <- numeric(S)
		# Chunks: S x n_i doubles at most 5e6 at a time.
		if (S * length(yi) <= 5e6) {
			acc <- rowSums(log1p(outer(loc, yi, "-")^2 * exp(-ls2) / nu))
		} else {
			size <- max(1L, floor(5e6 / length(yi)))
			for (b in split(seq_len(S), ceiling(seq_len(S) / size))) {
				z2 <- outer(loc[b], yi, "-")^2 * exp(-ls2[b])
				acc[b] <- rowSums(log1p(z2 / nu[b]))
			}
		}
		ll <- ll + st$n[i] * (norm_t - 0.5 * ls2) - (nu + 1) / 2 * acc
	}
	list(base = ll + lp_dir + lp_nu, th_loc = th_loc, lc = lc, lt = lt, s = s)
}

# Prior of theta under H1 given sigma, with g and mu integrated out. With
# d_i = sigma_i^2, Q0 = sum (theta_i - thetabar_w)^2 / d_i:
#   log p = log(r / sqrt(2 pi)) + lgamma(k/2) - (k/2) log((Q0 + r^2)/2)
#           - (k-1)/2 log(2 pi) - 1/2 sum log d_i - 1/2 log(sum 1/d_i)
.bay_t_prior_theta <- function(parts, r, k) {
	if (!is.null(parts$s)) {
		# One group: a_s / sigma_s ~ Cauchy(0, r sqrt(2)), mu flat.
		scale <- r * sqrt(2) * exp(0.5 * (parts$lc[, parts$s] - parts$lt))
		return(stats::dcauchy(parts$th_loc[, 2L], 0, scale, log = TRUE))
	}
	ld <- parts$lc - parts$lt
	w <- exp(-ld)
	W <- rowSums(w)
	tb <- rowSums(w * parts$th_loc) / W
	Q0 <- rowSums(w * (parts$th_loc - tb)^2)
	log(r) - 0.5 * log(2 * pi) + lgamma(k / 2) - (k / 2) * log((Q0 + r^2) / 2) -
		(k - 1) / 2 * log(2 * pi) - 0.5 * rowSums(ld) - 0.5 * log(W)
}

.bay_t_logpost <- function(x, st, h1, het, u, r, s = NULL) {
	X <- matrix(x, nrow = 1L)
	parts <- .bay_t_parts(X, st, h1, het, u, s)
	parts$base + if (h1) .bay_t_prior_theta(parts, r, st$k) else 0
}

.bay_t_start <- function(st, h1, het, robust, s = NULL) {
	k <- st$k
	locs <- if (robust) {
		vapply(st$y, stats::median, numeric(1), USE.NAMES = FALSE)
	} else {
		st$m
	}
	spread <- vapply(st$y, function(v) {
		if (length(v) < 2L) return(NA_real_)
		s <- if (robust) stats::mad(v) else stats::sd(v)
		if (s > 0) s^2 else NA_real_
	}, numeric(1), USE.NAMES = FALSE)
	pooled <- stats::median(spread, na.rm = TRUE)
	if (!is.finite(pooled)) pooled <- 1
	spread[is.na(spread)] <- pooled
	th <- if (!h1) {
		stats::median(unlist(st$y, use.names = FALSE))
	} else if (!is.null(s)) {
		mu0 <- stats::median(unlist(st$y[-s], use.names = FALSE))
		c(mu0, locs[s] - mu0)
	} else {
		locs
	}
	if (het) {
		prec <- 1 / spread
		lt <- log(mean(prec))
		z <- log(prec[-k]) - log(prec[k])
	} else {
		lt <- -log(pooled)
		z <- numeric(0)
	}
	c(th, lt, z, 0)
}

# Mode (two starts), then importance sampling around it: first a t at the mode, then,
# if the pilot shows that this fits badly (a far group, a posterior with a long
# ridge), a t fitted on the weighted pilot, up to three rounds. The t at the mode stays
# in the mixture. Only draws from the final mixture enter the estimate. Any failure of
# the mode search is returned as NA: the caller drops that law and says so.
.bay_t_ml <- function(st, h1, het, u, rs, draws, max_grow = 2L, s = NULL) {
	r0 <- rs[ceiling(length(rs) / 2)]
	ic <- ceiling(length(rs) / 2)
	fn <- function(x) {
		val <- -.bay_t_logpost(x, st, h1, het, u, r0, s)
		if (is.finite(val)) val else 1e10
	}
	best <- NULL
	for (robust in c(TRUE, FALSE)) {
		fit <- tryCatch(stats::optim(.bay_t_start(st, h1, het, robust, s), fn, method = "BFGS",
			control = list(maxit = 300L, reltol = 1e-9)), error = function(e) NULL)
		if (!is.null(fit) && (is.null(best) || fit$value < best$value)) best <- fit
	}
	fail <- list(log_ml = rep(NA_real_, length(rs)), rel_se = rep(NA_real_, length(rs)),
		ess = NA_real_)
	if (is.null(best) || best$value >= 1e10) return(fail)
	H <- tryCatch(stats::optimHess(best$par, fn), error = function(e) NULL)
	if (is.null(H) || anyNA(H)) return(fail)
	eig <- eigen((H + t(H)) / 2, symmetric = TRUE)
	vals <- pmax(eig$values, 1e-3)
	fit_mode <- .bay_mvt_make(best$par, 1 / vals, eig$vectors)
	comps <- list(list(fit = fit_mode, df = 5, f = 1, weight = 0.7),
		list(fit = fit_mode, df = 3, f = 2, weight = 0.3))
	log_weights <- function(X) {
		parts <- .bay_t_parts(X, st, h1, het, u, s)
		pr <- if (h1) {
			vapply(rs, function(r) .bay_t_prior_theta(parts, r, st$k), numeric(nrow(X)))
		} else {
			matrix(0, nrow(X), length(rs))
		}
		if (nrow(X) == 1L) pr <- matrix(pr, nrow = 1L)
		(parts$base + pr) - .bay_mix_ld(X, comps)
	}
	S0 <- max(300L, draws %/% 4L)
	X0 <- .bay_mix_draw(S0, comps)
	lw0 <- log_weights(X0)
	for (round in 1:3) {
		if (isTRUE(.bay_law_is(lw0[, ic])$ess >= S0 / 4)) break
		if (anyNA(lw0[, ic])) break
		fit_new <- .bay_mvt_fit(X0, lw0[, ic])
		comps <- list(list(fit = fit_new, df = 5, f = 1, weight = 0.5),
			list(fit = fit_new, df = 3, f = 2, weight = 0.2),
			list(fit = fit_mode, df = 5, f = 1, weight = 0.3))
		X0 <- .bay_mix_draw(S0, comps)
		lw0 <- log_weights(X0)
	}
	# The pilot chose the mixture: only fresh draws enter the estimate.
	lw <- NULL
	extra <- draws
	grown <- 0L
	repeat {
		lw <- rbind(lw, log_weights(.bay_mix_draw(extra, comps)))
		ess <- .bay_law_is(lw[, ic])$ess
		if (isTRUE(ess >= min(500, draws / 4)) || grown >= max_grow || nrow(lw) >= .bay_law_max_draws()) break
		extra <- min(3L * nrow(lw), .bay_law_max_draws() - nrow(lw))
		grown <- grown + 1L
	}
	sm <- lapply(seq_along(rs), function(i) {
		x <- .bay_law_is(lw[, i])
		if (!isTRUE(x$ess >= .bay_law_min_ess())) x$log_mean <- NA_real_
		x
	})
	list(log_ml = vapply(sm, `[[`, numeric(1), "log_mean"),
		rel_se = vapply(sm, `[[`, numeric(1), "rel_se"),
		ess = min(vapply(sm, `[[`, numeric(1), "ess")))
}


# H1 of a t law: the mixture of the two forms. Above 8 groups, only the 5 groups furthest
# from the weighted median are computed in the second form, for the time it takes. The
# others count for nothing there, which can only lower the Bayes factor, slightly.
.bay_t_h1_mix <- function(st, het, u, rs, draws) {
	zs <- .bay_t_ml(st, TRUE, het, u, rs, draws)
	zs$log_ml_zs <- zs$log_ml
	k <- st$k
	if (k == 2L) return(zs)
	groups <- seq_len(k)
	if (k > 8L) {
		ord <- order(st$m)
		mb <- st$m[ord][which(cumsum(st$n[ord]) / st$N >= 0.5)[1L]]
		groups <- order(-abs(st$m - mb) * sqrt(st$n))[1:5]
	}
	sp <- lapply(groups, function(s) .bay_t_ml(st, TRUE, het, u, rs, draws, s = s))
	lm <- rbind(zs$log_ml + log(0.5), do.call(rbind, lapply(sp, `[[`, "log_ml")) + log(0.5 / k))
	rel <- rbind(zs$rel_se, do.call(rbind, lapply(sp, `[[`, "rel_se")))
	total <- apply(lm, 2L, .bay_log_sum_exp)
	wts <- exp(sweep(lm, 2L, total, "-"))
	list(log_ml = total, log_ml_zs = zs$log_ml_zs,
		rel_se = sqrt(colSums((wts * rel)^2)),
		ess = min(zs$ess, vapply(sp, `[[`, numeric(1), "ess")))
}


# ---------------------------------------------------------------------------
# Eight-model table and averaging
# ---------------------------------------------------------------------------

.bay_law_names <- function() c("normal_eq", "normal_het", "t_eq", "t_het")

# Exact ties in a continuous model. A group with no variation, or mostly one repeated
# value, gives an infinite marginal likelihood under a variance per group as soon as
# (k - 1) u <= (n_i - 1) / 2, and a t likelihood is unbounded when ties reach about
# 2/3 of the data: the integral diverges, it is not a numerical accident (RELECTURE
# 2026-11-07, B2). These laws do not apply there. Returns the laws that do.
.bay_laws_allowed <- function(st, laws) {
	top_share <- function(v) max(table(v)) / length(v)
	big <- st$n >= 2
	by_group <- vapply(st$y[big], top_share, numeric(1))
	if (any(by_group > 1/2)) return(intersect(laws, "normal_eq"))
	if (top_share(unlist(st$y, use.names = FALSE)) > 1/2) {
		return(intersect(laws, c("normal_eq", "normal_het")))
	}
	laws
}

# Long table: one row per model and per r.
.bay_laws_table <- function(st, rs, u, draws, seed, laws = .bay_law_names()) {
	out <- list()
	add <- function(law, hyp, item) {
		out[[length(out) + 1L]] <<- data.frame(law = law, hypothesis = hyp, r = rs,
			log_ml = item$log_ml,
			log_ml_zs = if (is.null(item$log_ml_zs)) item$log_ml else item$log_ml_zs,
			rel_se = item$rel_se, ess = item$ess,
			stringsAsFactors = FALSE)
	}
	if (any(c("normal_eq", "normal_het") %in% laws)) {
		nm <- .bay_law_normal(st, rs, u, draws, seed)
		if ("normal_eq" %in% laws) {
			add("normal_eq", "H0", nm$eq_h0)
			add("normal_eq", "H1", nm$eq_h1)
		}
		if ("normal_het" %in% laws) {
			add("normal_het", "H0", nm$het_h0)
			add("normal_het", "H1", nm$het_h1)
		}
	}
	if (any(c("t_eq", "t_het") %in% laws)) {
		tt <- .bay_with_seed(seed + 1L, {
			lst <- list()
			if ("t_eq" %in% laws) {
				lst$t_eq_h0 <- .bay_t_ml(st, FALSE, FALSE, u, rs, draws)
				lst$t_eq_h1 <- .bay_t_h1_mix(st, FALSE, u, rs, draws)
			}
			if ("t_het" %in% laws) {
				lst$t_het_h0 <- .bay_t_ml(st, FALSE, TRUE, u, rs, draws)
				lst$t_het_h1 <- .bay_t_h1_mix(st, TRUE, u, rs, draws)
			}
			lst
		})
		if ("t_eq" %in% laws) {
			add("t_eq", "H0", tt$t_eq_h0)
			add("t_eq", "H1", tt$t_eq_h1)
		}
		if ("t_het" %in% laws) {
			add("t_het", "H0", tt$t_het_h0)
			add("t_het", "H1", tt$t_het_h1)
		}
	}
	do.call(rbind, out)
}

# Answer at one r from the table: P(H1 | D), inclusion BF across laws, P(law | D).
# prior_law: named weights on the laws. A model whose marginal likelihood failed
# (NA) is dropped with its whole law, so that the comparison stays balanced.
.bay_laws_average <- function(tab, r, prior_h1 = 1/2, prior_law = NULL, h1_form = c("mix", "zs")) {
	tab <- tab[tab$r == r, , drop = FALSE]
	if (match.arg(h1_form) == "zs") tab$log_ml <- tab$log_ml_zs
	laws <- unique(tab$law)
	bad <- unique(tab$law[is.na(tab$log_ml)])
	laws_ok <- setdiff(laws, bad)
	if (length(laws_ok) == 0L) {
		.exit("All sampling laws failed.",
			"Toutes les lois d'\u00e9chantillonnage ont \u00e9chou\u00e9.")
	}
	if (is.null(prior_law)) prior_law <- stats::setNames(rep(1, length(laws)), laws)
	pl <- prior_law[laws_ok]
	pl <- pl / sum(pl)
	tab <- tab[tab$law %in% laws_ok, , drop = FALSE]
	lw <- tab$log_ml + log(pl[tab$law])
	is1 <- tab$hypothesis == "H1"
	l1 <- .bay_log_sum_exp(lw[is1])
	l0 <- .bay_log_sum_exp(lw[!is1])
	log_bf <- l1 - l0
	# Delta method on the log BF: the error of each marginal, weighted.
	w1 <- exp(lw[is1] - l1)
	w0 <- exp(lw[!is1] - l0)
	err <- sqrt(sum((w1 * tab$rel_se[is1])^2) + sum((w0 * tab$rel_se[!is1])^2))
	odds <- prior_h1 / (1 - prior_h1)
	post <- lw + ifelse(is1, log(prior_h1), log(1 - prior_h1))
	post <- exp(post - .bay_log_sum_exp(post))
	p_law <- tapply(post, tab$law, sum)[laws_ok]
	bf_law <- vapply(laws_ok, function(l) {
		tab$log_ml[tab$law == l & is1] - tab$log_ml[tab$law == l & !is1]
	}, numeric(1))
	list(log_bf = unname(log_bf),
		p_h1 = .bay_post_prob(log_bf = log_bf, prior_odds = odds),
		p_h0 = .bay_post_prob(log_bf = -log_bf, prior_odds = 1 / odds),
		error_pct = 100 * err, p_law = p_law, log_bf_law = bf_law,
		laws_dropped = bad,
		models = data.frame(law = tab$law, hypothesis = tab$hypothesis,
			log_ml = tab$log_ml, posterior = post, rel_se = tab$rel_se, ess = tab$ess,
			stringsAsFactors = FALSE))
}
