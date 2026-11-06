# sys_bay_oneway.R
# Internal one-way Bayes factors for the future m.bay(). Nothing is exported.
#   .bay_bf_oneway_zs()     : common variance, Zellner-Siow prior, deterministic.
#   .bay_bf_oneway_hetero() : one variance per group, importance sampling for
#                             H1, exact one-dimensional integral for H0.
# No p-value is computed anywhere in this file.
#
# References (comments only, never in messages: dev/bp.log, 7.4 and BP-008).
#   Rouder, J. N., Morey, R. D., Speckman, P. L., & Province, J. M. (2012).
#     Default Bayes factors for ANOVA designs. Journal of Mathematical
#     Psychology, 56(5), 356-374. https://doi.org/10.1016/j.jmp.2012.08.001
#   Liang, F., Paulo, R., Molina, G., Clyde, M. A., & Berger, J. O. (2008).
#     Mixtures of g priors for Bayesian variable selection. JASA, 103(481),
#     410-423. (a verifier: form of the g-prior integrand)
#   Trefethen, L. N., & Weideman, J. A. C. (2014). The exponentially convergent
#     trapezoidal rule. SIAM Review, 56(3), 385-458. (a verifier: h choice)
#   Golub, G. H., & Welsch, J. H. (1969). Calculation of Gauss quadrature rules.
#     Mathematics of Computation, 23, 221-230. (a verifier)
#   Kong, A. (1992). A note on importance sampling using standardized weights.
#     Technical report, University of Chicago. (a verifier: ESS formula)


# ---------------------------------------------------------------------------
# Shared preparation
# ---------------------------------------------------------------------------

# Filters NA pairs with an announcement, checks the input and returns group
# statistics after an internal change of origin and of scale. Both Bayes factors
# are invariant under these changes: this only protects against overflow and
# loss of precision. A within-group sum of squares under the rounding noise of
# its own data is set to exactly 0, so that "variance zero" has one meaning.
.bay_oneway_prep <- function(y, g) {
	if (!is.numeric(y)) {
		.exit("y must be numeric.", "y doit \u00eatre num\u00e9rique.")
	}
	if (length(y) != length(g)) {
		.exit("y and g must have the same length.",
			"y et g doivent avoir la m\u00eame longueur.")
	}
	y <- as.numeric(y)
	g <- as.character(g)
	missing_pair <- is.na(y) | is.na(g)
	if (any(missing_pair)) {
		warning(.msg(paste0(sum(missing_pair), " observation(s) with NA in y or g removed."),
				paste0(sum(missing_pair), " observation(s) avec NA dans y ou g retir\u00e9e(s).")),
			call. = FALSE)
		y <- y[!missing_pair]
		g <- g[!missing_pair]
	}
	if (any(is.infinite(y))) {
		.exit("y contains infinite values.",
			"y contient des valeurs infinies.")
	}
	g <- factor(g)
	k <- nlevels(g)
	if (k < 2L) {
		.exit("At least 2 groups are required (k < 2).",
			"Au moins 2 groupes sont requis (k < 2).")
	}
	scale0 <- max(abs(y))
	if (scale0 > 0) y <- y / scale0
	by_group <- split(y, g)
	n <- as.numeric(lengths(by_group))
	m <- vapply(by_group, mean, numeric(1), USE.NAMES = FALSE)
	ss <- vapply(by_group, function(v) sum((v - mean(v))^2), numeric(1), USE.NAMES = FALSE)
	noise <- vapply(by_group, function(v) {
		length(v) * (.Machine$double.eps * max(abs(v)))^2
	}, numeric(1), USE.NAMES = FALSE)
	ss[ss <= noise] <- 0
	m <- m - sum(n * m) / sum(n)
	list(k = k, n = n, m = m, ss = ss, N = sum(n), levels = levels(g))
}

.bay_check_rscale <- function(r) {
	if (!is.numeric(r) || length(r) != 1L || !isTRUE(is.finite(r)) || !isTRUE(r > 0)) {
		.exit("r must be one finite, strictly positive number.",
			"r doit \u00eatre un seul nombre fini et strictement positif.")
	}
	as.numeric(r)
}


# ---------------------------------------------------------------------------
# Common variance: Zellner-Siow, deterministic
# ---------------------------------------------------------------------------

# Model (Rouder et al., 2012): y = mu + X a + e, with X = Z P, P an orthonormal
# k x (k-1) contrast basis, a | sigma2, g ~ N(0, sigma2 g I_{k-1}),
# g ~ InvGamma(1/2, r^2/2), p(mu, sigma2) prop 1/sigma2.
#
# With A = X_c'X_c = U diag(lambda) U' (X_c centred), c = U'X_c'y and
# s_j = c_j^2 / lambda_j (sum_j s_j = SSB), one has
#   log BF(g) = -1/2 sum_j log(1 + g lambda_j) - (N-1)/2 log(q(g) / SST)
#   q(g)      = SSE + sum_j s_j / (1 + g lambda_j)
# The usual form SST - c'(A + I/g)^{-1} c subtracts two almost equal numbers when
# F is huge; this one only adds positive numbers, so it keeps full precision.
#
# Integration over t = log g, trapezoids on the log scale. The integrand is
# analytic and decays doubly exponentially on the left and at rate k/2 on the
# right, so the rule converges exponentially (Trefethen & Weideman, 2014, a
# verifier). The bounds are set from the data: the peak lies near
# g lambda_min ~ (SSB/SSE) N. The step is halved until a two-step check passes.
#
# Returns one number. SSE = 0 with SSB > 0: Inf, with a warning and
# attr(., "signal") = "zero_residual". y constant: 0/0, refused rather than NA.
# n = k (no residual df) is refused too: the "perfect fit" there is not evidence.
.bay_bf_oneway_zs <- function(y, g, r = 1/2) {
	r <- .bay_check_rscale(r)
	st <- .bay_oneway_prep(y, g)
	k <- st$k
	n <- st$n
	N <- st$N
	if (N <= k) {
		.exit("No residual degree of freedom (N <= k).",
			"Aucun degr\u00e9 de libert\u00e9 r\u00e9siduel (N <= k).")
	}
	sse <- sum(st$ss)
	ssb <- sum(n * st$m^2)
	if (isTRUE(sse + ssb == 0)) {
		.exit("y is constant: the Bayes factor is undetermined.",
			"y est constant : le facteur de Bayes est ind\u00e9termin\u00e9.")
	}
	if (isTRUE(sse == 0)) {
		warning(.msg("Residual variance is zero: log BF set to Inf.",
				"Variance r\u00e9siduelle nulle : log BF fix\u00e9 \u00e0 Inf."),
			call. = FALSE)
		return(structure(Inf, signal = "zero_residual"))
	}

	# Orthonormal contrasts (normalised Helmert); the prior does not depend on the
	# choice because P P' = I - J/k.
	P <- stats::contr.helmert(k)
	P <- sweep(P, 2L, sqrt(colSums(P^2)), "/")
	A <- crossprod(P, (diag(n, k) - tcrossprod(n) / N) %*% P)
	eig <- eigen(A, symmetric = TRUE)
	lambda <- pmax(eig$values, .Machine$double.eps * max(eig$values))
	cc <- drop(crossprod(eig$vectors, crossprod(P, n * st$m)))
	s <- cc^2 / lambda
	sst <- sse + sum(s)

	log_prior0 <- 0.5 * log(r^2 / 2) - 0.5 * log(pi)
	log_f <- function(t) {
		gl <- outer(exp(t), lambda)
		q <- sse + rowSums(sweep(1 / (1 + gl), 2L, s, "*"))
		log_prior0 - 0.5 * t - r^2 / (2 * exp(t)) -
			0.5 * rowSums(log1p(gl)) - (N - 1) / 2 * (log(q) - log(sst))
	}
	# Left: the prior factor is below exp(-700). Right: rate k/2 >= 1, so + 60 units
	# past the peak leaves under exp(-60), whatever the size of F.
	t_lo <- log(r^2 / 1400)
	t_hi <- max(log(r^2), log(1 + sum(s) / sse) + log(N) - log(min(lambda))) + 60
	h <- min(0.05, 0.25 / sqrt(k))
	repeat {
		t_grid <- seq(t_lo, t_hi, by = h)
		lf <- log_f(t_grid)
		fine <- .bay_log_sum_exp(lf) + log(h)
		coarse <- .bay_log_sum_exp(lf[seq(1L, length(lf), by = 2L)]) + log(2 * h)
		if (isTRUE(abs(fine - coarse) < 1e-9) || h < 1e-3) break
		h <- h / 2
	}
	if (!is.finite(fine)) {
		# Not reachable for valid input (sse > 0): kept so that nothing is ever NA
		# without a message.
		.exit("Numerical failure of the Zellner-Siow integral.",
			"\u00c9chec num\u00e9rique de l'int\u00e9grale de Zellner-Siow.")
	}
	fine
}


# ---------------------------------------------------------------------------
# Group variances: quantities shared by H0 and H1
# ---------------------------------------------------------------------------

# Prior, common to H0 and H1: p(mu) prop 1, p(sigma_i^2) prop 1/sigma_i^2,
# independent. TO EXAMINE BY THE REVIEWER: these are improper priors on
# parameters shared by both models. The arbitrary constants cancel in the ratio
# only because mu and every sigma_i^2 are the same parameters, with the same
# measure, in H0 and H1. H1 adds a proper prior on a given sigma and g, in which
# sigma enters only through tau^2. The usual justification (Jeffreys; Berger,
# Pericchi & Varshavsky, 1998, a verifier) is written for a common scale, not for
# k independent scales: its extension here is a choice, not a theorem.
#
# Group i, with nu_i = n_i - 1 and SS_i = within-group sum of squares:
#   (1/s2) (2 pi s2)^{-n/2} exp(-(SS + n (m - th)^2) / (2 s2))
#     = (2 pi)^{-nu/2} n^{-1/2} (s2)^{-nu/2 - 1} exp(-SS / (2 s2)) N(m; th, s2/n)
# and the proposal q(s2): s2 = SS / chi2(nu), i.e. InvGamma(nu/2, SS/2), has
#   q(s2) = (SS/2)^{nu/2} / Gamma(nu/2) (s2)^{-nu/2 - 1} exp(-SS / (2 s2)).
# So the integrand is C_i q(s2) N(m; th, s2/n), with
#   log C_i = -nu/2 log(2 pi) - 1/2 log(n) + lgamma(nu/2) - nu/2 log(SS/2)
# Note: 2 pi and SS/2 at power nu/2, not n/2; the missing half power of 2 pi
# belongs to the normal density of m.
.bay_hetero_log_const <- function(st) {
	nu <- st$n - 1
	sum(-nu / 2 * log(2 * pi) - 0.5 * log(st$n) + lgamma(nu / 2) -
		nu / 2 * log(st$ss / 2))
}

# Log kernel of the means once mu is integrated out, flat: k-1 dimensional normal
# integral,
#   (2 pi)^{-(k-1)/2} |V|^{-1/2} (1'V^{-1}1)^{-1/2} exp(-Q/2)
# V: draws x k matrix of variances of the means. M: same shape, repeated means.
# Q is written as a sum of positive terms, never as a difference.
.bay_hetero_log_kernel <- function(V, M) {
	w <- 1 / V
	W <- rowSums(w)
	mbar <- rowSums(w * M) / W
	-0.5 * (rowSums(log(V)) + log(W) + rowSums(w * (M - mbar)^2) +
		(ncol(V) - 1) * log(2 * pi))
}

# Gauss-Legendre nodes and weights on [-1, 1] (Golub & Welsch, 1969, a verifier).
.bay_gauss_legendre <- function(n) {
	i <- seq_len(n - 1L)
	b <- i / sqrt(4 * i^2 - 1)
	J <- matrix(0, n, n)
	J[cbind(i, i + 1L)] <- b
	J[cbind(i + 1L, i)] <- b
	e <- eigen(J, symmetric = TRUE)
	ord <- order(e$values)
	list(x = e$values[ord], w = 2 * e$vectors[1L, ord]^2)
}

# log of the exact marginal likelihood under H0. Integrating out the variances:
#   prod_i (2 pi)^{-n_i/2} Gamma(n_i/2) (2 / (SS_i + n_i (m_i - mu)^2))^{n_i/2}
# then a one-dimensional integral over mu. The product of t kernels can have one
# peak per group, so: knots at each group mean, geometric on each side and at the
# scale of that group, then adaptive bisection with Gauss-Legendre, comparing an
# interval to its two halves. The tail is |mu|^{-N} with N >= 4: at 1e6 times the
# span of the data the neglected mass is under 1e-18 of the total.
.bay_hetero_h0_exact <- function(st, tol = 1e-10) {
	k <- st$k
	n <- st$n
	m <- st$m
	ss <- st$ss
	const <- sum(-n / 2 * log(2 * pi) + lgamma(n / 2))
	log_f <- function(mu) {
		res <- numeric(length(mu))
		for (i in seq_len(k)) {
			res <- res + n[i] / 2 * (log(2) - log(ss[i] + n[i] * (mu - m[i])^2))
		}
		res
	}
	s <- sqrt(ss / (n * pmax(n - 1, 1)))
	span <- diff(range(m)) + max(s)
	lower <- min(m) - 1e6 * span
	upper <- max(m) + 1e6 * span
	knots <- c(lower, upper)
	for (i in seq_len(k)) {
		steps <- s[i] * 0.5 * 2^(0:ceiling(log2(2e6 * span / s[i])))
		knots <- c(knots, m[i], m[i] - steps, m[i] + steps)
	}
	knots <- sort(unique(knots[knots >= lower & knots <= upper]))
	a <- knots[-length(knots)]
	b <- knots[-1L]
	gl <- .bay_gauss_legendre(16L)
	log_int <- function(a, b) {
		half <- (b - a) / 2
		mid <- (a + b) / 2
		l <- matrix(log_f(rep(mid, each = 16L) + rep(half, each = 16L) * gl$x),
			nrow = 16L)
		top <- apply(l, 2L, max)
		top + log(colSums(gl$w * exp(sweep(l, 2L, top, "-")))) + log(half)
	}
	done <- numeric(0)
	for (iter in seq_len(60L)) {
		mid <- (a + b) / 2
		whole <- log_int(a, b)
		left <- log_int(a, mid)
		right <- log_int(mid, b)
		total <- .bay_log_sum_exp(c(done, left, right))
		err <- abs(exp(whole - total) - exp(left - total) - exp(right - total))
		keep <- err <= tol / max(length(a), 1L)
		done <- c(done, left[keep], right[keep])
		if (all(keep)) break
		a_new <- c(a[!keep], mid[!keep])
		b <- c(mid[!keep], b[!keep])
		a <- a_new
	}
	if (!all(keep)) {
		done <- c(done, left[!keep], right[!keep])
		warning(.msg("H0 integral: bisection limit reached.",
				"Int\u00e9grale de H0 : limite de bissection atteinte."),
			call. = FALSE)
	}
	const + .bay_log_sum_exp(done)
}


# ---------------------------------------------------------------------------
# Group variances: importance sampling
# ---------------------------------------------------------------------------

.bay_hetero_draws <- function(st, n_draws) {
	s2 <- matrix(0, n_draws, st$k)
	for (i in seq_len(st$k)) {
		s2[, i] <- st$ss[i] / stats::rchisq(n_draws, st$n[i] - 1)
	}
	s2
}

# log of the normal part of H1 for each draw, g integrated over t = log g:
#   integral N-kernel(V(g)) IG(g; 1/2, r^2/2) g dt, V_i = s2_i/n_i + g tau^2.
# Replacing g tau^2 (I - J/k) by g tau^2 I: write a' ~ N(0, g tau^2 I_k) as
# c 1 + P b, with c the mean of a' and P b ~ N(0, g tau^2 (I - J/k)). Under a flat
# mu, the shift mu -> mu + c leaves the integral unchanged, with no extra factor,
# so H0 (g = 0) and H1 use the same measure on mu.
#
# Bounds: on the left the prior factor is under exp(-60). When g tau^2 >= SSm and
# g >= k, Q <= 1 and every V_i is dominated by g tau^2, so the integrand decays at
# rate k/2 >= 1: 25 further units leave under exp(-25).
.bay_hetero_log_h1_inner <- function(s2, st, r, h) {
	k <- st$k
	S <- nrow(s2)
	v <- sweep(s2, 2L, st$n, "/")
	tau2 <- rowMeans(s2)
	M <- matrix(rep(st$m, each = S), S, k)
	ssm <- sum((st$m - mean(st$m))^2)
	t_lo <- log(r^2 / 120)
	t_hi <- max(log(r^2), log(ssm) - log(min(tau2)), log(k)) + 25
	t_grid <- seq(t_lo, t_hi, by = h)
	log_prior0 <- 0.5 * log(r^2 / 2) - 0.5 * log(pi)
	acc <- rep(-Inf, S)
	for (t in t_grid) {
		g <- exp(t)
		val <- .bay_hetero_log_kernel(v + g * tau2, M) +
			log_prior0 - 0.5 * t - r^2 / (2 * g)
		hi <- pmax(acc, val)
		acc <- hi + log1p(exp(pmin(acc, val) - hi))
	}
	acc + log(h)
}

# Summary of log weights: log of the mean, ESS (Kong, 1992, a verifier), relative
# Monte Carlo standard error of the mean (delta method, iid draws).
.bay_is_summary <- function(lw) {
	S <- length(lw)
	top <- max(lw)
	w <- exp(lw - top)
	mw <- mean(w)
	list(log_mean = top + log(mw), ess = sum(w)^2 / sum(w^2),
		rel_se = stats::sd(w) / (mw * sqrt(S)))
}

# Sets the RNG for this call only: fixed kinds, so the result does not depend on
# RNGkind() of the user; with_preserve_seed restores seed and kinds on exit.
.bay_with_seed <- function(seed, expr) {
	withr::with_preserve_seed({
		set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
			sample.kind = "Rejection")
		expr
	})
}

# Control of the constants: H0 by the same importance sampling as H1. Compare to
# .bay_hetero_h0_exact(). The proposal ignores the means, so the ESS falls when
# F is large; this estimator is only meant for moderate F.
.bay_hetero_h0_is <- function(st, n_draws = 20000L, seed = 1L) {
	res <- .bay_with_seed(seed, {
		s2 <- .bay_hetero_draws(st, n_draws)
		M <- matrix(rep(st$m, each = n_draws), n_draws, st$k)
		.bay_is_summary(.bay_hetero_log_kernel(sweep(s2, 2L, st$n, "/"), M))
	})
	res$log_ml <- .bay_hetero_log_const(st) + res$log_mean
	res
}

# Step of the g grid: halved until two steps agree on a small pilot.
.bay_hetero_step <- function(s2, st, r) {
	h <- 0.25
	while (h > 1 / 64) {
		a <- .bay_hetero_log_h1_inner(s2, st, r, h)
		b <- .bay_hetero_log_h1_inner(s2, st, r, h / 2)
		if (isTRUE(max(abs(a - b)) < 1e-6)) break
		h <- h / 2
	}
	h
}

# Core, on the statistics. No standardisation here: the tests call it with
# rescaled statistics to check the invariance for real.
#
# ESS rule: min_ess = min(1000, n_draws / 4). ESS = 1000 means a relative error
# near 3 % on the BF, i.e. 0.03 on the log BF. Under this threshold the number of
# draws is multiplied by 4, at most twice; then low_ess = TRUE with a warning.
.bay_hetero_core <- function(st, r, n_draws, seed, max_grow = 2L,
                             min_ess = min(1000, n_draws / 4)) {
	log_h0 <- .bay_hetero_h0_exact(st)
	log_c <- .bay_hetero_log_const(st)
	res <- .bay_with_seed(seed, {
		h <- .bay_hetero_step(.bay_hetero_draws(st, 200L), st, r)
		lw <- numeric(0)
		extra <- n_draws
		grown <- 0L
		repeat {
			lw <- c(lw, .bay_hetero_log_h1_inner(.bay_hetero_draws(st, extra), st, r, h))
			sm <- .bay_is_summary(lw)
			if (isTRUE(sm$ess >= min_ess) || grown >= max_grow) break
			extra <- 3L * length(lw)
			grown <- grown + 1L
		}
		list(sm = sm, n = length(lw), h = h, grown = grown)
	})
	sm <- res$sm
	log_h1 <- log_c + sm$log_mean
	low_ess <- !isTRUE(sm$ess >= min_ess)
	if (low_ess) {
		warning(.msg(paste0("Low effective sample size (", round(sm$ess),
				"): the Bayes factor is imprecise."),
				paste0("Taille d'\u00e9chantillon effective faible (", round(sm$ess),
				") : le facteur de Bayes est impr\u00e9cis.")),
			call. = FALSE)
	}
	list(log_bf = log_h1 - log_h0, error_pct = 100 * sm$rel_se, ess = sm$ess,
		n_draws = res$n, rscale = r, tau_def = "mean_var", method = "hetero_is",
		log_ml_h1 = log_h1, log_ml_h0 = log_h0, low_ess = low_ess,
		grid_step = res$h, n_grow = res$grown)
}

# Checks specific to the variance-per-group model, then scale so that the mean
# within-group variance is 1 (invariance: see the tests).
.bay_hetero_prep <- function(y, g) {
	st <- .bay_oneway_prep(y, g)
	if (any(st$n < 2)) {
		.exit(paste0("Group(s) with n_i < 2, variance not estimable: ",
				paste(st$levels[st$n < 2], collapse = ", "), "."),
			paste0("Groupe(s) avec n_i < 2, variance non estimable : ",
				paste(st$levels[st$n < 2], collapse = ", "), "."))
	}
	if (any(st$ss == 0)) {
		.exit(paste0("Group(s) with zero variance (SS_i = 0): ",
				paste(st$levels[st$ss == 0], collapse = ", "), "."),
			paste0("Groupe(s) de variance nulle (SS_i = 0) : ",
				paste(st$levels[st$ss == 0], collapse = ", "), "."))
	}
	unit <- sqrt(mean(st$ss / (st$n - 1)))
	st$m <- st$m / unit
	st$ss <- st$ss / unit^2
	st
}

# Model: y_ij ~ N(mu + a_i, sigma_i^2).
# H1: a | sigma, g ~ N(0, g tau^2 (I - J/k)), tau^2 = unweighted mean of the
# sigma_i^2, g ~ InvGamma(1/2, r^2/2). H0: a = 0. Priors on mu and on the
# sigma_i^2: see .bay_hetero_log_const. With k = 2 and r = 1/2,
# delta = (a1 - a2)/tau ~ Cauchy(0, sqrt(2)/2).
# Why the proposal is safe: the normal kernel is bounded in sigma, so the weights
# are bounded and their variance is finite.
#
# error_pct: Monte Carlo error on the BF, in percent; H0 is deterministic, so it
# comes from H1 alone. It does not include any error of the model or of the grid.
# STATUS: on trial, must not be used by m.bay() as it stands. The computation has
# been checked, the prior has not held up: a small distant group or an unrelated
# dispersed group pushes the BF towards H0 (dev/mbay/RELECTURE_ADVERSE_2026-11-06.md,
# B1 and B2; reproduced by dev/mbay/scripts/agent_30_verif_coord_defauts.R).
.bay_bf_oneway_hetero <- function(y, g, r = 1/2, n_draws = 20000L, seed = 1L) {
	r <- .bay_check_rscale(r)
	if (!is.numeric(n_draws) || length(n_draws) != 1L || !isTRUE(n_draws >= 100)) {
		.exit("n_draws must be one number >= 100.",
			"n_draws doit \u00eatre un seul nombre >= 100.")
	}
	if (!is.numeric(seed) || length(seed) != 1L ||
		!isTRUE(abs(seed) <= .Machine$integer.max)) {
		.exit("seed must be one integer.",
			"seed doit \u00eatre un seul nombre entier.")
	}
	.bay_hetero_core(.bay_hetero_prep(y, g), r, as.integer(n_draws), as.integer(seed))
}
