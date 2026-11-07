# sys_bay_oneway.R
# Internal one-way Bayes factor with common variance (m.bay). Nothing is exported.
#   .bay_bf_oneway_zs() : common variance, Zellner-Siow prior, deterministic.
#   The variance-per-group trial of 6 Nov. 2026 was removed (defects B1, B2):
#   see R/sys_bay_laws.R for what replaced it.
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
	# Between-group variation under the rounding noise is no variation at all.
	if (isTRUE(ssb <= N * .Machine$double.eps^2)) ssb <- 0
	if (isTRUE(sse + ssb == 0)) {
		.exit("y is constant: the Bayes factor is undetermined.",
			"y est constant : le facteur de Bayes est ind\u00e9termin\u00e9.")
	}
	if (isTRUE(sse == 0)) {
		# N = k + 1: one residual df, the same argument as N = k applies.
		if (isTRUE(N == k + 1L)) {
			.exit("One residual degree of freedom and no residual variance.",
				"Un seul degr\u00e9 de libert\u00e9 r\u00e9siduel et aucune variance r\u00e9siduelle.")
		}
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
# Seed
# ---------------------------------------------------------------------------

# Sets the RNG for this call only: fixed kinds, so the result does not depend on
# RNGkind() of the user; with_preserve_seed restores seed and kinds on exit.
.bay_with_seed <- function(seed, expr) {
	withr::with_preserve_seed({
		set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
			sample.kind = "Rejection")
		expr
	})
}
