# sys_bay_fact_t.R
# Heavy tails with several factors between groups (m.bay, DECISIONS, D24). A Student t law
# with a common scale, for every hereditary model, on the measure of the normal laws
# (sys_bay_fact_laws.R): p(mu, sigma) prop 1 / sigma, effects of each term given sigma
# Zellner-Siow, that is beta_t / sigma ~ multivariate Cauchy of scale r once g is integrated
# (Rouder et al. 2012), t of standard deviation sigma with p(nu) = exp(-(nu - 2)), nu > 2, as
# with one factor (sys_bay_laws.R; Maier et al. 2025). The marginal likelihood is computed by
# importance sampling around the mode (multivariate t, df 4). The same code with the normal
# likelihood gives back the lmBF Bayes factors: that is the check of the measure
# (test-bay-fact-t.R). Nothing is exported.

# Design per observation, the columns of every term (orthonormal contrasts of BayesFactor).
.bay_ft_design <- function(factors, keys) {
	cl <- as.character(interaction(factors, drop = TRUE))
	cells <- levels(factor(cl))
	Xc <- .bay_fl_design(factors, keys)
	idx <- match(cl, cells)
	lapply(Xc, function(X) X[idx, , drop = FALSE])
}

# Log multivariate Cauchy density of the rows of B (S x k), scale r.
.bay_ft_lmvc <- function(B, r) {
	k <- ncol(B)
	lgamma((1 + k) / 2) - lgamma(1 / 2) - k / 2 * log(pi) - k * log(r) - (1 + k) / 2 * log1p(rowSums(B^2) / r^2)
}

# Log posterior kernel for S parameter rows: mu, beta (all columns), log sigma, log(nu - 2)
# (t only). Xs: the term matrices of the model.
.bay_ft_logpost <- function(P, y, Xs, r, t_law) {
	S <- nrow(P)
	p <- sum(vapply(Xs, ncol, integer(1)))
	mu <- P[, 1L]
	B <- P[, 1L + seq_len(p), drop = FALSE]
	ls <- P[, p + 2L]
	sig <- exp(ls)
	X <- if (p) do.call(cbind, Xs) else matrix(0, length(y), 0L)
	fit <- if (p) tcrossprod(B, X) else matrix(0, S, length(y))
	R <- -sweep(sweep(fit, 1L, mu, "+"), 2L, y, "-")
	N <- length(y)
	if (t_law) {
		nu <- 2 + exp(P[, p + 3L])
		sc <- sig * sqrt((nu - 2) / nu)
		Z <- R / sc
		ll <- rowSums(stats::dt(Z, df = matrix(nu, S, N), log = TRUE)) - N * log(sc)
		lpn <- -(nu - 2) + P[, p + 3L]
	} else {
		ll <- rowSums(stats::dnorm(R / sig, log = TRUE)) - N * ls
		lpn <- 0
	}
	lpb <- 0
	j <- 0L
	for (X1 in Xs) {
		k <- ncol(X1)
		Bt <- B[, j + seq_len(k), drop = FALSE] / sig
		lpb <- lpb + .bay_ft_lmvc(Bt, r) - k * ls
		j <- j + k
	}
	ll + lpb + lpn
}

# Log marginal likelihood of one model for each r, by importance sampling around the mode.
.bay_ft_ml <- function(y, Xs, rs, t_law, draws) {
	p <- sum(vapply(Xs, ncol, integer(1)))
	d <- p + 2L + as.integer(t_law)
	fit0 <- if (p) stats::lm.fit(cbind(1, do.call(cbind, Xs)), y) else list(coefficients = mean(y), residuals = y - mean(y))
	s0 <- sqrt(sum(fit0$residuals^2) / max(1, length(y) - p - 1))
	z0 <- c(fit0$coefficients, log(s0), if (t_law) log(5))
	z0[!is.finite(z0)] <- 0
	f <- function(z) {
		v <- .bay_ft_logpost(matrix(z, 1L), y, Xs, rs[1L], t_law)
		if (is.finite(v)) v else -1e10
	}
	opt <- stats::optim(z0, f, method = "BFGS", control = list(fnscale = -1, maxit = 500))
	H <- tryCatch(stats::optimHess(opt$par, f), error = function(e) NULL)
	ok <- !is.null(H) && all(is.finite(H))
	if (ok) {
		ev <- eigen(-H, symmetric = TRUE)
		ok <- all(ev$values > 0)
	}
	fit <- if (ok) .bay_mvt_make(opt$par, 1 / ev$values, ev$vectors) else .bay_mvt_make(opt$par, rep(0.1, d), diag(d))
	Z <- .bay_mvt_draw(draws, fit, 4, 1.1)
	lq <- .bay_mvt_ld(Z, fit, 4, 1.1)
	res <- lapply(rs, function(r) .bay_law_is(.bay_ft_logpost(Z, y, Xs, r, t_law) - lq))
	# + log 2: the normal laws of sys_bay_laws.R put 1 / taubar on the precision, which is twice
	# the measure 1 / sigma used here (checked on the null model, test-bay-fact-t.R).
	list(log_ml = log(2) + vapply(res, `[[`, numeric(1), "log_mean"), ess = min(vapply(res, `[[`, numeric(1), "ess")),
		rel_se = vapply(res, `[[`, numeric(1), "rel_se"))
}

# All the hereditary models, t law (t_law = TRUE) or normal (check). y standardised first:
# the ratios do not depend on the scale, the same transformation is applied to every law.
.bay_ft_all <- function(y, factors, sp, rs, draws, seed, t_law = TRUE) {
	# Same standardisation as .bay_laws_stats (sys_bay_fact_laws.R), so that the t law and the
	# two normal laws stand on the same data and the same measure.
	cl <- interaction(factors, drop = TRUE)
	y1 <- y / max(abs(y))
	unit <- sqrt(sum(tapply(y1, cl, function(v) sum((v - mean(v))^2))) / (length(y1) - nlevels(cl)))
	ys <- (y1 - mean(y1)) / unit
	Xall <- .bay_ft_design(factors, sp$keys)
	.bay_with_seed(seed, lapply(sp$models, function(m) .bay_ft_ml(ys, Xall[m], rs, t_law, draws)))
}
