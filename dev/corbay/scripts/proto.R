# Prototype: Bayes factors of cor.bay (Pearson, Kendall, Spearman on ranks).
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor))
# log BF10, log BF+0, log BF-0 for Pearson (correlationBF, stretched beta of width kappa).
lbf_pearson <- function(x, y, kappa = 1/3) {
	b <- correlationBF(y, x, rscale = kappa, nullInterval = c(0, 1))
	l <- extractBF(b, logbf = TRUE)$bf
	c(plus = l[1], minus = l[2], two = log((exp(l[1] - max(l)) + exp(l[2] - max(l))) / 2) + max(l))
}
# Kendall, van Doorn et al. 2018, eq. 6, 9, 11-13: T* normal, mean 3/2 tau sqrt(n), sd 1.
kendall_tstar <- function(x, y) {
	# Same as .corbay_tstar of the package: S standardised with its tie-corrected variance.
	n <- length(x); taub <- suppressWarnings(stats::cor(x, y, method = "kendall"))
	tx <- as.numeric(table(x)); ty <- as.numeric(table(y)); n0 <- n * (n - 1) / 2
	n1 <- sum(tx * (tx - 1)) / 2; n2 <- sum(ty * (ty - 1)) / 2; v0 <- n * (n - 1) * (2 * n + 5) / 18
	v <- (n * (n - 1) * (2 * n + 5) - sum(tx * (tx - 1) * (2 * tx + 5)) - sum(ty * (ty - 1) * (2 * ty + 5))) / 18 +
		sum(tx * (tx - 1)) * sum(ty * (ty - 1)) / (2 * n * (n - 1)) + sum(tx * (tx - 1) * (tx - 2)) * sum(ty * (ty - 1) * (ty - 2)) / (9 * n * (n - 1) * (n - 2))
	k <- sqrt((n0 - n1) * (n0 - n2))
	c(ts = taub * k / sqrt(v), m = 1.5 * sqrt(n) * (k / sqrt(v)) / (n0 / sqrt(v0)))
}
prior_tau <- function(tau, kappa) {
	a <- 1 / kappa
	rho <- sin(pi * tau / 2)
	stats::dbeta((rho + 1) / 2, a, a) / 2 * (pi / 2) * cos(pi * tau / 2)
}
lbf_kendall <- function(x, y, kappa = 1/3) {
	n <- length(x); tm <- kendall_tstar(x, y); ts <- tm[["ts"]]; m <- tm[["m"]]
	if (!is.finite(ts)) return(c(plus = NA, minus = NA, two = NA))
	l0 <- stats::dnorm(ts, 0, 1, log = TRUE)
	f <- function(t) exp(stats::dnorm(ts, m * t, 1, log = TRUE) - l0) * prior_tau(t, kappa)
	pk <- min(max(ts / m, -1), 1); w <- 12 / m
	brk <- sort(unique(c(-1, max(-1, pk - w), 0, min(1, pk + w), 1)))
	seg <- function(a, b) if (b > a) stats::integrate(f, a, b, rel.tol = 1e-8, subdivisions = 500L)$value else 0
	ints <- mapply(seg, brk[-length(brk)], brk[-1])
	pos <- sum(ints[brk[-1] > 0]); neg <- sum(ints[brk[-length(brk)] < 0])
	c(plus = log(2 * pos), minus = log(2 * neg), two = log(pos + neg))
}
lbf_spearman <- function(x, y, kappa = 1/3) lbf_pearson(rank(x), rank(y), kappa)
post3 <- function(l, prior_h1 = 1/2) {
	w <- c(none = log(1 - prior_h1), plus = log(prior_h1 / 2) + l[["plus"]], minus = log(prior_h1 / 2) + l[["minus"]])
	p <- exp(w - max(w)); p / sum(p)
}
