#' Multivariate normality tests without external package
#'
#' Mardia's skewness and kurtosis tests (Mardia 1970) and the Henze-Zirkler test
#' (Henze & Zirkler 1990), written from their closed forms. They replace the call to
#' the MVN package, which brought many dependencies for these two tests alone.
#' The values are checked against MVN in tests/testthat/test-mvn-tests.R.
#'
#' @param x Numeric matrix or data.frame, one row per individual, complete cases.
#' @return `.mvn_mardia()`: named vector `skew`, `kurt` of p-values.
#'   `.mvn_hz()`: the p-value. `NULL` when the test cannot be computed
#'   (fewer rows than columns + 2, or singular covariance).
#' @references
#' Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with
#' applications. Biometrika, 57(3), 519-530.
#' Henze, N., & Zirkler, B. (1990). A class of invariant consistent tests for
#' multivariate normality. Communications in Statistics - Theory and Methods,
#' 19(10), 3595-3617.
#' @keywords internal
#' @noRd
.mvn_prepare <- function(x) {
	x <- as.matrix(x)
	storage.mode(x) <- "double"
	x <- x[stats::complete.cases(x), , drop = FALSE]
	n <- nrow(x); p <- ncol(x)
	if (n < p + 2) return(NULL)
	xc <- scale(x, center = TRUE, scale = FALSE)
	# Maximum likelihood covariance (divisor n), as in both original papers
	S <- crossprod(xc) / n
	Sinv <- tryCatch(solve(S), error = function(e) NULL)
	if (is.null(Sinv)) return(NULL)
	list(D = xc %*% Sinv %*% t(xc), n = n, p = p)
}

.mvn_mardia <- function(x) {
	m <- .mvn_prepare(x)
	if (is.null(m)) return(NULL)
	n <- m$n; p <- m$p; D <- m$D
	b1 <- sum(D^3) / n^2
	b2 <- sum(diag(D)^2) / n
	df <- p * (p + 1) * (p + 2) / 6
	# Small-sample correction of the skewness statistic (Mardia 1974), under n = 20
	k <- (p + 1) * (n + 1) * (n + 3) / (n * ((n + 1) * (p + 1) - 6))
	skew <- if (n < 20) n * k * b1 / 6 else n * b1 / 6
	kurt <- (b2 - p * (p + 2)) * sqrt(n / (8 * p * (p + 2)))
	c(skew = stats::pchisq(skew, df, lower.tail = FALSE),
	  kurt = 2 * stats::pnorm(abs(kurt), lower.tail = FALSE))
}

.mvn_hz <- function(x) {
	m <- .mvn_prepare(x)
	if (is.null(m)) return(NULL)
	n <- m$n; p <- m$p; D <- m$D
	dj <- diag(D)
	djk <- outer(dj, dj, "+") - 2 * D
	b <- 1 / sqrt(2) * ((2 * p + 1) / 4)^(1 / (p + 4)) * n^(1 / (p + 4))
	hz <- n * (sum(exp(-b^2 / 2 * djk)) / n^2 -
		2 * (1 + b^2)^(-p / 2) * sum(exp(-b^2 / (2 * (1 + b^2)) * dj)) / n +
		(1 + 2 * b^2)^(-p / 2))
	wb <- (1 + b^2) * (1 + 3 * b^2)
	a <- 1 + 2 * b^2
	mu <- 1 - a^(-p / 2) * (1 + p * b^2 / a + p * (p + 2) * b^4 / (2 * a^2))
	si2 <- 2 * (1 + 4 * b^2)^(-p / 2) +
		2 * a^(-p) * (1 + 2 * p * b^4 / a^2 + 3 * p * (p + 2) * b^8 / (4 * a^4)) -
		4 * wb^(-p / 2) * (1 + 3 * p * b^4 / (2 * wb) + p * (p + 2) * b^8 / (2 * wb^2))
	# Log-normal approximation of the null distribution (Henze & Zirkler 1990)
	pmu <- log(sqrt(mu^4 / (si2 + mu^2)))
	psi <- sqrt(log((si2 + mu^2) / mu^2))
	stats::plnorm(hz, pmu, psi, lower.tail = FALSE)
}
