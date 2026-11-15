# sys_hetero_anova.R
# Internal one-way heteroscedastic ANOVA (m.test). Nothing is exported.
#   .bdm_oneway() : ANOVA-type statistic with Box approximation, variances free per group.
# It replaces fda.usc::fanova.hetero() for the one-factor case, with the same result,
# so that KefiR no longer installs fda.usc and its dependencies.
#
# Reference (comments only):
#   Brunner, E., Dette, H., & Munk, A. (1997). Box-type approximations in
#     nonparametric factorial designs. Journal of the American Statistical
#     Association, 92(440), 1494-1502. https://doi.org/10.1080/01621459.1997.10473671

# Returns c(statistic, df1, df2, p.value), or NULL when it cannot be computed
# (fewer than 2 groups, a group with fewer than 2 values, or no variance at all).
#
# With k groups of sizes n_i, means m_i and variances s2_i, N = sum(n_i),
# M = I - J/k and S = N diag(v_i), v_i = s2_i (n_i - 1) / n_i^2 (as in fda.usc):
#   F   = N m' M m / tr(D S)             D = diag(M) = (1 - 1/k) I
#   df1 = tr(D S)^2 / tr(M S M S)
#   df2 = tr(D S)^2 / tr(D^2 S^2 Delta)  Delta = diag(1 / (n_i - 1))
.bdm_oneway <- function(x, g) {
	keep <- !is.na(x) & !is.na(g)
	x <- as.numeric(x[keep])
	g <- droplevels(as.factor(g[keep]))
	k <- nlevels(g)
	if (k < 2L) return(NULL)
	n <- as.numeric(table(g))
	if (any(n < 2)) return(NULL)
	m <- as.numeric(tapply(x, g, mean))
	v <- as.numeric(tapply(x, g, stats::var)) * (n - 1) / n^2
	N <- sum(n)
	S <- N * v
	tr_ds <- (1 - 1 / k) * sum(S)
	if (!isTRUE(is.finite(tr_ds)) || tr_ds <= 0) return(NULL)
	M <- diag(k) - 1 / k
	MS <- M %*% diag(S, k)
	stat <- N * sum(m * drop(M %*% m)) / tr_ds
	df1 <- tr_ds^2 / sum(diag(MS %*% MS))
	df2 <- tr_ds^2 / ((1 - 1 / k)^2 * sum(S^2 / (n - 1)))
	c(statistic = stat, df1 = df1, df2 = df2,
		p.value = stats::pf(stat, df1, df2, lower.tail = FALSE))
}
