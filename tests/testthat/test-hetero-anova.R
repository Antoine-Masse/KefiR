# Internal heteroscedastic ANOVA (Brunner, Dette & Munk, 1997), which replaced
# fda.usc::fanova.hetero() in m.test. It follows the article (unbiased group variances),
# not fda.usc, which shrinks each variance by (n - 1) / n.

# The formulas of the article written with full matrices, as in the GFD package
# (calculate_statistics.R, ATS), to check the simplified code against them.
bdm_full_matrices <- function(x, g) {
	g <- factor(g)
	n <- as.numeric(table(g))
	k <- length(n)
	N <- sum(n)
	means <- as.numeric(tapply(x, g, mean))
	Sn <- N * diag(as.numeric(tapply(x, g, var)) / n)
	H <- diag(k) - matrix(1 / k, k, k)
	C <- t(H) %*% MASS::ginv(H %*% t(H)) %*% H
	D <- diag(diag(C))
	spur <- sum(diag(D %*% Sn))
	Lambda <- diag(1 / (n - 1))
	ats <- N / spur * drop(t(means) %*% C %*% means)
	df1 <- spur^2 / sum(diag(C %*% Sn %*% C %*% Sn))
	df2 <- spur^2 / sum(diag(D %*% D %*% Sn %*% Sn %*% Lambda))
	c(ats, df1, df2, pf(ats, df1, df2, lower.tail = FALSE))
}

test_that(".bdm_oneway equals Welch's t test with two groups", {
	set.seed(1)
	x <- c(rnorm(7, 0, 1), rnorm(19, 0.8, 4))
	g <- rep(c("a", "b"), c(7, 19))
	res <- KefiR:::.bdm_oneway(x, g)
	welch <- t.test(x ~ g, var.equal = FALSE)
	expect_equal(unname(res["statistic"]), unname(welch$statistic)^2)
	expect_equal(unname(res["df1"]), 1)
	expect_equal(unname(res["df2"]), unname(welch$parameter))
	expect_equal(unname(res["p.value"]), welch$p.value)
})

test_that(".bdm_oneway equals the full-matrix formulas of the article", {
	set.seed(2)
	for (i in 1:20) {
		k <- sample(3:7, 1)
		n <- sample(2:30, k, replace = TRUE)
		g <- rep(paste0("g", seq_len(k)), n)
		x <- rnorm(sum(n), rep(rnorm(k), n), rep(runif(k, 0.2, 5), n))
		expect_equal(unname(KefiR:::.bdm_oneway(x, g)), bdm_full_matrices(x, g))
	}
})

test_that(".bdm_oneway uses unbiased variances, unlike fda.usc", {
	# Equal group sizes: fda.usc's statistic is exactly n / (n - 1) times this one.
	# Reference: fda.usc 2.2.0, fanova.hetero(data.frame(x, g), x ~ g)$ans[1, 1] = 2.3877551
	x <- c(1, 2, 4, 3, 5, 9, 2, 8, 11)
	g <- rep(c("a", "b", "c"), each = 3)
	res <- KefiR:::.bdm_oneway(x, g)
	expect_equal(unname(res["statistic"]) * 3 / 2, 2.3877551, tolerance = 1e-7)
})

test_that(".bdm_oneway ignores missing pairs and the order of the data", {
	x <- c(5, 1, 7, 3, 9, 2, 8, 4, 6, 2)
	g <- c("b", "a", "c", "b", "c", "a", "c", "b", "a", "c")
	ref <- KefiR:::.bdm_oneway(x, g)
	expect_named(ref, c("statistic", "df1", "df2", "p.value"))
	o <- c(10, 3, 7, 1, 9, 2, 8, 4, 6, 5)
	expect_equal(KefiR:::.bdm_oneway(x[o], g[o]), ref)
	expect_equal(KefiR:::.bdm_oneway(c(x, NA, 3), c(g, "a", NA)), ref)
})

test_that(".bdm_oneway returns NULL when it cannot be computed", {
	expect_null(KefiR:::.bdm_oneway(1:5, c("a", "a", "b", "b", "c")))
	expect_null(KefiR:::.bdm_oneway(rep(1, 6), rep(c("a", "b"), 3)))
	expect_null(KefiR:::.bdm_oneway(1:4, rep("a", 4)))
})

test_that("m.test still runs its Welch route without fda.usc", {
	set.seed(42)
	g <- factor(rep(c("a", "b", "c"), each = 25))
	x <- rnorm(75, mean = rep(c(0, 0.3, 1), each = 25), sd = rep(c(1, 3, 6), each = 25))
	out <- capture.output(res <- m.test(x, g, verbose = TRUE, plot = FALSE))
	expect_true(any(grepl("Welch", out)))
})
