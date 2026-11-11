# Several responses, one factor (D22): stacked standardised responses, lmBF.

.manova_data <- function(seed, eff = c(0.8, 0, 0), n = 15L) {
	set.seed(seed)
	g <- factor(rep(c("a", "b", "c"), each = n))
	R <- matrix(0.3, 3, 3); diag(R) <- 1
	Y <- matrix(stats::rnorm(3L * n * 3L), ncol = 3) %*% chol(R)
	Y[g == "c", ] <- sweep(Y[g == "c", , drop = FALSE], 2L, eff, "+")
	data.frame(g = g, y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3])
}

test_that("MANOVA: same Bayes factors as lmBF on the stacked data", {
	skip_if_not_installed("BayesFactor")
	d <- .manova_data(1)
	r <- m.bay(cbind(y1, y2, y3) ~ g, data = d, verbose = FALSE, sensitivity = FALSE, draws = 2000)
	expect_identical(r$plan, "manova")
	long <- data.frame(y = as.vector(scale(as.matrix(d[2:4]))), g = rep(d$g, 3), vd = factor(rep(c("y1", "y2", "y3"), each = 45)),
		id = factor(rep(1:45, 3)))
	set.seed(2)
	b0 <- BayesFactor::lmBF(y ~ vd + id, data = long, whichRandom = "id", rscaleFixed = 0.5, iterations = 50000, progress = FALSE)
	b2 <- BayesFactor::lmBF(y ~ g + vd + g:vd + id, data = long, whichRandom = "id", rscaleFixed = 0.5, iterations = 50000, progress = FALSE)
	ref <- BayesFactor::extractBF(b2 / b0, logbf = TRUE)$bf
	expect_lt(abs(r$forms$log_bf[2L] - ref), 0.15)
	expect_equal(nrow(r$responses), 3L)
	expect_true(r$effects$p_h1 > 0 && r$effects$p_h1 < 1)
})

test_that("MANOVA: story, bilan and refusals", {
	d <- .manova_data(3, eff = c(1.5, 1, 0))
	txt <- capture.output(r <- m.bay(cbind(y1, y2, y3) ~ g, data = d, draws = 1000))
	expect_true(any(grepl("^1[)] Plan reconnu : 3 r", txt)))
	expect_true(any(grepl("^6[)] ", txt)))
	expect_true(any(grepl("^7[)] Bilan", txt)))
	expect_gt(r$effects$p_h1, 0.9)
	d$h <- factor(rep(1:3, 15))
	expect_error(m.bay(cbind(y1, y2) ~ g * h, data = d, verbose = FALSE), "un seul facteur|one factor")
})
