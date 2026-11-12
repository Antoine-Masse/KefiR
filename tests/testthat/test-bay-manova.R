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
	d$x <- stats::rnorm(45)
	expect_error(m.bay(cbind(y1, y2) ~ g + x, data = d, verbose = FALSE), "covariable|covariate")
})

test_that("MANOVA with two factors (D23): term answers, forms, story", {
	set.seed(1)
	d <- expand.grid(rep = 1:12, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
	d$y1 <- stats::rnorm(nrow(d)) + 1.2 * (d$A == "a2"); d$y2 <- 0.4 * d$y1 + stats::rnorm(nrow(d)); d$y3 <- stats::rnorm(nrow(d))
	txt <- capture.output(r <- m.bay(cbind(y1, y2, y3) ~ A * B, data = d, draws = 1000, sensitivity = FALSE))
	expect_identical(r$plan, "manova")
	expect_identical(r$effects$term, c("A", "B", "A:B"))
	expect_gt(r$effects$p_h1[1L], 0.9)
	expect_lt(r$effects$p_h1[3L], 0.5)
	expect_equal(nrow(r$models), 17L)
	expect_equal(sum(r$models$posterior), 1, tolerance = 1e-8)
	expect_true(any(grepl("a1 x b1", txt, fixed = TRUE)))
})

test_that("several factors with a yes: the bilan leads with the terms, never undecided (audit of 9c44239)", {
	set.seed(11)
	d <- expand.grid(rep = 1:12, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
	d$u <- stats::rnorm(nrow(d)) + 1.3 * (d$A == "a2"); d$v <- stats::rnorm(nrow(d)); d$w <- stats::rnorm(nrow(d))
	txt <- capture.output(r <- m.bay(cbind(u, v, w) ~ A * B, data = d, draws = 1000, sensitivity = FALSE))
	ttl <- txt[grep("^7[)] Bilan", txt)]
	expect_match(ttl, "^7[)] Bilan - A oui")
	expect_false(any(grepl("indécis|undecided", txt)))
	expect_true(any(grepl("u seule", txt, fixed = TRUE)))
	expect_true(all(c("u", "v", "w") %in% r$responses$response))
})
