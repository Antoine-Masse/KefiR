# Variances per cell with several factors (D19, sys_bay_fact_laws.R).

test_that("the common-variance kernel gives back lmBF and the one-factor Bayes factor", {
	skip_if_not_installed("BayesFactor")
	set.seed(11)
	d <- expand.grid(A = factor(1:3), B = factor(1:2), rep = 1:8)
	d$y <- stats::rnorm(nrow(d)) + 0.5 * (d$A == "2")
	f <- d[c("A", "B")]
	sp <- .bay_fact_space(c("A", "B", "A:B"))
	st <- .bay_laws_stats(d$y, interaction(f, drop = TRUE))
	Xs <- .bay_fl_design(f, sp$keys)
	eq0 <- .bay_fl_eq_ml(st, list(), 0.5, 1L)
	# A alone: exact one-factor computation.
	set.seed(1)
	expect_lt(abs(.bay_fl_eq_ml(st, Xs[1L], 0.5, 20000L) - eq0 - .bay_bf_oneway_zs(d$y, d$A, 0.5)), 0.05)
	set.seed(2)
	expect_lt(abs(.bay_fl_eq_ml(st, Xs[2L], 0.5, 20000L) - eq0 - .bay_bf_oneway_zs(d$y, d$B, 0.5)), 0.05)
	# A + B: lmBF.
	set.seed(3)
	lb <- BayesFactor::extractBF(BayesFactor::lmBF(y ~ A + B, data = d, rscaleFixed = 0.5, iterations = 50000, progress = FALSE), logbf = TRUE)$bf
	expect_lt(abs(.bay_fl_eq_ml(st, Xs[1:2], 0.5, 20000L) - eq0 - lb), 0.06)
})

test_that("unequal variances per cell are averaged, not capped", {
	set.seed(640001L)
	n <- c(5, 5, 5, 40); sdv <- c(4, 4, 4, 1); cell <- rep(1:4, n)
	d <- data.frame(A = factor(c("a1", "a2", "a1", "a2")[cell]), B = factor(c("b1", "b1", "b2", "b2")[cell]))
	d$y <- stats::rnorm(length(cell), 0, sdv[cell])
	r <- suppressMessages(m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE, draws = 2000))
	expect_true(r$het_used)
	expect_false(r$het_cells)
	expect_gt(r$laws$posterior[r$laws$law == "normal_het"], 0.9)
	expect_true(all(r$effects$p_h1 < 0.95))
	txt <- capture.output(print(r))
	expect_true(any(grepl("moyenn|averaged", txt)))
})

test_that("fewer than 3 values in a cell keeps the common variance alone", {
	set.seed(5)
	d <- expand.grid(A = factor(1:2), B = factor(1:2), rep = 1:2)
	d$y <- stats::rnorm(nrow(d))
	r <- suppressMessages(m.bay(y ~ A + B, data = d, verbose = FALSE, sensitivity = FALSE, draws = 1000))
	expect_false(r$het_used)
	expect_null(r$laws)
})

test_that("several factors: main effect levels, or cells when the interaction is probable", {
	set.seed(21)
	d <- expand.grid(rep = 1:12, A = factor(c("a1", "a2", "a3")), B = factor(c("b1", "b2")))
	d$y <- c(0, 1.5, 1.5)[d$A] + stats::rnorm(nrow(d))
	r <- suppressMessages(m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE, draws = 1000))
	ph <- r$posthoc
	expect_identical(ph$method, "factorial")
	expect_identical(vapply(ph$families, `[[`, character(1), "term"), "A")
	expect_false(ph$families[[1L]]$cells)
	pr <- ph$families[[1L]]$pairs
	expect_gt(pr$p_differ[pr$group1 == "a1" & pr$group2 == "a2"], 0.9)
	expect_lt(pr$p_differ[pr$group1 == "a2" & pr$group2 == "a3"], 0.5)
	expect_true("A:B" %in% ph$skipped$term)
	set.seed(22)
	d2 <- expand.grid(rep = 1:15, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
	d2$y <- 2 * (d2$A == "a2" & d2$B == "b2") + stats::rnorm(nrow(d2))
	r2 <- suppressMessages(m.bay(y ~ A * B, data = d2, verbose = FALSE, sensitivity = FALSE, draws = 1000))
	f <- r2$posthoc$families
	expect_length(f, 1L)
	expect_true(f[[1L]]$cells)
	expect_identical(f[[1L]]$k, 4L)
	txt <- capture.output(print(r2))
	expect_true(any(grepl("cases A x B|cells A x B", txt)))
})
