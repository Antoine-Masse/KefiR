# Heavy tails with several factors (D24, sys_bay_fact_t.R).

test_that("the importance sampling with the normal likelihood gives back lmBF and the closed form", {
	skip_if_not_installed("BayesFactor")
	set.seed(3)
	d <- expand.grid(rep = 1:10, A = factor(c("a1", "a2", "a3")), B = factor(c("b1", "b2")))
	d$y <- stats::rnorm(nrow(d)) + 0.6 * (d$A == "a2")
	f <- d[c("A", "B")]
	sp <- .bay_fact_space(c("A", "B", "A:B"))
	nm <- .bay_ft_all(d$y, f, sp, 0.5, 8000L, 1L, t_law = FALSE)
	ln <- vapply(nm, `[[`, numeric(1), "log_ml")
	ref <- .bay_fact_logbf(d$y, f, sp, 0.5, 50000L, 1L)$log_bf[, 1L]
	expect_lt(max(abs((ln - ln[1L]) - ref)), 0.1)
	# Same measure as the normal laws of sys_bay_laws.R (null model in closed form).
	st <- .bay_laws_stats(d$y, interaction(f, drop = TRUE))
	eq0 <- .bay_law_const(st) + .bay_fl_kernel(st, list(), matrix(1, 1L, st$k), matrix(0, 1L, 0L))
	expect_lt(abs(ln[1L] - eq0), 0.05)
})

test_that("several factors: the t law is averaged with the normal laws", {
	set.seed(709001)
	d <- expand.grid(rep = 1:12, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
	d$y <- stats::rt(nrow(d), 3) / sqrt(3) + 0.8 * (d$A == "a2")
	txt <- capture.output(r <- m.bay(y ~ A * B, data = d, draws = 2000, sensitivity = FALSE))
	expect_true("t_eq" %in% r$laws$law)
	expect_equal(sum(r$laws$posterior), 1, tolerance = 1e-8)
	expect_true(any(grepl("Queues|Tails", txt)))
})
