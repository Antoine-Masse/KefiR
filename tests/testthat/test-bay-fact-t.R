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

test_that("a yes above 95 % is never shown with an interval covering 0 without a word (audit of a206ba5)", {
	for (s in c(33L, 709001L, 5L)) {
		set.seed(s)
		d <- expand.grid(rep = 1:12, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
		d$y <- stats::rt(nrow(d), 3) / sqrt(3) + 1.2 * (d$A == "a2")
		txt <- capture.output(r <- m.bay(y ~ A * B, data = d, draws = 2000, sensitivity = FALSE))
		fam <- r$posthoc$families
		for (f in fam) {
			if (f$k == 2L && isTRUE(f$p_term > 0.95) && f$pairs$diff_low[1L] < 0 && f$pairs$diff_high[1L] > 0) {
				expect_true(any(grepl("couvre 0|covers 0", txt)))
			} else {
				expect_true(TRUE)
			}
		}
		# The gaps of step 6 come from the mix of laws when the Student law counts.
		if (isTRUE(r$p_t > 0.01)) expect_false(is.null(fam[[1L]]$t_weight))
	}
})
