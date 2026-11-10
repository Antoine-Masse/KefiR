# Repeated measures with replicates (D21): subject and subject x condition random.

.mixed_data <- function(seed, eff = c(0, 0.5, 0.6), n_id = 10L, n_rep = 4L) {
	set.seed(seed)
	d <- expand.grid(rep = seq_len(n_rep), cond = factor(c("c1", "c2", "c3")), id = factor(seq_len(n_id)))
	u0 <- stats::rnorm(n_id)
	uc <- matrix(stats::rnorm(3L * n_id, 0, 0.7), n_id)
	d$y <- u0[d$id] + uc[cbind(as.integer(d$id), as.integer(d$cond))] + eff[d$cond] + stats::rnorm(nrow(d))
	d
}

test_that("replicates: same Bayes factor as lmBF with subject and subject x condition", {
	skip_if_not_installed("BayesFactor")
	d <- .mixed_data(5)
	r <- m.bay(y ~ cond | id, data = d, verbose = FALSE, sensitivity = FALSE, draws = 2000)
	expect_identical(r$plan, "mixed")
	set.seed(1)
	ref <- BayesFactor::extractBF(BayesFactor::lmBF(y ~ cond + id + cond:id, data = d, whichRandom = "id", iterations = 50000, progress = FALSE) /
		BayesFactor::lmBF(y ~ id + cond:id, data = d, whichRandom = "id", iterations = 50000, progress = FALSE), logbf = TRUE)$bf
	expect_lt(abs(r$effects$log_bf[1L] - ref), 0.2)
	# The replicates are not counted as subjects: far below the model with the subject only.
	naive <- BayesFactor::extractBF(BayesFactor::lmBF(y ~ cond + id, data = d, whichRandom = "id", progress = FALSE) /
		BayesFactor::lmBF(y ~ id, data = d, whichRandom = "id", progress = FALSE), logbf = TRUE)$bf
	expect_lt(r$effects$log_bf[1L], naive - 2)
	expect_identical(r$posthoc$k, 3L)
})

test_that("replicates with a factor between subjects, story and post-hoc", {
	d <- .mixed_data(6)
	d$g <- factor(ifelse(as.integer(d$id) <= 5L, "A", "B"))
	txt <- capture.output(r <- m.bay(y ~ g * cond + Error(id / cond), data = d, draws = 1000))
	expect_identical(r$plan, "mixed")
	expect_true(any(grepl("avec r\u00e9pliques|with replicates", txt)))
	expect_true(any(grepl("R\u00e9pliques|Replicates", txt)))
	expect_identical(r$posthoc$method, "factorial")
})
