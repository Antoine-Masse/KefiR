# Internal heteroscedastic ANOVA (Brunner-Dette-Munk), which replaced fda.usc::fanova.hetero().
# Reference values: fda.usc 2.2.0, fanova.hetero(data.frame(x, g), x ~ g)$ans.

test_that(".bdm_oneway reproduces the fda.usc reference values", {
	x <- c(5, 1, 7, 3, 9, 2, 8, 4, 6, 2)
	g <- c("b", "a", "c", "b", "c", "a", "c", "b", "a", "c")
	res <- KefiR:::.bdm_oneway(x, g)
	expect_equal(unname(res), c(2.715667, 1.70996, 5.533132, 0.1524928), tolerance = 1e-6)
	expect_named(res, c("statistic", "df1", "df2", "p.value"))
})

test_that(".bdm_oneway ignores missing pairs and the order of the data", {
	x <- c(5, 1, 7, 3, 9, 2, 8, 4, 6, 2)
	g <- c("b", "a", "c", "b", "c", "a", "c", "b", "a", "c")
	ref <- KefiR:::.bdm_oneway(x, g)
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
