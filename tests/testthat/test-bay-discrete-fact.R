# Discrete responses with several factors (D25): logistic regression and counts, BIC, HOP.

test_that("yes/no with two factors: BIC of the glm models, HOP, answers per term", {
	set.seed(2)
	d <- expand.grid(rep = 1:25, tr = factor(c("temoin", "engrais")), sol = factor(c("argile", "sable")))
	d$germe <- factor(ifelse(stats::rbinom(nrow(d), 1, ifelse(d$tr == "engrais", 0.7, 0.4)) == 1, "oui", "non"), levels = c("non", "oui"))
	txt <- capture.output(r <- m.bay(germe ~ tr * sol, data = d, sensitivity = FALSE))
	expect_identical(r$plan, "discrete_fact")
	expect_identical(r$kind, "binary")
	# Independent check: -BIC/2 of the model with tr against the null, by hand.
	d$y <- as.integer(d$germe == "oui")
	lb_tr <- (-BIC(glm(y ~ tr, data = d, family = binomial)) + BIC(glm(y ~ 1, data = d, family = binomial))) / 2
	expect_equal(r$models$log_bf[r$models$model == "tr"], lb_tr, tolerance = 1e-8)
	expect_gt(r$effects$p_h1[1L], 0.75)
	expect_true(any(grepl("rapport des cotes|odds ratio", txt)))
	expect_true(any(grepl("^7[)] Bilan - tr oui", txt)))
})

test_that("counts with two factors: Poisson and negative binomial averaged", {
	set.seed(3)
	d <- expand.grid(rep = 1:25, tr = factor(c("a", "b")), sol = factor(c("argile", "sable")))
	d$n <- stats::rnbinom(nrow(d), mu = ifelse(d$sol == "sable", 6, 3), size = 2)
	r <- m.bay(n ~ tr * sol, data = d, verbose = FALSE, sensitivity = FALSE, response = "counts")
	expect_identical(r$plan, "discrete_fact")
	expect_identical(r$kind, "counts")
	expect_gt(r$laws$posterior[r$laws$law == "negbin"], 0.5)
	expect_gt(r$effects$p_h1[2L], 0.75)
	expect_true(all(c("ratio", "low", "high") %in% names(r$ratios)))
})

test_that("ordered response with two factors (D26): polr per model, BIC, HOP", {
	set.seed(4)
	d <- expand.grid(rep = 1:20, recette = factor(c("r1", "r2")), sexe = factor(c("F", "H")))
	d$note <- pmin(5, pmax(1, round(3 + 0.9 * (d$recette == "r2") + stats::rnorm(nrow(d), 0, 1))))
	r <- m.bay(note ~ recette * sexe, data = d, verbose = FALSE, sensitivity = FALSE)
	expect_identical(r$plan, "discrete_fact")
	expect_identical(r$kind, "ordinal")
	d$o <- factor(d$note, ordered = TRUE)
	lb <- (-BIC(MASS::polr(o ~ recette, data = d)) + BIC(MASS::polr(o ~ 1, data = d))) / 2
	expect_equal(r$models$log_bf[r$models$model == "recette"], lb, tolerance = 1e-6)
	expect_gt(r$effects$p_h1[1L], 0.9)
	expect_gt(r$ratios$ratio[1L], 1)
})
