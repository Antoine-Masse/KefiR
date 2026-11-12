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
	# D31: deviance Bayes factor (Zellner-Siow g, r = sqrt(2)/4) from the deviance gained by glm.
	z <- glm(y ~ 1, data = d, family = binomial)$deviance - glm(y ~ tr, data = d, family = binomial)$deviance
	lb_tr <- .bay_dev_lbf(z, 1, nrow(d), sqrt(2) / 4)
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
	z <- 2 * (as.numeric(logLik(MASS::polr(o ~ recette, data = d))) - as.numeric(logLik(MASS::polr(o ~ 1, data = d))))
	lb <- .bay_dev_lbf(z, 1, nrow(d), sqrt(2) / 4)
	expect_equal(r$models$log_bf[r$models$model == "recette"], lb, tolerance = 1e-6)
	expect_gt(r$effects$p_h1[1L], 0.9)
	expect_gt(r$ratios$ratio[1L], 1)
})

test_that("deviance Bayes factor with a Zellner-Siow g (D31): same integral by brute force, sane limits", {
	brute <- function(z, p, N, r) {
		b <- N * r^2 / 2
		f <- function(g) (1 + g)^(-p / 2) * exp(g / (1 + g) * z / 2) * sqrt(b) / gamma(0.5) * g^(-1.5) * exp(-b / g)
		log(stats::integrate(f, 0, Inf, rel.tol = 1e-10)$value)
	}
	for (cs in list(c(0, 1, 100), c(4, 1, 100), c(12, 2, 400), c(30, 3, 80))) {
		expect_equal(.bay_dev_lbf(cs[1], cs[2], cs[3], sqrt(2) / 4), brute(cs[1], cs[2], cs[3], sqrt(2) / 4), tolerance = 1e-5)
	}
	expect_lt(.bay_dev_lbf(0, 1, 100), 0)
	expect_gt(.bay_dev_lbf(30, 1, 100), .bay_dev_lbf(10, 1, 100))
	expect_equal(.bay_dev_lbf(5, 0, 100), 0)
})
