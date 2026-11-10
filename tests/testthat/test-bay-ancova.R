# m.bay : analyse de covariance (D16).

.anc_data <- function(n = 20, k = 3, eff = c(0, 0, 1), slope_diff = 0, seed = 1, shift = 0) {
	set.seed(seed)
	d <- data.frame(g = factor(rep(letters[seq_len(k)], each = n)), x = runif(n * k, 5, 15) + shift)
	d$y <- 10 + 0.8 * (d$x - shift) + eff[as.integer(d$g)] + slope_diff * (as.integer(d$g) == k) * (d$x - shift - 10) + rnorm(n * k)
	d
}

test_that("ANCOVA : BF du facteur = lmBF(y ~ g + xc) / lmBF(y ~ xc), rscaleCont sqrt(2)/4", {
	d <- .anc_data(seed = 2)
	r <- m.bay(y ~ g + x, data = d, verbose = FALSE, draws = 2000, sensitivity = FALSE)
	expect_identical(r$plan, "ancova")
	d$xc <- d$x - mean(d$x)
	ref <- BayesFactor::lmBF(y ~ g + xc, data = d, rscaleFixed = 0.5, rscaleCont = sqrt(2) / 4, iterations = 50000, progress = FALSE) /
		BayesFactor::lmBF(y ~ xc, data = d, rscaleCont = sqrt(2) / 4, progress = FALSE)
	expect_lt(abs(r$effects$log_bf - BayesFactor::extractBF(ref, logbf = TRUE)$bf), 0.05)
	expect_true(any(grepl("rscaleCont", r$story$detail_fr)))
	expect_true(any(grepl("centr", r$story$fr)))
})

test_that("ANCOVA : centrage, le resultat ne depend pas de l'origine de la covariable", {
	a <- m.bay(y ~ g + x, data = .anc_data(seed = 3), verbose = FALSE, draws = 1000, sensitivity = FALSE, posthoc = FALSE)
	b <- m.bay(y ~ g + x, data = .anc_data(seed = 3, shift = 1000), verbose = FALSE, draws = 1000, sensitivity = FALSE, posthoc = FALSE)
	expect_equal(a$effects$log_bf, b$effects$log_bf, tolerance = 0.02)
})

test_that("ANCOVA : ecarts ajustes proches d'emmeans, sens comme car::Anova type II", {
	d <- .anc_data(n = 60, seed = 4, eff = c(0, 0, 1))
	r <- m.bay(y ~ g + x, data = d, verbose = FALSE, draws = 4000, sensitivity = FALSE)
	fit <- stats::lm(y ~ g + x, data = d)
	em <- as.data.frame(pairs(emmeans::emmeans(fit, ~ g), reverse = TRUE, adjust = "none"))
	# emmeans : "b - a", "c - a", "c - b" ; m.bay : group2 - group1, memes paires
	expect_lt(max(abs(r$posthoc$pairs$diff - em$estimate)), 0.1)
	av <- car::Anova(fit, type = 2)
	expect_lt(av["g", "Pr(>F)"], 0.001)
	expect_gt(r$prob_effect, 0.99)
	expect_true(r$posthoc$pairs$p_differ[r$posthoc$pairs$group1 == "a" & r$posthoc$pairs$group2 == "c"] > 0.95)
})

test_that("ANCOVA : pentes differentes, ecarts aux quartiles comme emmeans(at = ), zone ou les groupes different", {
	d <- .anc_data(n = 40, k = 2, eff = c(0, 0), slope_diff = 0.6, seed = 5)
	r <- m.bay(y ~ g * x, data = d, verbose = FALSE, draws = 4000, sensitivity = FALSE)
	expect_identical(r$posthoc$method, "slopes")
	expect_true(any(grepl("Changement de strat", r$story$fr)))
	fit <- stats::lm(y ~ g * x, data = d)
	q <- stats::quantile(d$x, c(0.25, 0.5, 0.75), names = FALSE)
	em <- as.data.frame(pairs(emmeans::emmeans(fit, ~ g | x, at = list(x = q)), reverse = TRUE))
	expect_lt(max(abs(r$posthoc$gaps$diff - em$estimate)), 0.15)
	# zone : contient au moins une extremite de l'etendue, ou le vrai ecart est le plus grand
	z <- r$posthoc$zones[[1L]]
	expect_gt(nrow(z), 0L)
	expect_true(any(z$from <= min(d$x) + 0.5 | z$to >= max(d$x) - 0.5))
	court <- paste(capture.output(print(r)), collapse = "\n")
	expect_match(gsub("\\s+", " ", court), "dépend de x|depends on x")
})

test_that("ANCOVA : groupes qui se recouvrent peu sur la covariable, prudence et soutien plafonne", {
	set.seed(33)
	d <- data.frame(g = factor(rep(c("R1", "R2"), each = 15)), x = c(runif(15, 20, 35), runif(15, 32, 50)))
	d$y <- 60 + 0.4 * d$x + 4 * (d$g == "R2") + rnorm(30, 0, 1)
	r <- m.bay(y ~ g + x, data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE)
	expect_true(any(r$overlap$poor))
	expect_true(any(grepl("Prudence ajout", r$story$fr)))
	court <- paste(capture.output(print(r)), collapse = "\n")
	expect_match(court, "plafonn|capped")
	expect_false(grepl("soutien (mod|fort|tr)", court))
})

test_that("ANCOVA : plans hors champ refuses clairement", {
	d <- .anc_data(seed = 6)
	expect_error(m.bay(y ~ x, data = d, verbose = FALSE), "covariable seule|covariate alone")
	d$id <- factor(rep(1:20, 3))
	expect_error(m.bay(y ~ g + x + Error(id / g), data = d, verbose = FALSE), "covariables avec des sujets|Covariates with subjects|appari")
})
