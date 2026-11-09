# m.bay : entrees, retour et bilan.

set.seed(21)
d <- data.frame(y = c(rnorm(12, 10, 2), rnorm(12, 12, 2), rnorm(12, 10.5, 2)),
	variete = rep(c("a", "b", "c"), each = 12), stringsAsFactors = FALSE)
go <- function(...) m.bay(..., verbose = FALSE, draws = 500L, sensitivity = FALSE)

test_that("toutes les formes d'entree donnent le meme resultat", {
	ref <- go(y ~ variete, data = d)$prob_effect
	e <- d
	expect_identical(go(formula = y ~ variete, data = d)$prob_effect, ref)
	expect_identical(go(y ~ variete, d)$prob_effect, ref)
	expect_identical(go(d$y, d$variete)$prob_effect, ref)
	expect_identical(go("y", "variete", data = d)$prob_effect, ref)
	expect_identical(go(1, 2, data = d)$prob_effect, ref)
	expect_identical(go(y ~ factor(variete), data = d)$prob_effect, ref)
	expect_identical(go(e$y ~ e$variete)$prob_effect, ref)
	expect_identical(go(d$y, data.frame(v = d$variete))$prob_effect, ref)
	y <- d$y
	variete <- d$variete
	expect_identical(go(y ~ variete)$prob_effect, ref)
})

test_that("log(y) n'est transforme qu'une fois (BP-033)", {
	a <- go(log(y) ~ variete, data = d)
	e <- d
	e$ly <- log(e$y)
	expect_identical(a$prob_effect, go(ly ~ variete, data = e)$prob_effect)
})

test_that("les NA sont retires et annonces", {
	e <- d
	e$y[c(1, 20)] <- NA
	e$variete[30] <- NA
	r <- go(y ~ variete, data = e)
	expect_equal(r$n, nrow(e) - 3L)
	expect_true(any(grepl("3 observation", r$notes$en)))
})

test_that("structure de retour : sans p-value, a priori annonces", {
	r <- go(y ~ variete, data = d)
	expect_s3_class(r, "mbay")
	expect_true(all(c("term", "bf10", "log_bf", "p_h1", "p_h0", "prior_h1", "rscale",
		"error_pct", "analogue") %in% names(r$effects)))
	expect_false(any(grepl("p.value|pvalue", names(r$effects))))
	expect_false(any(grepl("p.value|pvalue", names(r))))
	expect_equal(sum(r$hypotheses$posterior), 1)
	expect_true(r$hypotheses$posterior[1L] >= r$hypotheses$posterior[2L])
	expect_equal(nrow(r$laws), 4L)
	expect_equal(r$priors$prior_h1, 1/2)
	expect_equal(r$priors$rscale, 1/2)
	expect_equal(r$descriptive_stats$n, c(12, 12, 12))
})

test_that("return = FALSE donne un nombre, sans impression", {
	expect_silent(p <- m.bay(y ~ variete, data = d, return = FALSE, draws = 500L, sensitivity = FALSE))
	expect_length(p, 1L)
	expect_true(p > 0 && p < 1)
})

test_that("prior_h1 change la probabilite selon la formule attendue, pas le BF", {
	a <- go(y ~ variete, data = d)
	b <- go(y ~ variete, data = d, prior_h1 = 0.2)
	expect_equal(a$effects$log_bf, b$effects$log_bf)
	expect_equal(b$prob_effect, plogis(a$effects$log_bf + qlogis(0.2)))
})

test_that("arguments de m.test sans objet : acceptes et signales", {
	r <- go(y ~ variete, data = d, alpha = 0.01, boot = FALSE)
	expect_true(any(grepl("alpha, boot", r$notes$en)))
})

test_that("plans pas encore valides : refus explicite", {
	e <- transform(d, x = rnorm(36), y2 = rnorm(36), id = rep(1:12, 3), B = rep(c("u", "v"), 18),
		b = rep(0:1, 18))
	expect_error(go(y ~ variete + x, data = e), "covariables")
	expect_error(go(cbind(y, y2) ~ variete, data = e), "MANOVA")
	expect_error(go(y ~ variete + Error(id/variete), data = e), "mesures r")
	expect_error(go(y ~ variete + (1 | id), data = e), "mesures r")
	expect_error(go(y ~ variete | id, data = e), "mesures r")
	expect_error(go(y ~ B, data = e, id = "id", paired = TRUE), "appari|mesures|mixte")
	expect_error(go(b ~ variete, data = e), "binaire")
	expect_error(go(y ~ absent, data = e), "introuvable")
	expect_error(go(y ~ 1, data = e), "Aucun facteur")
	expect_error(go(y ~ variete, data = e, prior_h1 = 0), "prior_h1")
	expect_error(go(y ~ variete, data = e, rscale = -1), "r doit")
	expect_error(go(y ~ variete, data = e, laws = "cauchy"), "laws")
})

test_that("bilan : annonce a priori et sensibilite, ne laisse pas croire a une p-value", {
	r <- m.bay(y ~ variete, data = d, verbose = FALSE, draws = 500L)
	out <- capture.output(print(r, detail = TRUE))
	txt <- paste(out, collapse = "\n")
	expect_match(txt, "r = ")
	expect_match(txt, "Dirichlet")
	expect_match(txt, "p-value")
	expect_true(any(grepl("priors|a priori", out)))
	expect_false(grepl("confiance|confidence", txt))
	# Un seul niveau de sensibilite annonce.
	expect_true(r$sensitivity$verdict %in% c("changes", "crosses", "stable"))
	expect_equal(dim(r$sensitivity$grid), c(3L, 3L))
	r0 <- go(y ~ variete, data = d)
	expect_match(paste(capture.output(print(r0)), collapse = "\n"), "sensitivity = FALSE")
	# Par defaut, court ; le detail est sur demande.
	court <- capture.output(print(r))
	# 16 lignes avant le post-hoc ; les paires en ajoutent k(k-1)/2 + 2 (D11).
	expect_lte(length(court), 16L + nrow(r$posthoc$pairs) + 2L)
	expect_false(any(grepl("Dirichlet", court)))
	expect_match(paste(court, collapse = " "), "p-value")
})

test_that("verdict de sensibilite : hypothese en tete, seuil, stable", {
	expect_equal(.bay_sensitivity_verdict(0.8, c(0.4, 0.9)), "changes")
	expect_equal(.bay_sensitivity_verdict(0.8, c(0.7, 0.9)), "crosses")
	expect_equal(.bay_sensitivity_verdict(0.97, c(0.96, 0.99)), "stable")
})

test_that("sources ASCII", {
	src <- c("m.bay.R", "sys_bay_input.R", "sys_bay_bilan.R", "sys_bay_laws.R", "sys_bay_oneway.R", "sys_bay_factorial.R", "sys_bay_story.R", "sys_bay_posthoc.R")
	root <- testthat::test_path("..", "..", "R")
	skip_if_not(dir.exists(root))
	for (f in src) {
		bytes <- readBin(file.path(root, f), "raw", file.info(file.path(root, f))$size)
		expect_true(all(as.integer(bytes) < 128L), info = f)
	}
})
