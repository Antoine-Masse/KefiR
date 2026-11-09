# m.bay : recit en etapes (verbose = TRUE), conclusion sans repetition, objet complet.

d1 <- local({
	set.seed(1)
	data.frame(y = c(rnorm(12, 10), rnorm(12, 11), rexp(12) + 10),
		g = factor(rep(c("a", "b", "c"), each = 12)))
})

test_that("verbose = FALSE : rien n'est imprime, le recit est garde dans l'objet", {
	out <- capture.output(r <- m.bay(y ~ g, data = d1, verbose = FALSE, draws = 1000))
	expect_length(out, 0L)
	expect_true(length(r$story$fr) > 0L)
	expect_equal(length(r$story$en), length(r$story$fr))
	titres <- grep("^[0-9]+\\) ", r$story$fr, value = TRUE)
	expect_equal(as.integer(sub("\\).*", "", titres)), seq_along(titres))
	expect_equal(length(titres), 6L)
	expect_true(any(grepl("==>", r$story$fr, fixed = TRUE)))
	expect_true(any(grepl("-->", r$story$fr, fixed = TRUE)))
})

test_that("verbose = TRUE : six etapes, puis la conclusion sans titre ni point de depart repetes", {
	txt <- capture.output(r <- m.bay(y ~ g, data = d1, draws = 1000))
	expect_true(any(grepl("^1\\) ", txt)))
	expect_true(any(grepl("^6\\) ", txt)))
	expect_false(any(grepl("Bilan de m.bay|Bilan of m.bay", txt)))
	expect_false(any(grepl("^\tPoint de d|^\tStarting point", txt)))
	# chaque hypothese est argumentee
	expect_equal(sum(grepl("^\t {5}[^ ]", txt)), 2L)
	# jamais de p-value chiffree
	expect_false(any(grepl("p[ -]?(value)? ?[=<>] ?[0-9]", txt)))
	# print(story = TRUE) rejoue le meme recit
	again <- capture.output(print(r, story = TRUE))
	expect_identical(again, txt)
	# print() seul : le bilan complet, avec titre et point de depart
	court <- capture.output(print(r))
	expect_true(any(grepl("Bilan de m.bay|Bilan of m.bay", court)))
	expect_true(any(grepl("Point de d|Starting point", court)))
})

test_that("changement de strategie dit : ex-aequo, lois ecartees", {
	set.seed(4)
	d <- data.frame(y = c(rep(1:3, 6), rep(2:4, 6)), g = factor(rep(c("a", "b"), each = 18)))
	r <- suppressWarnings(m.bay(y ~ g, data = d, verbose = FALSE, draws = 1000))
	if (nrow(r$laws) < 4L) {
		expect_true(any(grepl("Changement de strat|Change of strategy", c(r$story$fr, r$story$en))))
	} else {
		succeed()
	}
	expect_true(any(grepl("Ex-aequo", r$story$fr)))
})

test_that("plan refuse : l'etape 1 le nomme avant le refus", {
	d <- data.frame(y = rnorm(20), g = factor(rep(1:2, 10)), x = rnorm(20))
	txt <- capture.output(expect_error(m.bay(y ~ g + x, data = d)))
	expect_true(any(grepl("^1\\) ", txt)))
	expect_true(any(grepl("ANCOVA", txt)))
})

test_that("plusieurs facteurs : loi normale seule dite, a priori HOP, lmBF cite, sensibilite racontee", {
	set.seed(3)
	d <- expand.grid(A = factor(1:2), B = factor(1:3), rep = 1:5)
	d$y <- 1.2 * (d$A == "2") + rnorm(nrow(d))
	txt <- capture.output(r <- m.bay(y ~ A * B, data = d))
	all <- paste(txt, collapse = "\n")
	expect_match(all, "lmBF() {BayesFactor}", fixed = TRUE)
	expect_match(all, "HOP")
	expect_match(all, "variance commune|common variance")
	expect_match(all, "Sensibilit|Sensitivity")
	expect_equal(sum(grepl("^[0-9]+\\) ", txt)), 6L)
	expect_false(any(grepl("Bilan de m.bay|Bilan of m.bay", txt)))
	# un argument par terme et un pour la reponse globale
	expect_equal(sum(grepl("^\t {5}[^ ]", txt)), nrow(r$effects) + 1L)
	expect_identical(capture.output(print(r, story = TRUE)), txt)
})
