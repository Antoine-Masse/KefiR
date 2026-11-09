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

test_that("changement de strategie dit : groupe constant (integrale divergente), echantillonnage trop pauvre", {
	set.seed(2)
	d <- data.frame(y = c(rep(5, 8), round(rnorm(8, 6, 1), 2), round(rnorm(8, 6.5, 1), 2)), g = factor(rep(c("a", "b", "c"), each = 8)))
	r <- m.bay(y ~ g, data = d, verbose = FALSE, draws = 1000)
	expect_true(any(grepl("constant", r$story$fr)))
	expect_true(any(grepl("Changement de strat.*valeur finie", r$story$fr)))
	expect_false("normal_het" %in% r$laws$law)
	set.seed(1)
	n <- 30
	d3 <- data.frame(y = c(rnorm(n, 0, 1), rnorm(n, 0.5, 1), rt(n, 1.2)), g = factor(rep(c("a", "b", "c"), each = n)))
	r3 <- suppressWarnings(m.bay(y ~ g, data = d3, verbose = FALSE, draws = 100, sensitivity = FALSE))
	if (length(setdiff(c("normal_eq", "normal_het", "t_eq", "t_het"), r3$laws$law))) {
		expect_true(any(grepl("trop pauvre", r3$story$fr)))
		expect_true(any(grepl("Changement de strat", r3$story$fr)))
	} else {
		succeed()
	}
})

test_that("etape 4 : variances inegales annoncees, avec le groupe le plus disperse", {
	set.seed(11)
	d <- data.frame(y = c(rnorm(15, 10, 2), rnorm(15, 11, 2), rnorm(15, 12, 6)), g = factor(rep(c("a", "b", "c"), each = 15)))
	r <- m.bay(y ~ g, data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE)
	expect_true(any(grepl("Variances in.*le groupe c est le plus dispers", r$story$fr)))
	expect_true(any(grepl("Unequal variances.*group c", r$story$en)))
})

test_that("qualificatif : un facteur de Bayes entre 1/3 et 3 donne non tranche par les donnees, verdict indecis", {
	s <- .bay_support_bf(0.94, log(2))
	expect_equal(s$level, 1L)
	expect_match(s$fr, "non tranch")
	expect_equal(.bay_support_bf(0.94, log(20))$level, .bay_support(0.94)$level)
	expect_equal(.bay_support_bf(0.55, log(2))$level, 1L)
	# verdict plafonne : forte probabilite venue de l'a priori
	obj <- list(hypotheses = data.frame(posterior = c(0.9, 0.1), prior = c(0.8, 0.2)), prob_effect = 0.9,
		sensitivity = NULL, contradiction = FALSE)
	expect_match(.bay_verdict(obj)$fr, "ind")
})

test_that("session anglaise reelle : le recit sort en anglais", {
	old <- Sys.getlocale("LC_CTYPE")
	on.exit(Sys.setlocale("LC_CTYPE", old), add = TRUE)
	ok <- suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
	skip_if(identical(ok, "") || !identical(.msg("en", "fr"), "en"), "locale anglaise indisponible")
	txt <- capture.output(m.bay(y ~ g, data = d1, draws = 500, sensitivity = FALSE))
	expect_true(any(grepl("^1\\) Plan recognised", txt)))
	expect_false(any(grepl("Plan reconnu|Contr", txt)))
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
