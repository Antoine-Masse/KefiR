# m.bay : recit en etapes (verbose = TRUE), bilan visible, objet complet, mise en forme
# commune (relecture/lisibilite_2026-11-09).

d1 <- local({
	set.seed(1)
	data.frame(y = c(rnorm(12, 10), rnorm(12, 11), rexp(12) + 10),
		g = factor(rep(c("a", "b", "c"), each = 12)))
})

test_that("verbose = FALSE : rien n'est imprime, le recit est garde dans l'objet", {
	out <- capture.output(r <- m.bay(y ~ g, data = d1, verbose = FALSE, draws = 1000))
	expect_length(out, 0L)
	expect_true(length(r$story$fr) > 0L)
	# lignes repliees : leur nombre depend de la langue
	titres <- grep("^[0-9]+\\) ", r$story$fr, value = TRUE)
	expect_equal(as.integer(sub("\\).*", "", titres)), seq_along(titres))
	# six etapes racontees ; le bilan (7) est imprime a part
	expect_equal(length(titres), 6L)
	expect_equal(r$story$k, 6L)
	expect_true(any(grepl("==>", r$story$fr, fixed = TRUE)))
	expect_true(any(grepl("-->", r$story$fr, fixed = TRUE)))
})

test_that("verbose = TRUE : six etapes, bilan visible, une marque par ligne, rejouable", {
	txt <- capture.output(r <- m.bay(y ~ g, data = d1, draws = 1000))
	expect_true(any(grepl("^1\\) ", txt)))
	i7 <- grep("^7\\) Bilan", txt)
	expect_length(i7, 1L)
	expect_identical(txt[i7 - 1L], "")
	expect_false(any(grepl("Bilan de m.bay|Bilan of m.bay", txt)))
	expect_true(any(grepl("^\t==> 1\\. ", txt)) && any(grepl("^\t==> 2\\. ", txt)))
	# chaque ligne porte une marque : titre, ==>, -->, Note, ou suite d'une ligne repliee
	body <- txt[nzchar(txt) & !grepl("^[0-9]+\\) |^   ", txt)]
	expect_true(all(grepl("^\t(==> |--> |Note : |Note: | {4,})", body)))
	# jamais de p-value chiffree
	expect_false(any(grepl("p[ -]?(value)? ?[=<>] ?[0-9]", txt)))
	# print(story = TRUE) rejoue le meme recit
	again <- capture.output(print(r, story = TRUE))
	expect_identical(again, txt)
	# print() seul : le bilan, avec son titre
	court <- capture.output(print(r))
	expect_true(any(grepl("Bilan de m.bay|Bilan of m.bay", court)))
	expect_true(any(grepl("Verdict", court)))
	# le detail garde la methode (a priori, references)
	det <- paste(capture.output(print(r, detail = TRUE)), collapse = "\n")
	expect_match(det, "Rouder et al. 2012")
})

test_that("changement de strategie dit : groupe constant (integrale divergente), echantillonnage trop pauvre", {
	set.seed(2)
	d <- data.frame(y = c(rep(5, 8), round(rnorm(8, 6, 1), 2), round(rnorm(8, 6.5, 1), 2)), g = factor(rep(c("a", "b", "c"), each = 8)))
	r <- m.bay(y ~ g, data = d, verbose = FALSE, draws = 1000)
	fr <- paste(r$story$fr, collapse = " ")
	expect_match(fr, "constant")
	expect_match(gsub("\\s+", " ", fr), "Changement de stratégie.*valeur finie")
	expect_false("normal_het" %in% r$laws$law)
	set.seed(1)
	n <- 30
	d3 <- data.frame(y = c(rnorm(n, 0, 1), rnorm(n, 0.5, 1), rt(n, 1.2)), g = factor(rep(c("a", "b", "c"), each = n)))
	r3 <- suppressWarnings(m.bay(y ~ g, data = d3, verbose = FALSE, draws = 100, sensitivity = FALSE))
	if (length(setdiff(c("normal_eq", "normal_het", "t_eq", "t_het"), r3$laws$law))) {
		fr3 <- paste(r3$story$fr, collapse = " ")
		expect_match(fr3, "Trop peu de tirages")
		expect_match(fr3, "Changement de strat")
	} else {
		succeed()
	}
})

test_that("etape 3 : variances inegales annoncees, avec le groupe le plus disperse", {
	set.seed(11)
	d <- data.frame(y = c(rnorm(15, 10, 2), rnorm(15, 11, 2), rnorm(15, 12, 6)), g = factor(rep(c("a", "b", "c"), each = 15)))
	r <- m.bay(y ~ g, data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE)
	fr <- gsub("\\s+", " ", paste(r$story$fr, collapse = " ")); en <- gsub("\\s+", " ", paste(r$story$en, collapse = " "))
	expect_match(fr, "Variances inégales probables.*c le plus dispers")
	expect_match(en, "Unequal variances probable.*c the most dispersed")
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
	# ANCOVA acceptee depuis D16 : un plan encore refuse, la MANOVA
	d <- data.frame(y = rnorm(20), y2 = rnorm(20), g = factor(rep(1:2, 10)))
	txt <- capture.output(expect_error(m.bay(cbind(y, y2) ~ g, data = d)))
	expect_true(any(grepl("^1[)] ", txt)))
	expect_true(any(grepl("MANOVA", txt)))
})

test_that("plusieurs facteurs : loi normale seule dite, meme squelette, un terme par ligne, jargon au detail", {
	set.seed(3)
	d <- expand.grid(A = factor(1:2), B = factor(1:3), rep = 1:5)
	d$y <- 1.2 * (d$A == "2") + rnorm(nrow(d))
	txt <- capture.output(r <- m.bay(y ~ A * B, data = d))
	all <- paste(txt, collapse = "\n")
	expect_match(all, "lmBF() {BayesFactor}", fixed = TRUE)
	expect_match(all, "variance commune|common variance")
	expect_match(all, "Solidit")
	# meme squelette : 6 etapes racontees, dont l'etape 6 sans paires, puis le bilan
	expect_equal(sum(grepl("^[0-9]+\\) ", txt)), 7L)
	expect_true(any(grepl("^6\\) .*(non disponible|not available)", txt)))
	expect_false(any(grepl("Bilan de m.bay|Bilan of m.bay", txt)))
	expect_equal(sum(grepl("^\t==> (A|B|A:B) \\(", txt)), nrow(r$effects))
	det <- paste(capture.output(print(r, detail = TRUE)), collapse = "\n")
	expect_match(det, "HOP")
	expect_identical(capture.output(print(r, story = TRUE)), txt)
})

test_that("80 colonnes au plus, tabulation comptee 8 : 7 plans, francais et anglais, et width = 60", {
	set.seed(42)
	d2 <- data.frame(y = c(rnorm(20, 10, 2), rnorm(20, 11.5, 2)), g = factor(rep(c("A", "B"), each = 20)))
	d4 <- data.frame(y = c(rnorm(15, 10, 2), rnorm(15, 12, 2), rnorm(15, 10.5, 2), rnorm(15, 13, 2)), g = factor(rep(c("a", "b", "c", "d"), each = 15)))
	d22 <- expand.grid(rep = 1:10, A = c("A1", "A2"), B = c("B1", "B2", "B3"))
	d22$y <- 10 + 1.5 * (d22$A == "A2") + 1 * (d22$B == "B3") + rnorm(nrow(d22), 0, 2)
	d22$A <- factor(d22$A); d22$B <- factor(d22$B)
	ns <- 15; ind <- rnorm(ns, 0, 2)
	dp <- data.frame(id = factor(rep(1:ns, 2)), t = factor(rep(c("avant", "apres"), each = ns), levels = c("avant", "apres")))
	dp$y <- 20 + ind[as.integer(dp$id)] + 1.2 * (dp$t == "apres") + rnorm(2 * ns, 0, 1.2)
	dr <- data.frame(id = factor(rep(1:ns, 3)), t = factor(rep(c("T1", "T2", "T3"), each = ns)))
	dr$y <- 20 + ind[as.integer(dr$id)] + c(0, 1, 1.8)[as.integer(dr$t)] + rnorm(3 * ns, 0, 1.2)
	da <- data.frame(g = factor(rep(c("t1", "t2", "t3"), each = 15)), x = runif(45, 0, 10))
	da$y <- 5 + 0.8 * da$x + c(0, 1.5, 3)[as.integer(da$g)] + rnorm(45, 0, 1.5)
	dj <- data.frame(tr = factor(rep(c("temoin", "froid", "chaud"), each = 30)))
	dj$germe <- factor(ifelse(rbinom(90, 1, c(.5, .5, .8)[dj$tr]) == 1, "oui", "non"))
	cas <- list(function() m.bay(y ~ g, data = d2, draws = 1000), function() m.bay(y ~ g, data = d4, draws = 1000),
		function() m.bay(y ~ A * B, data = d22), function() m.bay(y ~ t, data = dp, id = "id", paired = TRUE),
		function() m.bay(y ~ t | id, data = dr), function() m.bay(y ~ g + x, data = da), function() m.bay(germe ~ tr, data = dj))
	vis <- function(x) nchar(gsub("\t", "        ", x, fixed = TRUE), type = "width")
	old_loc <- Sys.getlocale("LC_CTYPE"); old_w <- getOption("width")
	on.exit({ Sys.setlocale("LC_CTYPE", old_loc); options(width = old_w) }, add = TRUE)
	langs <- c("French_France.utf8", "English_United States.utf8")
	for (lc in langs) {
		if (identical(suppressWarnings(Sys.setlocale("LC_CTYPE", lc)), "")) next
		for (w in c(80, 120)) {
			options(width = w)
			for (f in cas) expect_lte(max(vis(capture.output(f()))), 80)
		}
	}
	options(width = 60)
	expect_lte(max(vis(capture.output(cas[[2]]()))), 60)
	expect_lte(max(vis(capture.output(cas[[6]]()))), 60)
})
