# g.bay : pendant bayesien du test G (dev/gbay/DECISIONS.md, G1 a G9).

tab23 <- matrix(c(20, 15, 25, 30, 10, 12), nrow = 2,
	dimnames = list(sexe = c("F", "M"), couleur = c("bleu", "vert", "rouge")))

test_that("independance (G12) : forme close = loi predictive a priori par Monte-Carlo ; meme nombre lignes, colonnes ou total fixes", {
	# Forme jointe (marges neutres) : DM(cases ; a) - DM(lignes ; c a) - DM(colonnes ; r a).
	joint <- function(t, a) .gbay_ldm(as.vector(t), rep(a, length(t))) - .gbay_ldm(rowSums(t), rep(ncol(t) * a, nrow(t))) - .gbay_ldm(colSums(t), rep(nrow(t) * a, ncol(t)))
	for (t in list(tab23, matrix(c(18, 7, 9, 16), 2), matrix(c(30, 10, 12, 8, 25, 9, 5, 6, 28), 3), matrix(c(3, 0, 1, 2, 5, 1), 2))) {
		t <- unname(unclass(t))
		for (a in c(1, 3)) {
			expect_equal(.gbay_lbf_ct(t, a), joint(t, a), tolerance = 1e-10)
			expect_equal(.gbay_lbf_ct(t(t), a), .gbay_lbf_ct(t, a), tolerance = 1e-10)
		}
		r <- g.bay(t, verbose = FALSE, sensitivity = FALSE)
		expect_equal(r$log_bf, .bay_dm_neutral_bf(t, 2), tolerance = 1e-12)
		expect_identical(.bay_dm_neutral_a(), 2)
		expect_equal(r$prob_effect, 1 / (1 + exp(-r$log_bf)), tolerance = 1e-12)
		expect_equal(r$G, unname(suppressWarnings(DescTools::GTest(t)$statistic)), tolerance = 1e-8)
	}
	# Controle independant : vraisemblances marginales par Monte-Carlo sur les a priori.
	set.seed(3)
	t <- matrix(c(4, 1, 2, 5), 2); a <- 3; S <- 200000
	rd <- function(S, al) { g <- matrix(stats::rgamma(S * length(al), rep(al, each = S)), S); g / rowSums(g) }
	p1 <- rd(S, rep(a, 4))
	m1 <- mean(exp(colSums(as.vector(t) * t(log(p1)))))
	pr <- rd(S, rep(2 * a, 2)); pc <- rd(S, rep(2 * a, 2))
	m0 <- mean(exp(rowSums(t)[1] * log(pr[, 1]) + rowSums(t)[2] * log(pr[, 2]) + colSums(t)[1] * log(pc[, 1]) + colSums(t)[2] * log(pc[, 2])))
	expect_equal(log(m1 / m0), .gbay_lbf_ct(t, a), tolerance = 0.02)
	# La table, x et y, et y ~ g donnent le meme nombre.
	set.seed(41)
	d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 30)))
	d$y <- factor(ifelse(rbinom(90, 1, c(.5, .5, .8)[d$g]) == 1, "oui", "non"))
	expect_equal(g.bay(d$g, d$y, verbose = FALSE, sensitivity = FALSE)$log_bf, g.bay(y ~ g, data = d, verbose = FALSE, sensitivity = FALSE)$log_bf)
})

test_that(".bay_dm_neutral_bf (fonction partagee avec m.bay) : forme jointe, symetrie, table de type table", {
	t <- matrix(c(12, 11, 6, 9, 10, 14, 7, 12, 19), 3, byrow = TRUE)
	j <- function(t, a) .gbay_ldm(as.vector(t), rep(a, length(t))) - .gbay_ldm(rowSums(t), rep(ncol(t) * a, nrow(t))) - .gbay_ldm(colSums(t), rep(nrow(t) * a, ncol(t)))
	for (a in c(1, 2, 4)) {
		expect_equal(.bay_dm_neutral_bf(t, a), j(t, a), tolerance = 1e-10)
		expect_equal(.bay_dm_neutral_bf(t(t), a), .bay_dm_neutral_bf(t, a), tolerance = 1e-10)
	}
	expect_equal(.bay_dm_neutral_bf(as.table(t)), .bay_dm_neutral_bf(t))
	expect_identical(.bay_dm_neutral_grid(), c(1, 2, 4))
	# Table de l'audit : BF 0,51 a a = 2 (expert, passe 2).
	expect_equal(exp(.bay_dm_neutral_bf(t, 2)), 0.51, tolerance = 0.01)
})

test_that("calibration : marges equilibrees sans lien (4 x 5, N = 50) ; un lien modere ne fait pas une independance forte", {
	set.seed(2026)
	X <- stats::rmultinom(300, 50, rep(1 / 20, 20))
	lb <- apply(X, 2L, function(x) { t <- matrix(x, 4); .gbay_lbf_ct(t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]) })
	expect_lt(mean(lb > log(19)), 0.01)
	# Table de l'audit : 3 x 3, N = 100, V de Cramer 0,21 ; le BIC disait l'independance a 99,6 %.
	t <- matrix(c(12, 11, 6, 9, 10, 14, 7, 12, 19), 3, byrow = TRUE)
	r <- g.bay(t, verbose = FALSE)
	expect_lt(1 - r$prob_effect, 0.95)
})

test_that("garde-fou de coherence : pas d'independance forte avec des cases qui s'en ecartent", {
	# Table ou l'independance est en tete mais une case s'ecarte (construite pour l'exemple).
	found <- NULL
	set.seed(8)
	for (i in 1:400) {
		t <- matrix(stats::rmultinom(1, 300, rep(1 / 12, 12)), 3)
		if (.gbay_lbf_ct(t) < -log(19)) {
			r <- g.bay(t, verbose = FALSE, sensitivity = FALSE, draws = 2000)
			if (any(r$cells$reading != "-")) { found <- r; break }
		}
	}
	skip_if(is.null(found))
	expect_true(found$guard)
	out <- capture.output(print(found))
	expect_true(any(grepl("ramen\u00e9 \u00e0 mod\u00e9r\u00e9|brought back to moderate", out)))
	expect_false(any(grepl("soutien fort|soutien tr\u00e8s fort|strong support", out)))
})

test_that("visibilite : verbose = TRUE rend un objet invisible, verbose = FALSE un objet visible", {
	invisible(capture.output(v <- withVisible(g.bay(tab23, sensitivity = FALSE, draws = 1000))))
	expect_false(v$visible)
	expect_true(withVisible(g.bay(tab23, verbose = FALSE, sensitivity = FALSE, draws = 1000))$visible)
	expect_length(capture.output(invisible(g.bay(tab23, verbose = FALSE, sensitivity = FALSE, draws = 1000))), 0L)
})

test_that("ajustement : forme close = integration numerique (2 modalites) et multibridge (proportions egales)", {
	for (cs in list(list(n = c(13, 4), p = c(0.6, 0.4), a = 1), list(n = c(7, 2), p = c(0.3, 0.7), a = 0.5), list(n = c(40, 60), p = c(0.5, 0.5), a = 2))) {
		num <- stats::integrate(function(t) stats::dbinom(cs$n[1], sum(cs$n), t) * stats::dbeta(t, cs$a, cs$a), 0, 1)$value
		expect_equal(.gbay_lbf_gof(cs$n, cs$p, cs$a), log(num / stats::dbinom(cs$n[1], sum(cs$n), cs$p[1])), tolerance = 1e-7)
	}
	skip_if_not_installed("multibridge")
	x <- c(3, 7, 12, 5)
	mb <- getExportedValue("multibridge", "mult_bf_equality")(x = x, a = rep(1, 4))
	expect_equal(.gbay_lbf_gof(x, rep(0.25, 4), 1), mb$bf$LogBFe0, tolerance = 1e-7)
})

test_that("formes d'entree equivalentes : table, x et y, ~ A + B, n ~ A + B, noms de colonnes", {
	set.seed(3)
	d <- data.frame(A = sample(c("a", "b", "c"), 80, TRUE), B = sample(c("u", "v"), 80, TRUE))
	ref <- g.bay(table(d$A, d$B), verbose = FALSE, sensitivity = FALSE)$log_bf
	expect_equal(g.bay(d$A, d$B, verbose = FALSE, sensitivity = FALSE)$log_bf, ref)
	expect_equal(g.bay(~ A + B, data = d, verbose = FALSE, sensitivity = FALSE)$log_bf, ref)
	expect_equal(g.bay("A", "B", data = d, verbose = FALSE, sensitivity = FALSE)$log_bf, ref)
	k <- as.data.frame(table(A = d$A, B = d$B))
	expect_equal(g.bay(Freq ~ A + B, data = k, verbose = FALSE, sensitivity = FALSE)$log_bf, ref)
	expect_identical(g.bay(~ A + B, data = d, verbose = FALSE, sensitivity = FALSE)$names, c("A", "B"))
	# Manquants retires et dits.
	d$A[2] <- NA
	r <- g.bay(d$A, d$B, verbose = FALSE, sensitivity = FALSE)
	expect_equal(r$n, 79)
	expect_true(any(grepl("1 observation", r$notes$en)))
})

test_that("ajustement : vecteur, facteur, ~ A, matrice a une ligne ; p, rescale.p ; refus", {
	n <- c(a = 30, b = 12, c = 8)
	ref <- g.bay(n, p = c(0.5, 0.3, 0.2), verbose = FALSE, sensitivity = FALSE)
	expect_identical(ref$test, "gof")
	expect_equal(ref$log_bf, .gbay_lbf_gof(n, c(0.5, 0.3, 0.2), 1))
	f <- factor(rep(names(n), n))
	expect_equal(g.bay(f, p = c(0.5, 0.3, 0.2), verbose = FALSE, sensitivity = FALSE)$log_bf, ref$log_bf)
	expect_equal(g.bay(~ f, p = c(0.5, 0.3, 0.2), verbose = FALSE, sensitivity = FALSE)$log_bf, ref$log_bf)
	r1 <- g.bay(matrix(n, 1, dimnames = list(NULL, names(n))), p = c(5, 3, 2), rescale.p = TRUE, verbose = FALSE, sensitivity = FALSE)
	expect_equal(r1$log_bf, ref$log_bf)
	expect_true(any(grepl("Change of strategy", r1$notes$en)))
	# Proportions egales par defaut, comme GTest.
	expect_equal(unname(g.bay(n, verbose = FALSE, sensitivity = FALSE)$p), rep(1 / 3, 3))
	expect_error(g.bay(n, p = c(0.5, 0.6, 0.1), verbose = FALSE))
	expect_error(g.bay(n, p = c(0.5, 0.5, 0), verbose = FALSE))
	expect_error(g.bay(c(1.5, 2), verbose = FALSE))
	expect_error(g.bay(array(1:8, c(2, 2, 2)), verbose = FALSE))
	# Une proportion nulle, jamais observee : modalite retiree.
	expect_equal(length(g.bay(c(10, 5, 0), p = c(0.5, 0.5, 0), verbose = FALSE, sensitivity = FALSE)$table), 2L)
})

test_that("ligne ou colonne vide retiree et dite ; refus clairs", {
	t0 <- cbind(tab23, noir = 0)
	r <- g.bay(t0, verbose = FALSE, sensitivity = FALSE)
	expect_equal(dim(r$table), c(2L, 3L))
	expect_equal(r$log_bf, g.bay(tab23, verbose = FALSE, sensitivity = FALSE)$log_bf)
	expect_true(any(grepl("noir", r$notes$fr)))
	expect_error(g.bay(y ~ A * B, data = data.frame(y = 1:4, A = 1:4, B = 1:4), verbose = FALSE))
})

test_that("2 x 2 : rapport des cotes proche de glm a grand effectif ; V de Cramer proche de l'observe", {
	set.seed(7)
	n <- 2000
	d <- data.frame(g = factor(rep(c("a", "b"), each = n)))
	d$y <- factor(rbinom(2 * n, 1, c(.3, .45)[d$g]), labels = c("non", "oui"))
	r <- g.bay(d$g, d$y, verbose = FALSE)
	lor <- unname(stats::coef(stats::glm(y ~ g, data = d, family = stats::binomial()))[2])
	expect_lt(abs(log(r$odds_ratio[["median"]]) - lor), 0.03)
	expect_gt(r$odds_ratio[["p_gt1"]], 0.999)
	t <- table(d$g, d$y)
	v_obs <- sqrt(suppressWarnings(stats::chisq.test(t, correct = FALSE))$statistic / sum(t))
	expect_lt(abs(r$cramer_v[["median"]] - v_obs), 0.01)
	expect_true(r$cramer_v[["low"]] < v_obs && v_obs < r$cramer_v[["high"]])
	# Cases : les quatre lues, de cote oppose deux a deux.
	expect_setequal(r$cells$reading, c("more", "less"))
})

test_that("cases : rapport observe / attendu retrouve a grand effectif", {
	set.seed(11)
	P <- matrix(c(.20, .10, .05, .05, .15, .45), 2)
	x <- stats::rmultinom(1, 5000, as.numeric(P))
	t <- matrix(x, 2, dimnames = list(c("r1", "r2"), c("a", "b", "c")))
	r <- g.bay(t, verbose = FALSE, sensitivity = FALSE)
	oe <- as.numeric(t / (outer(rowSums(t), colSums(t)) / sum(t)))
	expect_lt(max(abs(r$cells$ratio - oe)), 0.03)
	expect_true(all(r$cells$ratio_low < r$cells$ratio & r$cells$ratio < r$cells$ratio_high))
})

test_that("return = FALSE : un nombre, egal a prob_effect ; graine sans effet sur l'etat global", {
	set.seed(99); u <- stats::runif(1)
	set.seed(99)
	p <- g.bay(tab23, return = FALSE)
	expect_equal(stats::runif(1), u)
	expect_length(p, 1L)
	expect_equal(p, g.bay(tab23, verbose = FALSE)$prob_effect)
	expect_equal(g.bay(tab23, verbose = FALSE)$cells, g.bay(tab23, verbose = FALSE)$cells)
})

test_that("recit : 80 colonnes au plus, marques sur chaque ligne, francais et anglais, rejeu a l'identique", {
	set.seed(5)
	d <- data.frame(tr = factor(rep(c("temoin", "froid", "chaud"), each = 30)))
	d$germe <- factor(ifelse(rbinom(90, 1, c(.5, .5, .8)[d$tr]) == 1, "oui", "non"))
	big <- matrix(c(30, 10, 12, 8, 25, 9, 5, 6, 28), 3, dimnames = list(sol = c("argile", "sable", "limon"), espece = c("chene", "pin", "frene")))
	cas <- list(function() g.bay(tab23), function() g.bay(c(705, 224), p = c(3, 1) / 4, correct = "williams"),
		function() g.bay(matrix(c(18, 7, 9, 16), 2)), function() g.bay(germe ~ tr, data = d), function() g.bay(big),
		function() g.bay(c(rouge = 12, vert = 3, bleu = 1, jaune = 0)), function() g.bay(big, detail = TRUE))
	vis <- function(x) nchar(gsub("\t", "        ", x, fixed = TRUE), type = "width")
	old_loc <- Sys.getlocale("LC_CTYPE"); old_w <- getOption("width")
	on.exit({ Sys.setlocale("LC_CTYPE", old_loc); options(width = old_w) }, add = TRUE)
	ok_line <- function(s) s == "" || grepl("^([0-9]+\\) |   |\t(==> |--> |Note|    |      |       )|Method|M\u00e9thode|Bilan)", s)
	for (lc in c("French_France.utf8", "English_United States.utf8")) {
		if (identical(suppressWarnings(Sys.setlocale("LC_CTYPE", lc)), "")) next
		for (w in c(80, 120)) {
			options(width = w)
			for (f in cas) {
				out <- capture.output(f())
				expect_lte(max(vis(out)), 80)
				expect_true(all(vapply(out, ok_line, logical(1))))
				expect_false(any(grepl("p-value =|p = 0", out)))
			}
		}
		options(width = 60)
		expect_lte(max(vis(capture.output(cas[[5]]()))), 60)
		options(width = 80)
		r <- g.bay(big, verbose = FALSE)
		told <- capture.output(g.bay(big))
		expect_identical(capture.output(print(r, story = TRUE)), told)
	}
	expect_true(any(grepl("khi-deux|chi-squared", c(g.bay(big, verbose = FALSE)$story$fr, g.bay(c(3, 1, 0, 2), verbose = FALSE)$story$en))))
	Sys.setlocale("LC_CTYPE", "English_United States.utf8")
	expect_error(g.bay(c(1.5, 2), verbose = FALSE), "whole numbers")
})

