# Lois d'echantillonnage de m.bay, 2 et k groupes (R/sys_bay_laws.R). Comparaison a
# une quadrature ecrite separement : dev/mbay/scripts/agent_40_lois_calcul_independant.R ;
# ici, des controles plus courts.

sim <- function(n, m = 0, s = 1, seed = 1) {
	set.seed(seed)
	k <- length(n)
	m <- rep_len(m, k)
	s <- rep_len(s, k)
	list(y = unlist(lapply(seq_len(k), function(i) rnorm(n[i], m[i], s[i]))),
		g = factor(rep(seq_len(k), times = n)))
}

test_that("noyau normal : c = 1 et grille en g redonnent zs et anovaBF", {
	skip_if_not_installed("BayesFactor")
	old <- options(BFprogress = FALSE)
	on.exit(options(old), add = TRUE)
	for (des in list(c(15, 15), c(7, 15, 26), rep(10, 5))) {
		d <- sim(des, seq(0, 0.7, length.out = length(des)), seed = sum(des))
		st <- .bay_laws_stats(d$y, d$g)
		one <- matrix(1, 1L, st$k)
		grid <- .bay_law_gint(st, one, c(0.25, 0.5, 1)) - .bay_law_kernel(st, one, 0)
		for (i in 1:3) {
			r <- c(0.25, 0.5, 1)[i]
			expect_lt(abs(grid[1L, i] - .bay_bf_oneway_zs(d$y, d$g, r)), 1e-5)
		}
		ab <- BayesFactor::extractBF(BayesFactor::anovaBF(y ~ g, data.frame(y = d$y, g = d$g),
			rscaleFixed = 0.5, progress = FALSE), logbf = TRUE)$bf
		expect_lt(abs(grid[1L, 2L] - ab), 1e-4)
	}
})

test_that("la loi a variance commune est un cas limite de la loi par groupe (u grand)", {
	d <- sim(c(12, 18, 15), c(0, 0.6, 0), seed = 3)
	st <- .bay_laws_stats(d$y, d$g)
	tab <- .bay_laws_table(st, 0.5, 4000, 4000L, 1L, laws = c("normal_eq", "normal_het"))
	eq <- tab$log_ml[tab$law == "normal_eq"]
	het <- tab$log_ml[tab$law == "normal_het"]
	expect_lt(max(abs(eq - het)), 0.05)
})

test_that("BF et table : invariants par translation et changement d'echelle", {
	d <- sim(c(10, 14, 12), c(0, 0.5, 0.1), c(1, 2, 1), seed = 5)
	ref <- .bay_laws_average(.bay_laws_table(.bay_laws_stats(d$y, d$g), 0.5, 1.5, 1000L, 1L), 0.5)
	autre <- .bay_laws_average(.bay_laws_table(.bay_laws_stats(1000 + 37 * d$y, d$g), 0.5, 1.5, 1000L, 1L), 0.5)
	expect_lt(abs(autre$log_bf - ref$log_bf), 0.05)  # arrondi : memes tirages, recherche du mode a peine differente
	expect_lt(max(abs(as.numeric(autre$p_law) - as.numeric(ref$p_law))), 0.02)
})

test_that("les probabilites sommees font 1 et ne se contredisent pas", {
	d <- sim(c(15, 15), c(0, 0.8), seed = 8)
	a <- .bay_laws_average(.bay_laws_table(.bay_laws_stats(d$y, d$g), 0.5, 1.5, 1000L, 1L), 0.5)
	expect_equal(a$p_h1 + a$p_h0, 1)
	expect_equal(sum(a$p_law), 1)
	expect_equal(sum(a$models$posterior), 1)
	expect_equal(sum(a$models$posterior[a$models$hypothesis == "H1"]), a$p_h1)
})

test_that("etat du RNG utilisateur intact, resultat deterministe a seed fixee", {
	d <- sim(c(10, 10), c(0, 0.3), seed = 9)
	set.seed(42)
	before <- .Random.seed
	r1 <- m.bay(d$y, d$g, verbose = FALSE, draws = 500L, sensitivity = FALSE)
	expect_identical(.Random.seed, before)
	r2 <- m.bay(d$y, d$g, verbose = FALSE, draws = 500L, sensitivity = FALSE)
	expect_identical(r1$prob_effect, r2$prob_effect)
	r3 <- m.bay(d$y, d$g, verbose = FALSE, draws = 500L, sensitivity = FALSE, seed = 2L)
	expect_false(identical(r1$prob_effect, r3$prob_effect))
})

test_that("B1 : un petit groupe qui s'eloigne ne pousse jamais vers H0", {
	lb <- vapply(c(5, 50, 500), function(d) {
		set.seed(7)
		y <- c(rnorm(40), rnorm(3, d, 1))
		g <- rep(1:5, c(10, 10, 10, 10, 3))
		.bay_laws_average(.bay_laws_table(.bay_laws_stats(y, g), 0.5, 1.5, 2000L, 1L), 0.5)$log_bf
	}, numeric(1))
	expect_true(all(diff(lb) > 0))
	expect_gt(lb[3L], 5)
})

test_that("B2 : un groupe disperse sans rapport ne fait pas basculer la reponse", {
	set.seed(3)
	a <- rnorm(20)
	b <- rnorm(20)
	a <- (a - mean(a)) / sd(a)
	b <- (b - mean(b)) / sd(b) + 1.1
	lb <- vapply(c(1, 100, 1e4, 1e6), function(S) {
		set.seed(5)
		y <- c(a, b, rnorm(20, 0, S))
		.bay_laws_average(.bay_laws_table(.bay_laws_stats(y, rep(1:3, each = 20)), 0.5, 1.5, 2000L, 1L), 0.5)$log_bf
	}, numeric(1))
	# Un plancher : stable a grand S et du cote de H1.
	expect_gt(min(lb), 1)
	expect_lt(abs(lb[4L] - lb[3L]), 0.3)
})

test_that("cas limites : refus explicites, jamais de NA silencieux", {
	expect_error(m.bay(c(1, 1, 1, 1), c(1, 1, 2, 2), verbose = FALSE), "binaire|constant")
	expect_error(m.bay(c(.3, .3, .3, .1 + .2, .1 + .2, .1 + .2), rep(1:2, each = 3), verbose = FALSE), "binaire|constant")
	expect_error(m.bay(c(1, 1, 5, 5, 3, 3), c(1, 1, 2, 2, 3, 3), verbose = FALSE), "variation")
	expect_error(m.bay(c(1, 2, 4), c("a", "b", "c"), verbose = FALSE), "libert")
	expect_error(m.bay(c(1, 2.5, 3), c(1, 1, 1), verbose = FALSE), "niveaux")
	r <- m.bay(c(rnorm(10), 3), rep(1:2, c(10, 1)), verbose = FALSE, draws = 500L, sensitivity = FALSE)
	expect_true(is.finite(r$prob_effect))
	expect_true(any(grepl("n < 3", r$notes$en)))
})

test_that("une valeur aberrante n'efface pas l'effet ni ne le cree", {
	set.seed(11)
	y <- c(rnorm(30), rnorm(30, 1.2))
	y[1L] <- 40
	g <- rep(1:2, each = 30)
	r <- m.bay(y, g, verbose = FALSE, draws = 2000L, sensitivity = FALSE)
	# Le gaussien seul perd l'effet ; les queues lourdes le gardent.
	expect_lt(.bay_bf_oneway_zs(y, g, 0.5), 0)
	expect_gt(r$prob_effect, 0.9)
	expect_gt(sum(r$laws$posterior[r$laws$law %in% c("t_eq", "t_het")]), 0.9)
})

# ---- Relecture adverse du 7 novembre -----------------------------------------

test_that("B2 relecture : un groupe constant n'est ni un plantage ni un tirage au sort", {
	y <- c(3, 3, 3, 3, 3, 1, 2, 3, 4, 5)
	g <- rep(1:2, each = 5)
	a <- m.bay(y, g, verbose = FALSE, draws = 500L, sensitivity = FALSE)
	b <- m.bay(y, g, verbose = FALSE, draws = 500L, sensitivity = FALSE, seed = 2L)
	expect_identical(a$prob_effect, b$prob_effect)
	expect_equal(a$laws$law, "normal_eq")
	expect_true(any(grepl("Repeated values", a$notes$en)))
	expect_equal(a$effects$log_bf, .bay_bf_oneway_zs(y, g, 0.5))
	expect_silent(z <- m.bay(c(rep(0, 50), rnorm(10, 1)), rep(1:3, each = 20), verbose = FALSE,
		draws = 300L, sensitivity = FALSE))
	expect_true(is.finite(z$prob_effect))
})

test_that("B3 relecture : le noyau rapide reste juste a c extreme", {
	set.seed(1)
	st <- .bay_laws_stats(rnorm(20), rep(1:2, each = 10))
	for (c2 in c(1e10, 1e14, 1e16, 1e20)) {
		cc <- matrix(c(0.5, c2), 1L, 2L)
		a <- .bay_law_kernel(st, cc, 0.7)
		b <- .bay_law_kernel_grid(st, cc, 0.7)
		expect_lt(abs(a - drop(b)), 1e-6)
	}
})

test_that("I3 relecture : l'integrale en g couvre un rapport inter/intra extreme", {
	set.seed(1)
	y <- c(rnorm(10, 0, 1e-7), rnorm(10, 1, 1e-7), rnorm(10, 0.5, 100))
	g <- rep(1:3, each = 10)
	st <- .bay_laws_stats(y, g)
	cc <- matrix(c(1 / (3 * 0.5), 1 / (3 * 0.5), 1e8), 1L, 3L)
	a <- .bay_law_gint(st, cc, 0.5)
	# reference : meme integrale, grille large et fine, ecrite ici
	ts <- seq(-20, 80, by = 0.05)
	lk <- vapply(ts, function(t) .bay_law_kernel(st, cc, exp(t)), numeric(1))
	lp <- 0.5 * log(0.5^2 / 2) - 0.5 * log(pi) - 0.5 * ts - 0.5^2 / (2 * exp(ts))
	expect_lt(abs(drop(a) - (.bay_log_sum_exp(lk + lp) + log(0.05))), 1e-3)
})

test_that("B4 relecture : entrees mal comprises, refus ou bon plan", {
	d <- data.frame(y = rnorm(30), temps = rep(c("a", "b", "c"), 10), dose = rep(c(1, 2, 5), 10))
	# wt / within sans id : l'ordre des lignes apparie, comme dans m.test (D15)
	expect_identical(m.bay(y ~ temps, data = d, wt = "temps", verbose = FALSE, draws = 500L, sensitivity = FALSE, posthoc = FALSE)$plan, "repeated")
	expect_identical(m.bay(y ~ temps, data = d, within = "temps", verbose = FALSE, draws = 500L, sensitivity = FALSE, posthoc = FALSE)$plan, "repeated")
	expect_error(m.bay(y ~ I(dose > 2), data = d, verbose = FALSE), "variables")
	expect_error(m.bay(y ~ ifelse(dose > 2, "a", "b"), data = d, verbose = FALSE), "variables")
	expect_error(m.bay(d$y, rep(c("a", "b"), length.out = 15), verbose = FALSE), "longueur")
	expect_error(m.bay(y ~ temps + (1 | id / temps), data = transform(d, id = rep(1:10, each = 3)), verbose = FALSE), "seule variable")
	expect_error(m.bay(d$y, as.list(d$temps), verbose = FALSE), "vecteur")
	expect_error(m.bay(y ~ temps, data = d, draws = 1L, verbose = FALSE), "draws")
	expect_error(m.bay(y ~ temps, data = d, seed = NA, verbose = FALSE), "seed")
	expect_error(m.bay(y ~ temps, data = d, rscale = 1e-300, verbose = FALSE), "rscale")
})

test_that("I1 et I2 relecture : loi imposee hors de la plage, cotes exactes", {
	set.seed(11)
	y <- c(rnorm(30), rnorm(30, 1.2))
	y[1L] <- 40
	r <- m.bay(y, rep(1:2, each = 30), verbose = FALSE, draws = 800L)
	expect_lt(unname(r$sensitivity$by_law["normal_eq"]), 0.5)
	expect_false(identical(r$sensitivity$verdict, "changes"))
	r2 <- m.bay(y, rep(1:2, each = 30), verbose = FALSE, draws = 500L, prior_h1 = 0.2)
	expect_equal(r2$sensitivity$odds, colnames(r2$sensitivity$grid))
	expect_match(paste(capture.output(print(r2, detail = TRUE)), collapse = " "), "1:8")
})

test_that("echelles extremes de y : invariance, sans erreur brute", {
	set.seed(4)
	y <- c(rnorm(10), rnorm(10, 1))
	g <- rep(1:2, each = 10)
	go <- function(v) m.bay(v, g, verbose = FALSE, draws = 500L, sensitivity = FALSE)$prob_effect
	expect_lt(abs(go(y * 1e200) - go(y)), 0.01)
	expect_lt(abs(go(y * 1e-200) - go(y)), 0.01)
})

# ---- Deux formes d'effet (RELECTURE 2026-11-07, B1) ----------------------------

test_that("k = 2 : une seule forme, identique a zs et ttestBF", {
	skip_if_not_installed("BayesFactor")
	old <- options(BFprogress = FALSE)
	on.exit(options(old), add = TRUE)
	d <- sim(c(12, 18), c(0, 0.5), seed = 12)
	tab <- .bay_laws_table(.bay_laws_stats(d$y, d$g), 0.5, 1.5, 500L, 1L, laws = "normal_eq")
	lb <- tab$log_ml[tab$hypothesis == "H1"] - tab$log_ml[tab$hypothesis == "H0"]
	tb <- BayesFactor::extractBF(BayesFactor::ttestBF(d$y[d$g == 1], d$y[d$g == 2], rscale = sqrt(2) / 2), logbf = TRUE)$bf
	expect_lt(abs(lb - tb), 1e-6)
	expect_identical(tab$log_ml, tab$log_ml_zs)
})

test_that("la seconde forme est le BF de 's contre le reste', a anovaBF pres", {
	skip_if_not_installed("BayesFactor")
	old <- options(BFprogress = FALSE)
	on.exit(options(old), add = TRUE)
	d <- sim(c(8, 15, 11, 12), c(0, 0.7, 0.1, -0.2), seed = 13)
	st <- .bay_laws_stats(d$y, d$g)
	one <- matrix(1, 1L, 4L)
	for (s in 1:4) {
		ours <- drop(.bay_law_gint(st, one, 0.5 * sqrt(2), s = s)) - .bay_law_kernel(st, one, 0)
		merged <- factor(ifelse(d$g == levels(d$g)[s], "s", "rest"))
		ref <- BayesFactor::extractBF(BayesFactor::anovaBF(y ~ m, data.frame(y = d$y, m = merged),
			rscaleFixed = 0.5, progress = FALSE), logbf = TRUE)$bf
		expect_lt(abs(ours - ref), 1e-4)
	}
	# Le noyau ligne a ligne, avec ou sans s, concorde avec lui-meme.
	tab <- .bay_laws_table(st, 0.5, 1.5, 500L, 1L, laws = "normal_eq")
	h1 <- tab[tab$hypothesis == "H1", ]
	expect_lt(abs(h1$log_ml_zs - (tab$log_ml[tab$hypothesis == "H0"] + .bay_bf_oneway_zs(d$y, d$g, 0.5))), 1e-8)
})

test_that("B1 relecture : k grand, petit groupe eloigne, la preuve croit avec d", {
	for (des in list(c(rep(10, 7), 3), c(rep(5, 11), 2))) {
		lb <- vapply(c(5, 50, 1000), function(d) {
			set.seed(7)
			k <- length(des)
			y <- c(rnorm(sum(des[-k])), rnorm(des[k], d))
			g <- rep(seq_len(k), des)
			.bay_laws_average(.bay_laws_table(.bay_laws_stats(y, g), 0.5, 1.5, 400L, 1L,
				laws = c("normal_eq", "normal_het")), 0.5)$log_bf
		}, numeric(1))
		expect_true(all(diff(lb) > 0))
		expect_gt(lb[3L], 5)
	}
})

test_that("sensibilite : la forme 'tous les groupes' est montree a part", {
	d <- sim(c(10, 10, 10), c(0, 0.8, 0), seed = 14)
	r <- m.bay(d$y, d$g, verbose = FALSE, draws = 400L)
	expect_true(r$sensitivity$all_groups > 0 && r$sensitivity$all_groups < 1)
	expect_match(paste(capture.output(print(r, detail = TRUE)), collapse = " "), "forme")
})

test_that("formes differentes, moyennes egales : le bilan pr\u00e9vient", {
	set.seed(43)
	s <- 0.8
	y <- c((rlnorm(400, 0, s) - exp(s^2 / 2)) / sqrt((exp(s^2) - 1) * exp(s^2)), rnorm(400))
	r <- m.bay(y, rep(1:2, each = 400), verbose = FALSE, draws = 500L, sensitivity = FALSE)
	if (r$prob_effect > 0.75) {
		expect_true(r$contradiction)
		out <- paste(capture.output(print(r)), collapse = " ")
		expect_match(out, "faible|weak")
		expect_match(out, "centres")
	}
})
