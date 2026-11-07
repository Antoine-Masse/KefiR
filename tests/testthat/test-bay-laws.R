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
	expect_equal(autre$log_bf, ref$log_bf, tolerance = 1e-6)
	expect_equal(as.numeric(autre$p_law), as.numeric(ref$p_law), tolerance = 1e-6)
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
