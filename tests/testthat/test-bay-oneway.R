# Facteurs de Bayes a un facteur : variance commune (zs) et variances par groupe.
# Les messages de .exit() sortent en francais quelle que soit la locale (voir
# R/sys_exit.R) ; ceux de warning(.msg()) suivent la locale : motifs a deux langues.

sim_oneway <- function(n, eff, sds = 1, seed = 1) {
	set.seed(seed)
	k <- length(n)
	g <- factor(rep(seq_len(k), times = n))
	y <- rnorm(sum(n), mean = rep(eff, times = n), sd = rep(rep_len(sds, k), times = n))
	list(y = y, g = g)
}

stats_of <- function(y, g) .bay_hetero_prep(y, g)

# ---- zs contre BayesFactor --------------------------------------------------

test_that("zs : anovaBF, equilibre et desequilibre, k = 2, 3, 5, r = 0.5 et 1", {
	skip_if_not_installed("BayesFactor")
	old <- options(BFprogress = FALSE)
	on.exit(options(old), add = TRUE)
	designs <- list(c(15, 15), c(8, 21), c(12, 12, 12), c(7, 15, 26),
		rep(10, 5), c(5, 9, 14, 20, 6))
	for (i in seq_along(designs)) {
		n <- designs[[i]]
		for (r in c(0.5, 1)) {
			d <- sim_oneway(n, seq(0, 0.7, length.out = length(n)), seed = 100 + i)
			ab <- BayesFactor::extractBF(BayesFactor::anovaBF(y ~ g,
				data.frame(y = d$y, g = d$g), rscaleFixed = r, progress = FALSE),
				logbf = TRUE)
			expect_lt(abs(.bay_bf_oneway_zs(d$y, d$g, r) - ab$bf), 1e-5)
		}
	}
})

test_that("zs : ttestBF(rscale = sqrt(2)/2) a k = 2", {
	skip_if_not_installed("BayesFactor")
	for (n in list(c(20, 20), c(11, 30))) {
		d <- sim_oneway(n, c(0, 0.6), seed = sum(n))
		x1 <- d$y[d$g == 1]
		x2 <- d$y[d$g == 2]
		tb <- BayesFactor::extractBF(BayesFactor::ttestBF(x1, x2, rscale = sqrt(2) / 2),
			logbf = TRUE)
		expect_lt(abs(.bay_bf_oneway_zs(d$y, d$g, 1/2) - tb$bf), 1e-5)
	}
})

# ---- zs : stabilite ---------------------------------------------------------

# Groupes aux residus centres : F est fixe exactement, a l'arrondi pres.
data_with_F <- function(Fval, n = 10, k = 3, seed = 3) {
	set.seed(seed)
	e <- rnorm(n * k)
	g <- factor(rep(seq_len(k), each = n))
	e <- e - ave(e, g)
	e <- e / sqrt(sum(e^2) / (n * k - k))
	mus <- seq(-1, 1, length.out = k)
	mus <- mus * sqrt(Fval * (k - 1) / (n * sum(mus^2)))
	list(y = mus[g] + e, g = g)
}

test_that("zs : F enorme, valeur finie et croissante, jamais NA", {
	Fs <- c(1e2, 1e4, 1e6, 1e8, 1e10, 1e12, 1e14)
	vals <- vapply(Fs, function(Fv) {
		d <- data_with_F(Fv)
		.bay_bf_oneway_zs(d$y, d$g)
	}, numeric(1))
	expect_false(anyNA(vals))
	expect_true(all(is.finite(vals)))
	expect_true(all(diff(vals) > 0))
})

test_that("zs : pente asymptotique en log F, valeur theorique (N - 1 - k)/2", {
	# Pour F >> 1 : (SST/q(g))^{(N-1)/2} apporte F^{(N-1)/2} ; l'integrale en t de
	# exp(-k t/2) jusqu'au pic, situe en g ~ F, retire F^{k/2}. N = 30, k = 3 : 13.
	d1 <- data_with_F(1e10)
	d2 <- data_with_F(1e12)
	slope <- (.bay_bf_oneway_zs(d2$y, d2$g) - .bay_bf_oneway_zs(d1$y, d1$g)) / log(100)
	expect_equal(slope, (30 - 1 - 3) / 2, tolerance = 1e-6)
})

test_that("zs : la ou l'integrale de BayesFactor echoue, nous restons finis", {
	skip_if_not_installed("BayesFactor")
	old <- options(BFprogress = FALSE)
	on.exit(options(old), add = TRUE)
	# anovaBF a F = 1e4 : ancrage, la route est encore commune.
	bf_of <- function(d) {
		out <- NULL
		utils::capture.output(out <- suppressMessages(suppressWarnings(tryCatch(
			BayesFactor::extractBF(BayesFactor::anovaBF(y ~ g, data.frame(y = d$y, g = d$g),
				progress = FALSE), logbf = TRUE)$bf, error = function(e) NA_real_))),
			type = "message")
		out
	}
	d4 <- data_with_F(1e4)
	expect_lt(abs(.bay_bf_oneway_zs(d4$y, d4$g) - bf_of(d4)), 1e-5)
	d6 <- data_with_F(1e6)
	ab6 <- bf_of(d6)
	# BayesFactor 0.9.12 : NA a F = 1e6 (integrate, "roundoff error"). Une version
	# corrigee doit alors concorder avec nous, pas casser ce test.
	v6 <- .bay_bf_oneway_zs(d6$y, d6$g)
	expect_true(is.finite(v6))
	if (!is.na(ab6)) expect_lt(abs(v6 - ab6), 1e-5)
})

test_that("zs : residu nul, Inf avec signal ; y constant et N = k refuses", {
	g <- factor(rep(1:3, each = 4))
	y <- c(1, 2, 4)[g]
	expect_warning(res <- .bay_bf_oneway_zs(y, g), "zero|nulle")
	expect_identical(as.numeric(res), Inf)
	expect_identical(attr(res, "signal"), "zero_residual")
	expect_error(.bay_bf_oneway_zs(rep(2, 12), g), "constant")
	expect_error(.bay_bf_oneway_zs(c(1, 2, 3), factor(1:3)), "N <= k")
})

test_that("zs : F = 0 exact, valeur finie et negative", {
	g <- factor(rep(1:3, each = 5))
	y <- rep(c(-2, -1, 0, 1, 2), 3)
	v <- .bay_bf_oneway_zs(y, g)
	expect_true(is.finite(v))
	expect_lt(v, 0)
})

# ---- hetero : constantes ----------------------------------------------------

test_that("hetero : H0 par IS contre integrale exacte, a moins de 3 erreurs MC", {
	cases <- list(
		list(n = c(30, 30, 30), eff = c(0, 0.2, 0.4), sds = c(1, 1, 1)),
		list(n = c(8, 20, 45), eff = c(0, 0.3, 0), sds = c(0.5, 1, 3)),
		list(n = c(5, 5), eff = c(0, 0.5), sds = c(1, 2)),
		list(n = c(2, 3, 4, 6), eff = c(0, 0, 0.1, 0), sds = 1))
	for (i in seq_along(cases)) {
		cs <- cases[[i]]
		d <- sim_oneway(cs$n, cs$eff, cs$sds, seed = 40 + i)
		st <- stats_of(d$y, d$g)
		exact <- .bay_hetero_h0_exact(st)
		is <- .bay_hetero_h0_is(st, n_draws = 200000L, seed = 2L)
		expect_lt(abs(exact - is$log_ml), 3 * is$rel_se)
	}
})

test_that("hetero : H0 exacte contre integrate() sur un cas regulier", {
	d <- sim_oneway(c(12, 15, 9), c(0, 0.4, 0.2), c(1, 1.5, 0.8), seed = 8)
	st <- stats_of(d$y, d$g)
	f <- function(mu) {
		vapply(mu, function(u) {
			exp(sum(st$n / 2 * (log(2) - log(st$ss + st$n * (u - st$m)^2))))
		}, numeric(1))
	}
	val <- integrate(f, -3, 3, subdivisions = 500L, rel.tol = 1e-10)$value
	ref <- log(val) + sum(-st$n / 2 * log(2 * pi) + lgamma(st$n / 2))
	expect_equal(.bay_hetero_h0_exact(st), ref, tolerance = 1e-7)
})

test_that("hetero : H0 exacte stable quand les groupes sont tres separes", {
	# Plusieurs pics, ou un pic etroit loin de toute moyenne.
	g <- factor(rep(1:3, each = 40))
	set.seed(5)
	y <- rnorm(120, sd = 0.01) + c(0, 1e3, 5e3)[g]
	expect_true(is.finite(.bay_hetero_h0_exact(stats_of(y, g))))
})

test_that("hetero : BF -> 1 quand r -> 0 (H1 se referme sur H0, meme constante)", {
	d <- sim_oneway(c(20, 25, 15), c(0, 0.3, 0.1), c(1, 1.4, 0.8), seed = 21)
	st <- stats_of(d$y, d$g)
	res <- .bay_hetero_core(st, r = 1e-4, n_draws = 5000L, seed = 1L)
	expect_lt(abs(res$log_bf), 1e-2)
})

test_that("hetero : k = 2, H1 contre integration directe en delta a sigma fixe", {
	# Controle du g integre : a variances fixees, delta ~ Cauchy(0, sqrt(2)*r) rend
	# la marginale de m1 - m2 calculable par integrate(), sans passer par g.
	s2 <- matrix(c(1.2, 0.7), 1L, 2L)
	st <- list(k = 2L, n = c(10, 14), m = c(0.5, -0.3), ss = c(8, 9))
	inner <- .bay_hetero_log_h1_inner(s2, st, r = 0.5, h = 0.125)
	v <- s2[1, ] / st$n
	tau <- sqrt(mean(s2))
	diff_m <- st$m[1] - st$m[2]
	f <- function(delta) {
		dnorm(diff_m, delta * tau, sqrt(sum(v))) * dcauchy(delta, 0, sqrt(2) / 2)
	}
	val <- integrate(f, -Inf, Inf, rel.tol = 1e-10)$value
	# Le noyau apres mu plat vaut la densite de m1 - m2, sans facteur additionnel.
	expect_equal(inner, log(val), tolerance = 1e-7)
})

# ---- hetero : invariances, convergence, RNG ---------------------------------

test_that("hetero : invariance par translation et par changement d'echelle", {
	d <- sim_oneway(c(15, 22, 30), c(0, 0.5, 0.2), c(1, 2, 0.7), seed = 11)
	base <- .bay_bf_oneway_hetero(d$y, d$g, n_draws = 5000L)
	expect_equal(.bay_bf_oneway_hetero(d$y + 1000, d$g, n_draws = 5000L)$log_bf,
		base$log_bf, tolerance = 1e-6)
	expect_equal(.bay_bf_oneway_hetero(d$y * 1e4, d$g, n_draws = 5000L)$log_bf,
		base$log_bf, tolerance = 1e-6)
	expect_equal(.bay_bf_oneway_hetero(d$y * 1e-6, d$g, n_draws = 5000L)$log_bf,
		base$log_bf, tolerance = 1e-6)
})

test_that("hetero : invariance d'echelle reelle, sans la standardisation interne", {
	d <- sim_oneway(c(15, 22, 30), c(0, 0.5, 0.2), c(1, 2, 0.7), seed = 12)
	st <- stats_of(d$y, d$g)
	cc <- 37
	st2 <- st
	st2$m <- st$m * cc
	st2$ss <- st$ss * cc^2
	a <- .bay_hetero_core(st, 0.5, 5000L, 1L)
	b <- .bay_hetero_core(st2, 0.5, 5000L, 1L)
	expect_equal(b$log_bf, a$log_bf, tolerance = 1e-6)
	# Chaque marginale change de c^{-(N-1)} : mu plat rend un facteur c.
	expect_equal(b$log_ml_h0 - a$log_ml_h0, -(st$N - 1) * log(cc), tolerance = 1e-6)
	expect_equal(b$log_ml_h1 - a$log_ml_h1, -(st$N - 1) * log(cc), tolerance = 1e-6)
})

test_that("hetero : convergence quand n_draws croit", {
	d <- sim_oneway(c(12, 18, 25), c(0, 0.6, 0.1), c(1, 2, 0.6), seed = 13)
	r1 <- .bay_bf_oneway_hetero(d$y, d$g, n_draws = 2000L, seed = 1L)
	r2 <- .bay_bf_oneway_hetero(d$y, d$g, n_draws = 8000L, seed = 2L)
	r3 <- .bay_bf_oneway_hetero(d$y, d$g, n_draws = 64000L, seed = 3L)
	expect_lt(r3$error_pct, r1$error_pct)
	expect_lt(abs(r1$log_bf - r3$log_bf), 3 * (r1$error_pct + r3$error_pct) / 100)
	expect_lt(abs(r2$log_bf - r3$log_bf), 3 * (r2$error_pct + r3$error_pct) / 100)
})

test_that("hetero : le pas de grille n'influence pas le resultat", {
	d <- sim_oneway(c(12, 18, 25), c(0, 3, 0.1), c(1, 2, 0.6), seed = 14)
	st <- stats_of(d$y, d$g)
	s2 <- withr::with_seed(9, .bay_hetero_draws(st, 500L))
	a <- .bay_hetero_log_h1_inner(s2, st, 0.5, 0.25)
	b <- .bay_hetero_log_h1_inner(s2, st, 0.5, 0.03125)
	expect_lt(max(abs(a - b)), 1e-6)
})

test_that("hetero : deterministe a seed fixee, sensible a la seed", {
	d <- sim_oneway(c(10, 12), c(0, 0.5), c(1, 2), seed = 15)
	a <- .bay_bf_oneway_hetero(d$y, d$g, n_draws = 3000L, seed = 7L)
	b <- .bay_bf_oneway_hetero(d$y, d$g, n_draws = 3000L, seed = 7L)
	cc <- .bay_bf_oneway_hetero(d$y, d$g, n_draws = 3000L, seed = 8L)
	expect_identical(a$log_bf, b$log_bf)
	expect_false(identical(a$log_bf, cc$log_bf))
})

test_that("hetero : etat du RNG utilisateur intact, independant de RNGkind", {
	d <- sim_oneway(c(10, 12), c(0, 0.5), c(1, 2), seed = 16)
	set.seed(123)
	before <- get(".Random.seed", envir = globalenv())
	a <- .bay_bf_oneway_hetero(d$y, d$g, n_draws = 1000L)
	expect_identical(get(".Random.seed", envir = globalenv()), before)
	kinds <- RNGkind()
	withr::with_seed(5, .rng_kind = "L'Ecuyer-CMRG", {
		before2 <- get(".Random.seed", envir = globalenv())
		b <- .bay_bf_oneway_hetero(d$y, d$g, n_draws = 1000L)
		expect_identical(get(".Random.seed", envir = globalenv()), before2)
		expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")
	})
	expect_identical(RNGkind(), kinds)
	expect_identical(a$log_bf, b$log_bf)
})

test_that("hetero : ESS sous le seuil, deux augmentations puis signal", {
	d <- sim_oneway(c(10, 12), c(0, 0.5), c(1, 2), seed = 17)
	st <- stats_of(d$y, d$g)
	expect_warning(res <- .bay_hetero_core(st, 0.5, 500L, 1L, min_ess = Inf), "effective")
	expect_identical(res$n_draws, 8000L)
	expect_identical(res$n_grow, 2L)
	expect_true(res$low_ess)
})

# ---- hetero : modele voisin -------------------------------------------------

test_that("hetero : variances voisines et grands groupes, proche de zs (pas identique)", {
	# Ecart documente : les modeles different (k variances contre une, tau^2 aleatoire).
	# Seuil 0.5 sur le log BF, a comparer a l'erreur MC, tres inferieure.
	for (s in 1:3) {
		d <- sim_oneway(c(300, 300, 300), c(0, 0.2, 0.1), c(1, 1.02, 0.98), seed = 60 + s)
		zs <- .bay_bf_oneway_zs(d$y, d$g)
		he <- .bay_bf_oneway_hetero(d$y, d$g, n_draws = 5000L)
		expect_lt(abs(zs - he$log_bf), 0.5)
	}
})

test_that("hetero : k = 2 ne se retourne pas contre l'evidence", {
	d <- sim_oneway(c(40, 40), c(0, 2), c(1, 3), seed = 71)
	expect_gt(.bay_bf_oneway_hetero(d$y, d$g, n_draws = 3000L)$log_bf, 2)
	d0 <- sim_oneway(c(100, 100), c(0, 0), c(1, 3), seed = 72)
	expect_lt(.bay_bf_oneway_hetero(d0$y, d0$g, n_draws = 3000L)$log_bf, 0)
})

test_that("hetero : F tres grand, sortie finie", {
	g <- factor(rep(1:3, each = 20))
	set.seed(73)
	y <- rnorm(60, sd = c(0.001, 0.002, 0.0015)[g]) + c(0, 10, 25)[g]
	res <- .bay_bf_oneway_hetero(y, g, n_draws = 5000L)
	expect_true(is.finite(res$log_bf))
	expect_true(res$log_bf > 50)
})

test_that("hetero : sorties annoncees, sans p-value", {
	d <- sim_oneway(c(10, 12, 14), c(0, 0.5, 0), c(1, 2, 1), seed = 74)
	res <- .bay_bf_oneway_hetero(d$y, d$g, r = 1, n_draws = 2000L)
	expect_true(all(c("log_bf", "error_pct", "ess", "n_draws", "rscale", "tau_def",
		"method") %in% names(res)))
	expect_identical(res$rscale, 1)
	expect_identical(res$tau_def, "mean_var")
	expect_identical(res$method, "hetero_is")
	expect_false(any(grepl("^p", names(res)) | grepl("value", names(res))))
})

# ---- refus ------------------------------------------------------------------

test_that("refus : k < 2, n_i < 2, SS_i = 0 ; NA annonces et filtres", {
	y <- c(1.1, 2.3, 3.2, 4.8, 5.1, 6.9)
	expect_error(.bay_bf_oneway_hetero(y, rep("a", 6)), "k < 2")
	expect_error(.bay_bf_oneway_zs(y, rep("a", 6)), "k < 2")
	expect_error(.bay_bf_oneway_hetero(y, c("a", "b", "b", "b", "b", "b")), "n_i < 2")
	expect_error(.bay_bf_oneway_hetero(c(2, 2, 2, 1, 3, 5), rep(c("a", "b"), each = 3)),
		"SS_i = 0")
	# Le filtrage peut faire tomber un groupe sous n_i = 2.
	expect_error(suppressWarnings(.bay_bf_oneway_hetero(c(1, NA, 3, 4, 7, 2),
		rep(c("a", "b", "c"), each = 2))), "n_i < 2")
	expect_error(.bay_bf_oneway_hetero(c(1, 2, Inf, 4), c(1, 1, 2, 2)), "infini")
	expect_error(.bay_bf_oneway_hetero(y, rep(1:2, 3), r = -1), "r ")
	expect_error(.bay_bf_oneway_hetero(y, rep(1:2, 3), n_draws = 10), "n_draws")
	d <- sim_oneway(c(10, 12), c(0, 0.5), c(1, 2), seed = 80)
	y_na <- c(d$y, NA, 1)
	g_na <- c(as.character(d$g), "1", NA)
	expect_warning(a <- .bay_bf_oneway_hetero(y_na, g_na, n_draws = 1000L), "NA")
	b <- .bay_bf_oneway_hetero(d$y, d$g, n_draws = 1000L)
	expect_identical(a$log_bf, b$log_bf)
	expect_warning(a2 <- .bay_bf_oneway_zs(y_na, g_na), "NA")
	expect_equal(a2, .bay_bf_oneway_zs(d$y, d$g))
})
